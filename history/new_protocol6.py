#!/usr/bin/env python3
"""
Tiny-TSS-PV ▸ benchmark with separate Trace / TrVer / R-box counters
────────────────────────────────────────────────────────────────────
"""

import secrets, hashlib, time, random, curve25519_python as cp

# ────────── normalise to bytes ──────────
E    = lambda v: v if isinstance(v, (bytes, bytearray)) else bytes(v)

# ────────── global counters ──────────
CNT     = 0                   # all group muls outside the R-box
RB_CNT  = 0                   # muls consumed *inside* the reconstruction box
bump    = lambda: globals().__setitem__("CNT", CNT + 1)
rb_bump = lambda d: globals().__setitem__("RB_CNT", RB_CNT + d)
reset   = lambda: globals().__setitem__("CNT", 0)
rb_reset= lambda: globals().__setitem__("RB_CNT", 0)

# ────────── curve wrappers ──────────
B   = b'\x58' + b'\x66'*31
P   = lambda A,B: E(cp.point_addition(A,B))
Pm  = lambda s,Pt: bump() or E(cp.point_multiply(s,Pt))
Sm  = lambda s:    bump() or E(cp.scalar_multiply(s))
S2  = lambda a,b:  E(cp.scalar_multiply_scalar(a,b))
Inv = lambda z:    E(cp.scalar_inverse(z))
Sub = lambda a,b:  E(cp.scalar_subtraction(a,b))
rand= lambda: secrets.token_bytes(32)

ONE = b'\x01' + b'\0'*31
H1, H2 = Pm(b'\x02'+b'\0'*31, B), Pm(b'\x03'+b'\0'*31, B)

# ────────── polynomial & sharing helpers ──────────
gpoly = lambda t: [Sm(rand()) for _ in range(t)]

def hp(x,g):
    acc,xp = ONE,x
    for gk in g:
        acc = P(acc, Pm(xp,gk));  xp = S2(xp,x)
    return acc                      # (k-1) muls

cm  = lambda r,s: P(Pm(r,H1), Pm(s,H2))
rho = lambda h,r,s: P(Pm(r,h ), Pm(s,B ))

def proof(h,ρ,cm_,r,s):
    kr,ks = rand(),rand()
    Acm,Ar= P(Pm(kr,H1),Pm(ks,H2)), P(Pm(kr,h),Pm(ks,B))
    c     = hashlib.sha256(b''.join([cm_,ρ,Acm,Ar])).digest()[:32]
    zr,zs = Sub(kr,S2(c,r)), Sub(ks,S2(c,s))
    return c,zr,zs                 # 6 muls

def check(h,ρ,cm_,π):
    c,zr,zs = π
    Mc = P(Pm(c,cm_), P(Pm(zr,H1), Pm(zs,H2)))
    Mr = P(Pm(c,ρ ),  P(Pm(zr,h ), Pm(zs,B )))
    return hashlib.sha256(b''.join([cm_,ρ,Mc,Mr])).digest()[:32]==c

def recon(shs):                    # naive O(k²) Lagrange
    acc=ONE
    for j,(xj,yj) in enumerate(shs):
        l=ONE
        for m,(xm,_) in enumerate(shs):
            if m==j: continue
            l=S2(l,S2(xm,Inv(Sub(xm,xj))))
        acc=P(acc,Pm(l,yj))
    return acc                     # exactly k muls

# ────────── verification helpers ──────────
ShD = lambda T: int(all(check(*line) for line in T))

def ShS(sh,T,g):
    x,y=sh; fx=hp(x,g)
    return int(any(fx==f and y==ρ for f,ρ,*_ in T))

_sample = lambda x,g,r,s: (x, rho(hp(x,g),r,s))

# ────────── reconstruction oracle with own counter ──────────
def Rbox(k, embeds):
    def R(extra):
        extra=[extra] if isinstance(extra,tuple) else list(extra)
        uniq={}; [uniq.setdefault(x,y) for x,y in embeds+extra]
        if len(uniq)<k: return None
        before,t0 = CNT, time.perf_counter()
        res = recon(list(uniq.items())[:k])
        rb_bump(CNT-before)
        globals().__setitem__("RB_TIME",
            globals().get("RB_TIME",0.0) + (time.perf_counter()-t0))
        return res
    return R

# ────────── Trace / TrVer (exclude R-box muls) ──────────
def Trace(tk,T,g,f,t,R):
    ζ,shv=tk; r,s,S=ζ; I=set(range(len(shv)))
    boxX={x for x,_ in (R([]) or [])}
    for ℓ,(xl,_) in enumerate(shv):
        ban, fresh = boxX|{xl}, []
        while len(fresh)<t-f-1:
            x=rand()
            if x not in ban: ban.add(x); fresh.append(_sample(x,g,r,s))
        if R(fresh+[shv[ℓ]])==S: I.discard(ℓ)
    return sorted(I), [shv[i] for i in I]

def TrVer(vk,I,T,π,g,f,t,R):
    ζ,shv=vk; r,s,S=ζ
    if any(_sample(x,g,r,s)!=(x,ρ) for x,ρ in shv): return 0
    for idx in I:
        ban,fresh={x for x,_ in (R([]) or [])}|{shv[idx][0]},[]
        while len(fresh)<t-f-1:
            x=rand()
            if x not in ban: ban.add(x); fresh.append(_sample(x,g,r,s))
        if R(fresh+[shv[idx]])==S and ShS(shv[idx],T,g): return 0
    return 1

# ────────── benchmark harness ──────────
def run(n,k,f):
    assert 0<f<k-1
    print(f"\n== n={n} k={k} f={f} ==")
    g=gpoly(k-1); xs=[(i+1).to_bytes(32,'little') for i in range(n)]

    reset(); t0=time.perf_counter()
    hs=[hp(x,g) for x in xs]
    print(f"[share ] {1e3*(time.perf_counter()-t0)/n:6.1f} ms  {CNT//n:3d} mul/ply")

    r,s=rand(),rand(); cm_=cm(r,s)
    reset(); t0=time.perf_counter()
    ρs,T=[],[]
    for h in hs:
        ρ=rho(h,r,s); π=proof(h,ρ,cm_,r,s)
        ρs.append(ρ); T.append((h,ρ,cm_,π))
    print(f"[dealer] {1e3*(time.perf_counter()-t0):6.1f} ms  {CNT:4d} mul")

    reset(); t0=time.perf_counter(); ok=ShD(T)
    print(f"[ShD   ] {1e3*(time.perf_counter()-t0):6.1f} ms  {CNT//n:3d} mul  ok={ok}")

    reset(); t0=time.perf_counter()
    allok=all(ShS((x,ρ),T,g) for x,ρ in zip(xs[:k],ρs[:k]))
    print(f"[ShS   ] {1e3*(time.perf_counter()-t0):6.1f} ms  {CNT//k:3d} mul  all={allok}")

    rb_reset(); globals().__setitem__("RB_TIME",0.0)
    R=Rbox(k, random.sample(list(zip(xs,ρs)), f))

    reset(); t0=time.perf_counter()
    ok = recon(list(zip(xs,ρs))[:k]) == Sm(s)
    print(f"[Recon ] {1e3*(time.perf_counter()-t0):6.1f} ms  {CNT:4d} mul  ok={ok}")

    tk=((r,s,Sm(s)), list(zip(xs,ρs)))

    # -------- Trace (exclude oracle) --------
    rb_before,c_before,t0=RB_CNT,CNT,time.perf_counter()
    I,π=Trace(tk,T,g,f,k,R)
    dt  = time.perf_counter()-t0
    trace_mul = CNT - c_before - (RB_CNT - rb_before)
    trace_ms  = 1e3*(dt - globals()["RB_TIME"] )
    print(f"[Trace ] {trace_ms:6.1f} ms  {trace_mul:4d} mul  |I|={len(I)}")

    # -------- TrVer (exclude oracle) --------
    rb_before,c_before,t0=RB_CNT,CNT,time.perf_counter()
    ok=TrVer(tk,I,T,π,g,f,k,R)
    dt  = time.perf_counter()-t0
    trv_mul  = CNT - c_before - (RB_CNT - rb_before)
    trv_ms   = 1e3*(dt - (globals()["RB_TIME"] - globals().get("RB_TIME_USED",0.0)))
    globals()["RB_TIME_USED"]=globals()["RB_TIME"]
    print(f"[TrVer ] {trv_ms:6.1f} ms  {trv_mul:4d} mul  ok={ok}")

    # -------- R-box totals --------
    print(f"[R-box ] {1e3*globals()['RB_TIME']:6.1f} ms  {RB_CNT:4d} mul  (all queries)")

# ────────── entry ──────────
if __name__=="__main__":
    for n,k,f in [(32,22,4),(64,44,10),(128,88,20)]:
        run(n,k,f)
