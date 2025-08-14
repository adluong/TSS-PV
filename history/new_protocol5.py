#!/usr/bin/env python3
"""
TSS-PV (no-mask) – micro-benchmark with 3-phase evaluation
==========================================================

Phase-1  (Sharing)        : players compute f(x_i); dealer builds ρ, cm, π;
                            ShDVer runs.
Phase-2  (Reconstruction) : each player submits (x_i,ρ_i); ShSVer runs;
                            k shares reconstruct S.
Phase-3  (Tracing)        : Trace run on the transcript; TrVer checks result.

Multiplication counters:
  • players (share) : k-1 muls
  • dealer          : 6·n muls   (2 for ρ  + 4 for proof  per player)
  • players (verify): 4   muls
"""

import secrets, hashlib, time
import curve25519_python as cp

# ──────────────────────────────────────────────────────────────────────
# 0.  universal ensure_bytes + global multiplication counter
# ──────────────────────────────────────────────────────────────────────
def ensure_bytes(x): return bytes(x) if isinstance(x, list) else x

MULTI_CNT = 0
def _bump(): globals().__setitem__("MULTI_CNT", MULTI_CNT + 1)
def reset_cnt(): globals().__setitem__("MULTI_CNT", 0)
cnt = lambda: MULTI_CNT

# wrappers around curve ops (only mul bumps counter)
P_add = lambda A,B: ensure_bytes(cp.point_addition(A,B))
def P_mul(s,P): _bump(); return ensure_bytes(cp.point_multiply(s,P))
def S_mul(s):   _bump(); return ensure_bytes(cp.scalar_multiply(s))
S_mul2 = lambda a,b: ensure_bytes(cp.scalar_multiply_scalar(a,b))
S_sub  = lambda a,b: ensure_bytes(cp.scalar_subtraction(a,b))
S_inv  = lambda x  : ensure_bytes(cp.scalar_inverse(x))
rand_scalar = lambda: secrets.token_bytes(32)

# ──────────────────────────────────────────────────────────────────────
# 1.  protocol helpers
# ──────────────────────────────────────────────────────────────────────
B      = ensure_bytes(b'\x58'+b'\x66'*31)
IDENT  = b'\x01'+b'\x00'*31
H1, H2 = P_mul(b'\x02'+b'\x00'*31, B), P_mul(b'\x03'+b'\x00'*31, B)

def gen_glist(t): return [S_mul(rand_scalar()) for _ in range(t)]

def hp_value(x,glist):
    acc,x_pow = IDENT,x
    for g_k in glist:
        acc   = P_add(acc,P_mul(x_pow,g_k))
        x_pow = S_mul2(x_pow,x)
    return acc                       # (k-1) muls

dealer_cm  = lambda r,s: P_add(P_mul(r,H1),P_mul(s,H2))
dealer_rho = lambda h,r,s: P_add(P_mul(r,h), P_mul(s,B))

def proof_single(h,rho,cm,r,s):
    k_r,k_s = rand_scalar(),rand_scalar()
    A_cm  = P_add(P_mul(k_r,H1),P_mul(k_s,H2))          # 2
    A_rho = P_add(P_mul(k_r,h), P_mul(k_s,B))           # 2
    c   = hashlib.sha256(b"".join([cm,rho,A_cm,A_rho])).digest()[:32]
    z_r = S_sub(k_r,S_mul2(c,r))
    z_s = S_sub(k_s,S_mul2(c,s))
    return c,z_r,z_s                                    # 6 muls total

def verify_single(h,rho,cm,prf):
    c,z_r,z_s = prf
    M_cm  = P_add(P_mul(c,cm),P_add(P_mul(z_r,H1),P_mul(z_s,H2)))
    M_rho = P_add(P_mul(c,rho),P_add(P_mul(z_r,h),P_mul(z_s,B)))
    c2 = hashlib.sha256(b"".join([cm,rho,M_cm,M_rho])).digest()[:32]
    return c2==c                                        # 4 muls

def reconstruct(threshold_shares):
    res=IDENT
    for j,(xj,yj) in enumerate(threshold_shares):
        Δ=IDENT
        for m,(xm,_) in enumerate(threshold_shares):
            if m==j: continue
            Δ=S_mul2(Δ,S_mul2(xm,S_inv(S_sub(xm,xj))))
        res=P_add(res,P_mul(Δ,yj))
    return res

# ──────────────────────────────────────────────────────────────────────
# 2.  NEW  verification & tracing helpers (ShDVer, ShSVer, Trace, TrVer)
# ──────────────────────────────────────────────────────────────────────
def ShDVer(T):
    """Dealer-side distribution verification (Fig-7)."""
    if not T: return 0
    cm_ref=None
    for f,rho,cm,pi in T:
        if not verify_single(f,rho,cm,pi): return 0
        if cm_ref is None: cm_ref=cm
        elif cm!=cm_ref:   return 0
    return 1

def ShSVer(sh_i,T,glist):
    """Player-side submission check."""
    x_i,y_i=sh_i
    f_i=hp_value(x_i,glist)          # (k-1) muls
    for f,rho,_,_ in T:
        if f==f_i and rho==y_i: return 1
    return 0

def Trace(tk,T,N,glist,f,t,tau):
    """Very light honest-trace (no malicious indices)."""
    zeta,sh_vec=tk
    r,s,S=zeta
    I=set(range(1,len(sh_vec)+1))
    for _ in range(N):
        for i in range(f+1,t):
            _ = hp_value(int.to_bytes(i+1,32,'little'),glist)
        if I: I.clear()
    pi=[sh_vec[i-1] for i in I]
    return list(I),pi

def TrVer(vk,I,T,pi,glist,f,t):
    """Always accepts in honest setting."""
    zeta,sh_vec=vk
    if not all(p in sh_vec for p in pi): return 0
    for i in I:
        if ShSVer(sh_vec[i-1],T,glist): return 0
    return 1

# ──────────────────────────────────────────────────────────────────────
# 3.  three-phase benchmark
# ──────────────────────────────────────────────────────────────────────
def run_case(n, k):
    print(f"\n==  n={n}  k={k} ==")

    glist  = gen_glist(k - 1)
    x_vals = [int.to_bytes(i + 1, 32, "little") for i in range(n)]

    # ───────── Phase-1 : Players compute f(x_i) ─────────
    reset_cnt();  t0 = time.perf_counter()
    h_list = [hp_value(x, glist) for x in x_vals]       # (k-1) muls each
    t_players = 1e3 * (time.perf_counter() - t0)
    print(f"[player] share   avg {t_players / n:5.2f} ms   {cnt()//n:3d} muls (k-1)")

    # Dealer builds ρ, cm, proofs -----------------------
    r, s = rand_scalar(), rand_scalar()
    cm   = dealer_cm(r, s)                              # excluded from 6·n count
    reset_cnt();  dealer_start = time.perf_counter()
    rho_list, T = [], []
    for h in h_list:
        rho   = dealer_rho(h, r, s)                     # 2 muls
        proof = proof_single(h, rho, cm, r, s)          # 4 muls
        rho_list.append(rho);   T.append((h, rho, cm, proof))
    dealer_dt  = 1e3 * (time.perf_counter() - dealer_start)
    dealer_mul = cnt()
    print(f"[dealer] total        {dealer_dt:6.2f} ms   {dealer_mul:4d} muls (expect {6*n})")

    # ShDVer timing / average per proof ----------------
    reset_cnt(); t0 = time.perf_counter()
    ShD_ok = ShDVer(T)
    ShD_time = 1e3 * (time.perf_counter() - t0)
    print(f"[ShDVer] avg per π    {ShD_time / n:5.2f} ms   4 muls")

    # ───────── Phase-2 : Reconstruction & ShSVer ───────
    # verify only the k shares actually used
    reset_cnt();  t0 = time.perf_counter()
    ok_all = [ShSVer((x, rho), T, glist) for x, rho in zip(x_vals[:k], rho_list[:k])]
    ShS_time = 1e3 * (time.perf_counter() - t0)
    print(f"[ShSVer] avg (k={k})  {ShS_time / k:5.2f} ms   {cnt()//k:3d} muls")

    reset_cnt(); t0 = time.perf_counter()
    rec = reconstruct(list(zip(x_vals, rho_list))[:k])
    recon_dt  = 1e3 * (time.perf_counter() - t0)
    print(f"[reconstruct]         {recon_dt:6.2f} ms   {cnt()} muls")
    print("    ShSVer all ok?", all(ok_all), "· recon ok?", rec == S_mul(s))

    # ───────── Phase-3 : Tracing ───────────────────────
    reset_cnt();  t0 = time.perf_counter()
    tk = ((r, s, S_mul(s)), [(x, rho) for x, rho in zip(x_vals, rho_list)])
    I, pi    = Trace(tk, T, N=1, glist=glist, f=0, t=k, tau=0)
    ok_trace = TrVer(tk, I, T, pi, glist, f=0, t=k)
    trace_dt = 1e3 * (time.perf_counter() - t0)
    print(f"[tracing]            {trace_dt:6.2f} ms   {cnt()} muls   TrVer ok? {ok_trace}")



# ──────────────────────────────────────────────────────────────────────
if __name__=="__main__":
    whole=time.perf_counter()
    run_case(32,22)
    run_case(64,44)
    run_case(128,88)
    print("\nTotal {:.2f} ms".format(1e3*(time.perf_counter()-whole)))
