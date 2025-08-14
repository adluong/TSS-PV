#tiny_tss_pv_c25519/1.3/fix missing generator+identity; finish bench(); minor tidy
#!/usr/bin/env python3
"""TSS-PV micro-benchmark – Curve25519 backend (rev 1.3)"""

import secrets, hashlib, time, random, gc, tracemalloc, curve25519_python as cp

# ───── helpers ─────
tracemalloc.start()
heap = lambda: print(f"   [heap] {tracemalloc.get_traced_memory()[0]//1024:6d} KiB")
def E(v, l=32):
    if isinstance(v, (bytes, bytearray)) and len(v) == l:
        return v
    if isinstance(v, list) and len(v) == l and all(0 <= b < 256 for b in v):
        return bytes(v)
    raise TypeError("expected 32-byte bytes or list")

# ───── global counters ─────
CNT = RB_CNT = 0
bump     = lambda: globals().__setitem__("CNT", CNT + 1)
rb_bump  = lambda d: globals().__setitem__("RB_CNT", RB_CNT + d)
reset    = lambda: globals().__setitem__("CNT", 0)
rb_reset = lambda: globals().__setitem__("RB_CNT", 0)

# ───── constants (32 B) ─────
ZERO = b"\0" * 32
ONE  = b"\x01" + b"\0" * 31
B    = E(cp.scalar_multiply(ONE))   # basepoint
ID   = E(cp.scalar_multiply(ZERO))  # point at infinity

# ───── curve wrappers ─────
P   = lambda A, B_: E(cp.point_addition(E(A), E(B_)))
Pm  = lambda s, P_: bump() or E(cp.point_multiply(E(s), E(P_)))
Sm  = lambda s: bump() or E(cp.scalar_multiply(E(s)))
S2  = lambda a, b: E(cp.scalar_multiply_scalar(E(a), E(b)))
Inv = lambda z: E(cp.scalar_inverse(E(z)))
Sub = lambda a, b: E(cp.scalar_subtraction(E(a), E(b)))
rand = lambda: secrets.token_bytes(32)

H1, H2 = Pm(b"\x02" + b"\0" * 31, B), Pm(b"\x03" + b"\0" * 31, B)

# ───── protocol helpers ─────
gpoly = lambda t: [ID] + [Sm(rand()) for _ in range(t - 1)]
def hp(x, g):
    acc, xp = ID, ONE
    for gk in g:
        acc = P(acc, Pm(xp, gk))
        xp  = S2(xp, x)
    return acc

cm  = lambda r, s: P(Pm(r, H1), Pm(s, H2))
rho = lambda h, r, s: P(Pm(r, h), Pm(s, B))

def proof(h, ρ, cm_, r, s):
    kr, ks = rand(), rand()
    πc     = hashlib.sha256(
               b"".join([cm_, ρ,
                         P(Pm(kr, H1), Pm(ks, H2)),
                         P(Pm(kr, h),  Pm(ks, B))])
             ).digest()[:32]
    return πc, Sub(kr, S2(πc, r)), Sub(ks, S2(πc, s))

def check(h, ρ, cm_, π):
    c, zr, zs = π
    Mc = P(Pm(c, cm_), P(Pm(zr, H1), Pm(zs, H2)))
    Mr = P(Pm(c, ρ ), P(Pm(zr, h ), Pm(zs, B )))
    return hashlib.sha256(b"".join([cm_, ρ, Mc, Mr])).digest()[:32] == c

def recon(shs):                                  # naive O(k²)
    acc = ID
    for j, (xj, yj) in enumerate(shs):
        l = ONE
        for m, (xm, _) in enumerate(shs):
            if m == j: continue
            l = S2(l, S2(xm, Inv(Sub(xm, xj))))
        acc = P(acc, Pm(l, yj))
    return acc

ShD = lambda T: int(all(check(*row) for row in T))

def ShS(sh, T, g):
    x, y = sh; fx = hp(x, g)
    return int(any(fx == f and y == ρ for f, ρ, *_ in T))

_sample = lambda x, g, r, s: (x, rho(hp(x, g), r, s))

# ───── reconstruction oracle ─────
def Rbox(k, emb):
    def R(extra):
        extra = [extra] if isinstance(extra, tuple) else list(extra)
        uniq  = dict(emb + extra)
        if len(uniq) < k: return None
        before, t0 = CNT, time.perf_counter()
        res = recon(list(uniq.items())[:k])
        rb_bump(CNT - before)
        globals().__setitem__(
            "RB_TIME",
            globals().get("RB_TIME", 0.0) + (time.perf_counter() - t0),
        )
        return res
    return R

def Trace(tk, T, g, f, t, R):
    ζ, shv = tk
    r, s, S = ζ
    I       = set(range(len(shv)))
    banned  = {x for x, _ in (R([]) or [])} | {x for x, _ in shv}
    F = []
    while len(F) < t - f - 1:
        x = rand()
        if x not in banned:
            banned.add(x); F.append(_sample(x, g, r, s))
    for idx, sh in enumerate(shv):
        if R(F + [sh]) == S:
            I.discard(idx)
    return sorted(I), [shv[i] for i in I]

def TrVer(vk, I, T, π, g, f, t, R):
    ζ, shv = vk
    r, s, S = ζ
    if any(_sample(x, g, r, s) != (x, ρ) for x, ρ in shv): return 0
    banned = {x for x, _ in (R([]) or [])} | {x for x, _ in shv}
    F = []
    while len(F) < t - f - 1:
        x = rand()
        if x not in banned:
            banned.add(x); F.append(_sample(x, g, r, s))
    for idx in I:
        if R(F + [shv[idx]]) == S and ShS(shv[idx], T, g):
            return 0
    return 1

# ───── benchmark driver ─────
def bench(n, k, f):
    assert 0 < f < k - 1
    print(f"\n== n={n} k={k} f={f} ==")
    g  = gpoly(k)
    xs = [(37821381 + i).to_bytes(32, "big") for i in range(n)]

    reset(); t0 = time.perf_counter()
    hs = [hp(x, g) for x in xs]
    print(f"[share ] {1e3*(time.perf_counter()-t0)/n:7.1f} ms  {CNT//n:3d} mul/ply"); heap()

    r, s = rand(), rand()
    cm_  = cm(r, s)
    reset(); t0 = time.perf_counter()
    ρs, T = [], []
    for h in hs:
        ρ = rho(h, r, s); π = proof(h, ρ, cm_, r, s)
        ρs.append(ρ); T.append((h, ρ, cm_, π))
    print(f"[dealer] {1e3*(time.perf_counter()-t0):7.1f} ms  {CNT:4d} mul"); heap()

    reset(); t0 = time.perf_counter(); ok = ShD(T)
    print(f"[ShD   ] {1e3*(time.perf_counter()-t0):7.1f} ms  {CNT//n:3d} mul  ok={ok}"); heap()

    reset(); t0 = time.perf_counter()
    allok = all(ShS((x, ρ), T, g) for x, ρ in zip(xs[:k], ρs[:k]))
    print(f"[ShS   ] {1e3*(time.perf_counter()-t0):7.1f} ms  {CNT//k:3d} mul-per-share  all={allok}"); heap()

    rb_reset(); globals().__setitem__("RB_TIME", 0.0)
    R = Rbox(k, random.sample(list(zip(xs, ρs)), f))

    reset(); t0 = time.perf_counter(); ok = recon(list(zip(xs, ρs))[:k]) == Sm(s)
    print(f"[Recon ] {1e3*(time.perf_counter()-t0):7.1f} ms  {CNT:4d} mul  ok={ok}"); heap()

    tk = ((r, s, Sm(s)), list(zip(xs, ρs)))

    rb0, c0, t0 = RB_CNT, CNT, time.perf_counter()
    I, π = Trace(tk, T, g, f, k, R)
    dt = time.perf_counter() - t0
    print(f"[Trace ] {1e3*(dt-globals()['RB_TIME']):7.1f} ms  {CNT-c0-(RB_CNT-rb0):4d} mul  |I|={len(I)}")
    heap(); tmp = globals()["RB_TIME"]

    rb0, c0, t0 = RB_CNT, CNT, time.perf_counter()
    ok = TrVer(tk, I, T, π, g, f, k, R)
    dt = time.perf_counter() - t0
    print(f"[TrVer ] {1e3*(dt-(globals()['RB_TIME']-tmp)):7.1f} ms  {CNT-c0-(RB_CNT-rb0):4d} mul  ok={ok}"); heap()

    print(f"[R-box ] {1e3*globals()['RB_TIME']:7.1f} ms  {RB_CNT:4d} mul"); heap()

# ───── main ─────
if __name__ == "__main__":
    for n, k, f in [(32, 17, 4), (64, 33, 10), (128, 65, 20), (256, 129, 40), (512, 257, 80)]:
        bench(n, k, f)
