#!/usr/bin/env python3
"""
Tiny-TSS-PV v13 ▸ benchmark with time / muls / heap (tracemalloc) usage
────────────────────────────────────────────────────────────────────
Implements TSS-PV Ver 1.0 with ElGamal encryption for dummy shares.
Fixes TrVer and ShD to correctly check cm consistency using T[i][7].
"""

import secrets, hashlib, time, random, gc, tracemalloc, curve25519_python as cp

# ────────── Heap tracing ──────────
tracemalloc.start()
def heap_kib():
    cur, peak = tracemalloc.get_traced_memory()
    return cur // 1024, peak // 1024
def banner():
    gc.collect()
    cur, peak = heap_kib()
    print(f"   [heap] {cur:6d} KiB  (peak {peak} KiB)")

# ────────── Normalise to bytes ──────────
def E(v):
    return v if isinstance(v, (bytes, bytearray)) else bytes(v)

# ────────── Global counters ──────────
CNT = RB_CNT = 0
def bump():
    global CNT
    CNT += 1
def rb_bump(d):
    global RB_CNT
    RB_CNT += d
def reset():
    global CNT
    CNT = 0
def rb_reset():
    global RB_CNT
    RB_CNT = 0

# ────────── Curve wrappers ──────────
B = b'\x58' + b'\x66'*31
def P(A, B):
    return E(cp.point_addition(A, B))
def Pm(s, P_):
    bump()
    return E(cp.point_multiply(s, P_))
def Sm(s):
    bump()
    return E(cp.scalar_multiply(s))
def S2(a, b):
    return E(cp.scalar_multiply_scalar(a, b))
def Inv(z):
    return E(cp.scalar_inverse(z))
def Sub(a, b):
    return E(cp.scalar_subtraction(a, b))
def rand():
    return secrets.token_bytes(32)
def cm(r, s):
    return P(Pm(r, H1), Pm(s, H2))
def rho(h, r, S):
    return P(Pm(r, h), S)
def _sample(x, g, r, S):
    return (x, rho(hp(x, g), r, S))

ONE = b'\x01' + b'\0'*31
H1, H2 = Pm(b'\x02'+b'\0'*31, B), Pm(b'\x03'+b'\0'*31, B)
ID = Sm(b"\0" * 32)

# ────────── ElGamal encryption ──────────
def elgamal_encrypt(m, pk, sk=None):
    """Encrypt m (point) under pk (public key point) on Curve25519."""
    k = rand()
    c1 = Pm(k, B)
    c2 = P(m, Pm(k, pk))
    return (c1, c2, sk) if sk else (c1, c2)

# ────────── Polynomial & sharing helpers ──────────
def gpoly(t):
    return [Sm(rand()) for _ in range(t)]
def hp(x, g):
    acc, xp = ONE, x
    for gk in g:
        acc = P(acc, Pm(xp, gk))
        xp = S2(xp, x)
    return acc

def proof(h, p, cm_, r, s):
    kr, ks = rand(), rand()
    Acm, Ar = P(Pm(kr, H1), Pm(ks, H2)), P(Pm(kr, h), Pm(ks, B))
    c = hashlib.sha256(b''.join([cm_, p, Acm, Ar])).digest()[:32]
    zr, zs = Sub(kr, S2(c, r)), Sub(ks, S2(c, s))
    return c, zr, zs

def check(h, p, cm_, π):
    c, zr, zs = π
    Mc = P(Pm(c, cm_), P(Pm(zr, H1), Pm(zs, H2)))
    Mr = P(Pm(c, p), P(Pm(zr, h), Pm(zs, B)))
    return hashlib.sha256(b''.join([cm_, p, Mc, Mr])).digest()[:32] == c

def montgomery_batch_invert(values):
    k = len(values)
    partials = [ONE]
    for v in values:
        partials.append(S2(partials[-1], v))
    inv_total = Inv(partials[-1])
    invs = [None] * k
    for i in range(k-1, -1, -1):
        invs[i] = S2(partials[i], inv_total)
        inv_total = S2(inv_total, values[i])
    return invs

def recon(shs):
    acc = ID
    k = len(shs)
    denoms = [ONE] * k
    for j in range(k):
        for m in range(k):
            if m != j:
                denoms[j] = S2(denoms[j], Sub(shs[m][0], shs[j][0]))
    inv_denoms = montgomery_batch_invert(denoms)
    for j, (xj, yj) in enumerate(shs):
        l = inv_denoms[j]
        for m in range(k):
            if m != j:
                l = S2(l, shs[m][0])
        acc = P(acc, Pm(l, yj))
    return acc

# ────────── Verification helpers ──────────
def ShD(T, cm):
    if any(t[7] != cm for t in T):
        return False
    for fx_i, fz_i, px_i, pz_i, π_i, πp_i, _, _ in T:
        if not (check(fx_i, px_i, cm, π_i) and check(fz_i, pz_i, cm, πp_i)):
            return False
    return True

def ShS(sh, T, g):
    x, y = sh
    fx = hp(x, g)
    for fx_j, fz_j, px_j, pz_j, π_j, πp_j, _, cm_j in T:
        if (fx == fx_j and y == px_j):
            if check(fx, px_j, cm_j, π_j):
                return True
        elif (fx == fz_j and y == pz_j):
            if check(fx, pz_j, cm_j, πp_j):
                return True
    return False  # No match found

# ────────── Reconstruction oracle ──────────
def Rbox(k, embeds, T, g, cm):
    def R(shares):
        shares = [shares] if isinstance(shares, tuple) else list(shares)
        valid_shares = embeds[:]
        for x, y in shares: #t-f shares * t+5 muls
            if not ShS((x, y), T, g):
                return None
            valid_shares.append((x, y))
        uniq = dict(valid_shares)
        if len(uniq) < k:
            return None
        before, t0 = CNT, time.perf_counter()
        res = recon(list(uniq.items())[:k]) #t shares
        # rb_bump(CNT - before)
        globals()['RB_TIME'] = globals().get('RB_TIME', 0.0) + (time.perf_counter() - t0)
        return res
    return R

# ────────── Helper for dummy share selection ──────────
def select_dummy_shares(dsh, shv, t, f):
    banned = {x for x, _ in shv}
    return [(z, pz) for z, pz in random.sample(dsh, t - f - 1) if z not in banned]

# ────────── Trace / TrVer ──────────
def Trace(tk, T, g, f, t, R, cm):
    ζ, dsh, shv = tk
    r, s, S = ζ
    I = set(range(len(shv)))
    DSH = select_dummy_shares(dsh, shv, t, f)
    for idx, sh in enumerate(shv):
        if R(DSH + [sh]) == S:
            I.discard(idx)
    return sorted(I), [shv[i] for i in I]

def TrVer(vk, I, T, π, g, f, t, R, cm_):
    ζ, dsh, shv = vk
    r, s, S = ζ
    if any(t[7] != cm_ for t in T):
        return False
    if cm_ != cm(r, s) or S != Sm(s): #3 muls
        return False
    if not all(sh in shv for sh in π):
        return False
    for x, px in shv:
        if px != rho(hp(x, g), r, S): #nt muls
            return False
    for z, pz in dsh:
        if pz != rho(hp(z, g), r, S): #nt muls
            return False
    DSH = select_dummy_shares(dsh, shv, t, f)
    for idx in I:
        x, y = shv[idx]
        if not any((y == px_j) or (y == pz_j) for _, _, px_j, pz_j, _, _, _, _ in T):
            return False
        if R(DSH + [shv[idx]]) == S: #(t+5)muls per ShS on 1 share. There are t-f shares fed to R
            return False
    return True

# ────────── Benchmark harness ──────────
def run(n, k, f):
    assert 0 < f < k - 1
    print(f"\n== n={n} k={k} f={f} ==")
    g = gpoly(k - 1)
    xs = [(i + 1).to_bytes(32, 'little') for i in range(n)]

    # Player share generation
    reset()
    t0 = time.perf_counter()
    fx = [hp(x, g) for x in xs]
    print(f"[share ] {1e3*(time.perf_counter()-t0)/n:7.1f} ms {CNT//n:3d} mul/ply")
    banner()

    # Dealer operations
    r, s = rand(), rand()
    S = Sm(s)
    sk = rand()
    pk = Pm(sk, B)
    reset()
    cm_ = cm(r, s)
    t0 = time.perf_counter()
    zs = [rand() for _ in range(n)]
    fz = [hp(z, g) for z in zs]
    T = []
    dsh = []
    for i in range(n):
        px_i = rho(fx[i], r, S)
        π_i = proof(fx[i], px_i, cm_, r, s)
        pz_i = rho(fz[i], r, S)
        πp_i = proof(fz[i], pz_i, cm_, r, s)
        T.append((fx[i], fz[i], px_i, pz_i, π_i, πp_i, (), cm_))
        dsh.append((zs[i], pz_i))
    shv = [(x, t[2]) for x, t in zip(xs, T)]
    print(f"[dealer] {1e3*(time.perf_counter()-t0):7.1f} ms {CNT:4d} mul")
    banner()

    # VerifySD
    reset()
    t0 = time.perf_counter()
    ok = ShD(T, cm_)
    print(f"[ShD ] {1e3*(time.perf_counter()-t0):7.1f} ms {CNT:4d} mul ok={ok}")
    banner()

    # VerifySS
    reset()
    t0 = time.perf_counter()
    allok = all(ShS(sh, T, g) for sh in shv[:k])
    print(f"[ShS ] {1e3*(time.perf_counter()-t0):7.1f} ms {CNT:4d} mul all={allok}")
    banner()

    # Reconstruction
    rb_reset()
    globals()['RB_TIME'] = 0.0
    R = Rbox(k, random.sample(shv, f), T, g, cm_)
    reset()
    t0 = time.perf_counter()
    reconstructed = recon(shv[:k])
    ok = reconstructed == S
    print(f"[Recon ] {1e3*(time.perf_counter()-t0):7.1f} ms {CNT:4d} mul ok={ok}")
    banner()

    # Trace
    reset()
    tk = ((r, s, S), dsh, shv)
    t0 = time.perf_counter()
    I, π = Trace(tk, T, g, f, k, R, cm_)
    trace_mul = CNT
    trace_ms = 1e3 * (time.perf_counter() - t0)
    print(f"[Trace ] {trace_ms:7.1f} ms {trace_mul:4d} mul |I|={len(I)}")
    banner()

    # TrVer
    reset()
    t0 = time.perf_counter()
    ok = TrVer(tk, I, T, π, g, f, k, R, cm_)
    trv_mul = CNT
    trv_ms = 1e3 * (time.perf_counter() - t0)
    print(f"[TrVer ] {trv_ms:7.1f} ms {trv_mul:4d} mul ok={ok}")
    banner()

# ────────── Entry ──────────
if __name__ == "__main__":
    for n, k, f in [(32, 17, 11), (64, 33, 22), (128, 65, 43), (256, 129, 85)]:
        run(n, k, f)