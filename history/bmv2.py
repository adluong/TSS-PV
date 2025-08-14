#!/usr/bin/env python3
"""
Tiny-TSS-PV v13 ▸ benchmark with time / muls / heap (tracemalloc) usage
────────────────────────────────────────────────────────────────────
Implements TSS-PV Ver 1.0 with ElGamal encryption for dummy shares.
Fixes TrVer and ShD to correctly check cm consistency using T[i][7].
"""

import secrets, hashlib, time, random, gc, tracemalloc, curve25519_python as cp

# ────────── start heap tracing (very cheap: 1-2 µs / snapshot) ──────────
tracemalloc.start()
def heap_kib():
    cur, peak = tracemalloc.get_traced_memory()
    return cur // 1024, peak // 1024
def banner():
    gc.collect()
    cur, peak = heap_kib()
    print(f"   [heap] {cur:6d} KiB  (peak {peak} KiB)")

# ────────── normalise to bytes ──────────
E = lambda v: v if isinstance(v, (bytes, bytearray)) else bytes(v)

# ────────── global counters ──────────
CNT = RB_CNT = RBOX_CNT = 0  # Added RBOX_CNT for Rbox-specific multiplications
bump = lambda: globals().__setitem__("CNT", CNT + 1)
rb_bump = lambda d: globals().__setitem__("RB_CNT", RB_CNT + d)
rbox_bump = lambda d: globals().__setitem__("RBOX_CNT", RBOX_CNT + d)
reset = lambda: globals().__setitem__("CNT", 0)
rb_reset = lambda: globals().__setitem__("RB_CNT", 0)
rbox_reset = lambda: globals().__setitem__("RBOX_CNT", 0)

# ────────── curve wrappers ──────────
B = b'\x58' + b'\x66'*31
P = lambda A, B: E(cp.point_addition(A, B))
Pm = lambda s, P_: bump() or E(cp.point_multiply(s, P_))
Sm = lambda s: bump() or E(cp.scalar_multiply(s))
S2 = lambda a, b: E(cp.scalar_multiply_scalar(a, b))
Inv = lambda z: E(cp.scalar_inverse(z))
Sub = lambda a, b: E(cp.scalar_subtraction(a, b))
rand = lambda: secrets.token_bytes(32)
cm = lambda r, s: P(Pm(r, H1), Pm(s, H2))
rho = lambda h, r, s: P(Pm(r, h), Pm(s, B))
_sample = lambda x, g, r, s: (x, rho(hp(x, g), r, s))

ONE = b'\x01' + b'\0'*31
H1, H2 = Pm(b'\x02'+b'\0'*31, B), Pm(b'\x03'+b'\0'*31, B)

# ────────── ElGamal encryption ──────────
def elgamal_encrypt(m, pk, sk=None):
    k = rand()
    c1 = Pm(k, B)
    c2 = P(m, Pm(k, pk))
    if sk:
        return (c1, c2, sk)
    return (c1, c2)

def elgamal_decrypt(c1, c2, sk):
    temp = Pm(sk, c1)
    return P(c2, Pm(Inv(sk), temp))

# ────────── polynomial & sharing helpers ──────────
gpoly = lambda t: [Sm(rand()) for _ in range(t)]
def hp(x, g):
    acc, xp = ONE, x
    for gk in g:
        acc = P(acc, Pm(xp, gk)); xp = S2(xp, x)
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
    acc = E(cp.scalar_multiply(b"\0" * 32))
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

# ────────── verification helpers ──────────
def ShD(T, cm):
    if any(t[7] != cm for t in T):
        return 0
    for fx_i, fz_i, px_i, pz_i, π_i, πp_i, ct_i, _ in T:
        if not (check(fx_i, px_i, cm, π_i) and check(fz_i, pz_i, cm, πp_i)):
            return 0
    return 1

def ShS(sh, T, g):
    x, y = sh
    fx = hp(x, g)
    return int(any((fx == fx_j and y == px_j) or (fx == fz_j and y == pz_j) for fx_j, fz_j, px_j, pz_j, _, _, _, _ in T))

# ────────── reconstruction oracle with ShS and ShD checks ──────────
def Rbox(k, embeds, T, g, cm):
    def R(shares):
        shares = [shares] if isinstance(shares, tuple) else list(shares)
        all_shares = shares
        valid_shares = [] + embeds
        for x, y in all_shares:
            if not ShS((x, y), T, g):
                return None
            valid_shares.append((x, y))
        uniq = {}; [uniq.setdefault(x, y) for x, y in valid_shares]
        if len(uniq) < k:
            return None
        before, t0 = CNT, time.perf_counter()
        res = recon(list(uniq.items())[:k])
        rbox_bump(CNT - before)  # Track Rbox-specific multiplications
        globals().__setitem__("RBOX_TIME", globals().get("RBOX_TIME", 0.0) + (time.perf_counter() - t0))  # Track Rbox-specific time
        return res
    return R

# ────────── Trace / TrVer ──────────
def Trace(tk, T, g, f, t, R, cm):
    ζ, dsh, shv = tk; r, s, S = ζ; n = len(shv); I = set(range(n))
    DSH = []
    banned = {x for x, _ in shv}
    for z, pz in random.sample(dsh, t - f - 1):
        if z not in banned:
            DSH.append((z, pz))
            banned.add(z)
    for idx, sh in enumerate(shv):
        if R(DSH + [sh]) == S:
            I.discard(idx)
    return sorted(I), [shv[i] for i in I]

def TrVer(vk, I, T, π, g, f, t, R, cm):
    ζ, dsh, shv = vk; r, s, S = ζ
    if any(t[7] != cm for t in T):
        return 0
    if not all(sh in shv for sh in π):
        return 0
    for x, px in shv:
        fx = hp(x, g)
        if px != rho(fx, r, s):
            return 0
    for z, pz in dsh:
        fz = hp(z, g)
        if pz != rho(fz, r, s):
            return 0
    DSH = []
    banned = {x for x, _ in shv}
    for z, pz in random.sample(dsh, t - f - 1):
        if z not in banned:
            DSH.append((z, pz))
            banned.add(z)
    for idx in I:
        x, y = shv[idx]
        if not any((y == px_j and hp(x, g) == fx_j) or (y == pz_j and hp(x, g) == fz_j) for fx_j, fz_j, px_j, pz_j, _, _, _, _ in T):
            return 0
        rbox_reset()  # Reset Rbox counter for this call
        rbox_t0 = time.perf_counter()  # Start Rbox timing
        if R(DSH + [shv[idx]]) == S:
            return 0
        globals().__setitem__("RBOX_TIME_TRVER", globals().get("RBOX_TIME_TRVER", 0.0) + (time.perf_counter() - rbox_t0))  # Accumulate Rbox time
        rbox_bump(CNT)  # Track Rbox multiplications for this call
    return 1

# ────────── benchmark harness ──────────
def run(n, k, f):
    assert 0 < f < k - 1
    print(f"\n== n={n} k={k} f={f} ==")
    g = gpoly(k - 1)
    xs = [(i + 1).to_bytes(32, 'little') for i in range(n)]
    zs = [rand() for _ in range(n)]

    reset(); t0 = time.perf_counter()
    fx = [hp(x, g) for x in xs]
    fz = [hp(z, g) for z in zs]
    print(f"[share ] {1e3*(time.perf_counter()-t0)/n:7.1f} ms {CNT//n:3d} mul/ply"); banner()

    r, s = rand(), rand(); S = Sm(s); cm_ = cm(r, s)

    sk = rand(); pk = Pm(sk, B)
    reset(); t0 = time.perf_counter()
    px, pz, T = [], [], []
    dsh = []
    for i in range(n):
        px_i = rho(fx[i], r, s)
        π_i = proof(fx[i], px_i, cm_, r, s)
        pz_i = rho(fz[i], r, s)
        ct_i = elgamal_encrypt(fz[i], pk, sk)
        πp_i = proof(fz[i], pz_i, cm_, r, s)
        T.append((fx[i], fz[i], px_i, pz_i, π_i, πp_i, ct_i[:2], cm_))
        px.append(px_i)
        pz.append(pz_i)
        dsh.append((zs[i], pz_i))
    shv = [(xs[i], T[i][2]) for i in range(n)]
    print(f"[dealer] {1e3*(time.perf_counter()-t0):7.1f} ms {CNT:4d} mul"); banner()

    reset(); t0 = time.perf_counter()
    ok = ShD(T, cm_)
    print(f"[ShD ] {1e3*(time.perf_counter()-t0):7.1f} ms {CNT//n:3d} mul ok={ok}"); banner()

    reset(); t0 = time.perf_counter()
    allok = all(ShS(sh, T, g) for sh in shv[:k])
    print(f"[ShS ] {1e3*(time.perf_counter()-t0):7.1f} ms {CNT//k:3d} mul all={allok}"); banner()

    rb_reset(); globals().__setitem__("RB_TIME", 0.0)
    R = Rbox(k, random.sample(shv, f), T, g, cm_)
    reset(); t0 = time.perf_counter()
    secret = S
    reconstructed = recon(shv[:k])
    ok = reconstructed == secret
    print(f"[Recon ] {1e3*(time.perf_counter()-t0):7.1f} ms {CNT:4d} mul ok={ok}"); banner()

    tk = ((r, s, S), dsh, shv)
    rb_before, c_before, t0 = RB_CNT, CNT, time.perf_counter()
    I, π = Trace(tk, T, g, f, k, R, cm_)
    dt = time.perf_counter() - t0
    trace_mul = CNT - c_before
    trace_ms = 1e3 * dt
    rbox_trace_mul = RBOX_CNT  # Rbox-specific multiplications for Trace
    rbox_trace_time = globals().get("RBOX_TIME", 0.0)  # Rbox-specific time for Trace
    print(f"[Trace ] {trace_ms:7.1f} ms {trace_mul:4d} mul |I|={len(I)}"); banner()
    print(f"[Rbox_in_Trace ] {1e3*rbox_trace_time:7.1f} ms {rbox_trace_mul:4d} mul"); banner()

    rb_before, c_before, t0 = RB_CNT, CNT, time.perf_counter()
    ok = TrVer(tk, I, T, π, g, f, k, R, cm_)
    dt = time.perf_counter() - t0
    trv_mul = CNT - c_before
    trv_ms = 1e3 * dt
    rbox_trver_mul = RBOX_CNT  # Rbox-specific multiplications for TrVer
    rbox_trver_time = globals().get("RBOX_TIME_TRVER", 0.0)  # Rbox-specific time for TrVer
    print(f"[TrVer ] {trv_ms:7.1f} ms {trv_mul:4d} mul ok={ok}"); banner()
    print(f"[Rbox_in_TrVer ] {1e3*rbox_trver_time:7.1f} ms {rbox_trver_mul:4d} mul"); banner()

# ────────── entry ──────────
if __name__ == "__main__":
    for n, k, f in [
        (32, 17, 11),
        (64, 33, 22), (128, 65, 43),
        (256, 129, 86)
        ]:
        run(n, k, f)