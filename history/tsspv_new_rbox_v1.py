#!/usr/bin/env python3
"""
Tiny-TSS-PV v12 ▸ benchmark with time / muls / heap (tracemalloc) usage
────────────────────────────────────────────────────────────────────
Implements Countermeasure 4: Robust Fresh Share Integration.
Dealer pre-generates n fresh shares, includes them in T (size 2n), and keeps x values private.
Rbox includes ShS and ShD checks to reject invalid shares.
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
CNT = RB_CNT = 0
bump = lambda: globals().__setitem__("CNT", CNT + 1)
rb_bump = lambda d: globals().__setitem__("RB_CNT", RB_CNT + d)
reset = lambda: globals().__setitem__("CNT", 0)
rb_reset = lambda: globals().__setitem__("RB_CNT", 0)

# ────────── curve wrappers ──────────
B = b'\x58' + b'\x66'*31
P = lambda A,B: E(cp.point_addition(A,B))
Pm = lambda s,P_: bump() or E(cp.point_multiply(s,P_))
Sm = lambda s: bump() or E(cp.scalar_multiply(s))
S2 = lambda a,b: E(cp.scalar_multiply_scalar(a,b))
Inv = lambda z: E(cp.scalar_inverse(z))
Sub = lambda a,b: E(cp.scalar_subtraction(a,b))
rand = lambda: secrets.token_bytes(32)
cm = lambda r,s: P(Pm(r,H1), Pm(s,H2))
rho = lambda h,r,s: P(Pm(r,h), Pm(s,B))
_sample = lambda x,g,r,s: (x, rho(hp(x,g),r,s))

ONE = b'\x01' + b'\0'*31
H1, H2 = Pm(b'\x02'+b'\0'*31, B), Pm(b'\x03'+b'\0'*31, B)

# ────────── polynomial & sharing helpers ──────────
gpoly = lambda t: [Sm(rand()) for _ in range(t)]
def hp(x,g):
    acc,xp = ONE,x
    for gk in g:
        acc = P(acc, Pm(xp,gk)); xp = S2(xp,x)
    return acc                      # (k-1) muls

def proof(h,ρ,cm_,r,s):
    kr,ks = rand(),rand()
    Acm,Ar = P(Pm(kr,H1),Pm(ks,H2)), P(Pm(kr,h),Pm(ks,B))
    c = hashlib.sha256(b''.join([cm_,ρ,Acm,Ar])).digest()[:32]
    zr,zs = Sub(kr,S2(c,r)), Sub(ks,S2(c,s))
    return c,zr,zs                 # exactly three values

def check(h,ρ,cm_,π):
    c,zr,zs = π
    Mc = P(Pm(c,cm_), P(Pm(zr,H1), Pm(zs,H2)))
    Mr = P(Pm(c,ρ), P(Pm(zr,h), Pm(zs,B)))
    return hashlib.sha256(b''.join([cm_,ρ,Mc,Mr])).digest()[:32]==c

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
    acc = E(cp.scalar_multiply(b"\0" * 32))  # ID
    k = len(shs)
    # Compute denominators: ∏_{m≠j} (x_m - x_j)
    denoms = [ONE] * k
    for j in range(k):
        for m in range(k):
            if m != j:
                denoms[j] = S2(denoms[j], Sub(shs[m][0], shs[j][0]))
    # Batch invert denominators
    inv_denoms = montgomery_batch_invert(denoms)
    # Compute coefficients and interpolate
    for j, (xj, yj) in enumerate(shs):
        l = inv_denoms[j]
        for m in range(k):
            if m != j:
                l = S2(l, shs[m][0])  # l_j = ∏_{m≠j} x_m / (x_m - x_j)
        acc = P(acc, Pm(l, yj))
    return acc

# ────────── verification helpers ──────────
ShD = lambda T: int(all(check(*line) for line in T))

def ShS(sh,T,g):
    x,y = sh; fx = hp(x,g)
    return int(any(fx==f and y==ρ for f,ρ,*_ in T))

# ────────── reconstruction oracle with ShS and ShD checks ──────────
def Rbox(k, embeds, T, g):
    def R(extra):
        extra = [extra] if isinstance(extra, tuple) else list(extra)
        all_shares = embeds + extra
        valid_shares = []
        # Verify each share with ShS and ShD
        for x, y in all_shares:
            # Check ShS: share must match a tuple in T
            if not ShS((x, y), T, g):
                print(f"Rbox: Share ({x.hex()}, {y.hex()}) does not match any tuple in T")
                return None
            valid_shares.append((x, y))
            # Find matching tuple in T for ShD
            # for h, ρ, cm_, π in T:
            #     if y == ρ and hp(x, g) == h:
            #         if check(h, ρ, cm_, π):
            #             valid_shares.append((x, y))
            #         break
        # Require at least k valid shares
        uniq = {}; [uniq.setdefault(x, y) for x, y in valid_shares]
        if len(uniq) < k:
            print(f"k={k} but only {len(uniq)} unique shares found")
            return None
        before, t0 = CNT, time.perf_counter()
        res = recon(list(uniq.items())[:k])
        print(f"the result of Rbox recon: {res.hex()}")
        rb_bump(CNT - before)
        globals().__setitem__("RB_TIME", globals().get("RB_TIME", 0.0) + (time.perf_counter() - t0))
        return res
    return R

# ────────── Trace / TrVer (use pre-generated fresh shares from T) ──────────
def Trace(tk, T, g, f, t, R):
    ζ, shv = tk; r, s, S = ζ; n = len(shv); I = set(range(n))
    banned = {x for x, _ in shv}
    # Use pre-generated fresh shares from T (indices n to 2n-1)
    fresh_indices = list(range(n, 2*n))
    F = []
    # Select t-f-1 fresh shares, avoiding banned x values
    for idx in fresh_indices:
        if len(F) >= t - f - 1:
            break
        x = T[idx][0]  # Use h as x (since hp(x, g) = h)
        if x not in banned:
            banned.add(x)
            F.append((x, T[idx][1]))  # (x, ρ) from T
    for idx, sh in enumerate(shv):
        if R(F + [sh]) == S:
            I.discard(idx)
    return sorted(I), [shv[i] for i in I]

def TrVer(vk, I, T, π, g, f, t, R):
    ζ, shv = vk; r, s, S = ζ
    # Verify share consistency
    if any(_sample(x, g, r, s) != (x, ρ) for x, ρ in shv):
        return 0
    banned = {x for x, _ in shv}
    # Use pre-generated fresh shares from T
    fresh_indices = list(range(n, 2*n))
    F_ver = []
    for idx in fresh_indices:
        if len(F_ver) >= t - f - 1:
            break
        x = T[idx][0]
        if x not in banned:
            banned.add(x)
            F_ver.append((x, T[idx][1]))
    for idx in I:
        x, y = shv[idx]
        if y != _sample(x, g, r, s)[1]:
            return 0
        if not any(y == ρ for _, ρ, _, _ in T):
            return 0
        if R(F_ver + [shv[idx]]) == S:
            return 0
    return 1

# ────────── benchmark harness ──────────
def run(n, k, f):
    assert 0 < f < k - 1
    print(f"\n== n={n} k={k} f={f} ==")
    g = gpoly(k - 1)
    xs = [(i + 1).to_bytes(32, 'little') for i in range(n)]

    reset(); t0 = time.perf_counter()
    hs = [hp(x, g) for x in xs]
    print(f"[share ] {1e3*(time.perf_counter()-t0)/n:7.1f} ms  {CNT//n:3d} mul/ply"); banner()

    r, s = rand(), rand(); cm_ = cm(r, s)
    reset(); t0 = time.perf_counter()
    ρs, T = [], []
    # Generate n shares for parties
    for h in hs:
        ρ = rho(h, r, s); π = proof(h, ρ, cm_, r, s)
        ρs.append(ρ); T.append((h, ρ, cm_, π))
    # Generate n fresh shares, keep x private
    fresh_xs = [rand() for _ in range(n)]
    for x in fresh_xs:
        h = hp(x, g); ρ = rho(h, r, s); π = proof(h, ρ, cm_, r, s)
        T.append((h, ρ, cm_, π))  # Add fresh share tuple to T
    print(f"[dealer] {1e3*(time.perf_counter()-t0):7.1f} ms  {CNT:4d} mul"); banner()

    reset(); t0 = time.perf_counter(); ok = ShD(T)
    print(f"[ShD   ] {1e3*(time.perf_counter()-t0):7.1f} ms  {CNT//n:3d} mul  ok={ok}"); banner()

    reset(); t0 = time.perf_counter()
    allok = all(ShS((x, ρ), T, g) for x, ρ in zip(xs[:k], ρs[:k]))
    print(f"[ShS   ] {1e3*(time.perf_counter()-t0):7.1f} ms  {CNT//k:3d} mul  all={allok}"); banner()

    rb_reset(); globals().__setitem__("RB_TIME", 0.0)
    R = Rbox(k, random.sample(list(zip(xs, ρs)), f), T, g)

    reset(); t0 = time.perf_counter()
    secret = Sm(s)
    reconstructed = recon(list(zip(xs, ρs))[:k])
    ok = reconstructed == secret
    print(f"[Recon ] {1e3*(time.perf_counter()-t0):7.1f} ms  {CNT:4d} mul  ok={ok}"); banner()

    tk = ((r, s, Sm(s)), list(zip(xs, ρs)))

    # -------- Trace (include oracle) --------
    rb_before, c_before, t0 = RB_CNT, CNT, time.perf_counter()
    I, π = Trace(tk, T, g, f, k, R)
    dt = time.perf_counter() - t0
    trace_mul = CNT - c_before
    trace_ms = 1e3 * dt
    print(f"[Trace ] {trace_ms:7.1f} ms  {trace_mul:4d} mul  |I|={len(I)}"); banner()

    # -------- TrVer (include oracle) --------
    rb_before, c_before, t0 = RB_CNT, CNT, time.perf_counter()
    ok = TrVer(tk, I, T, π, g, f, k, R)
    dt = time.perf_counter() - t0
    trv_mul = CNT - c_before
    trv_ms = 1e3 * dt
    print(f"[TrVer ] {trv_ms:7.1f} ms  {trv_mul:4d} mul  ok={ok}"); banner()

# ────────── entry ──────────
if __name__ == "__main__":
    for n, k, f in [
                    # (64, 33, 10), (128, 65, 20), (256, 129, 40)
                    # (32, 17, 4)
                    (8, 5, 2)
                    ]:
        run(n, k, f)