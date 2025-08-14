#!/usr/bin/env python3
"""
Tiny-TSS-PV v13 ▸ benchmark with time / muls / heap (tracemalloc) usage
────────────────────────────────────────────────────────────────────
Implements TSS-PV Ver 1.0 with ElGamal encryption for dummy shares.
Fixes TrVer and ShD to correctly check cm consistency using T[i][7].
"""

import time, random
from utils import reset, rb_reset, banner, rand, CNT, RB_CNT
from curve_ops import rho, Sm, cm, B
from poly_helpers import gpoly, hp, proof
from elgamal import elgamal_encrypt, Pm
from verification import ShD, ShS
from trace import Trace, TrVer, Rbox
from reconstruction import recon

def run(n, k, f):
    assert 0 < f < k - 1
    print(f"\n== n={n} k={k} f={f} ==")
    g = gpoly(k - 1)
    xs = [(i + 1).to_bytes(32, 'little') for i in range(n)]
    zs = [rand() for _ in range(n)]  # Dummy share x values

    # Player share generation
    reset()
    t0 = time.perf_counter()
    fx = [hp(x, g) for x in xs]
    fz = [hp(z, g) for z in zs]
    print(f"[share ] {1e3*(time.perf_counter()-t0)/n:7.1f} ms {CNT//n:3d} mul/ply")
    banner()

    # Dealer operations
    r, s = rand(), rand()
    S = Sm(s)
    cm_ = cm(r, s)

    # Generate ElGamal key pair for dealer
    sk = rand()
    pk = Pm(sk, B)
    reset()
    t0 = time.perf_counter()

    px, pz, T = [], [], []
    dsh = []
    for i in range(n):
        # Party share
        px_i = rho(fx[i], r, s)  # fx_i^r * S
        π_i = proof(fx[i], px_i, cm_, r, s)
        # Dummy share
        pz_i = rho(fz[i], r, s)  # fz_i^r * S
        ct_i = elgamal_encrypt(fz[i], pk, sk)  # (c1, c2, sk)
        πp_i = proof(fz[i], pz_i, cm_, r, s)
        T.append((fx[i], fz[i], px_i, pz_i, π_i, πp_i, ct_i[:2], cm_))  # Include cm_
        px.append(px_i)
        pz.append(pz_i)
        dsh.append((zs[i], pz_i))
        
    # Ensure shv uses exact px from T
    shv = [(xs[i], T[i][2]) for i in range(n)]  # Explicitly take px_i from T
    print(f"[dealer] {1e3*(time.perf_counter()-t0):7.1f} ms {CNT:4d} mul")
    banner()

    # VerifySD
    reset()
    t0 = time.perf_counter()
    ok = ShD(T, cm_)
    print(f"[ShD ] {1e3*(time.perf_counter()-t0):7.1f} ms {CNT//n:3d} mul ok={ok}")
    banner()

    # VerifySS
    reset()
    t0 = time.perf_counter()
    allok = all(ShS(sh, T, g) for sh in shv[:k])
    print(f"[ShS ] {1e3*(time.perf_counter()-t0):7.1f} ms {CNT//k:3d} mul all={allok}")
    banner()

    # Reconstruction
    rb_reset()
    R = Rbox(k, random.sample(shv, f), T, g, cm_)  # Use dummy shares for embeds
    reset()
    t0 = time.perf_counter()
    secret = S
    reconstructed = recon(shv[:k])
    ok = reconstructed == secret
    print(f"[Recon ] {1e3*(time.perf_counter()-t0):7.1f} ms {CNT:4d} mul ok={ok}")
    banner()

    # Trace
    tk = ((r, s, S), dsh, shv)
    rb_before, c_before, t0 = RB_CNT, CNT, time.perf_counter()
    I, π = Trace(tk, T, g, f, k, R, cm_)
    dt = time.perf_counter() - t0
    trace_mul = CNT - c_before
    trace_ms = 1e3 * dt
    print(f"[Trace ] {trace_ms:7.1f} ms {trace_mul:4d} mul |I|={len(I)}")
    banner()

    # TrVer
    rb_before, c_before, t0 = RB_CNT, CNT, time.perf_counter()
    ok = TrVer(tk, I, T, π, g, f, k, R, cm_)
    dt = time.perf_counter() - t0
    trv_mul = CNT - c_before
    trv_ms = 1e3 * dt
    print(f"[TrVer ] {trv_ms:7.1f} ms {trv_mul:4d} mul ok={ok}")
    banner()

if __name__ == "__main__":
    for n, k, f in [
        (32, 17, 11),
        (64, 33, 22), (128, 65, 43),
        (256, 129, 86)
        # (8, 5, 2)
        ]:
        run(n, k, f)
