#!/usr/bin/env python3
"""
TSS-PV (no mask) – micro-benchmark
Counts EC multiplications so that
  • players:  k-1  muls  to compute  f(x_i)
  • dealer:   6·n muls  for shares + proofs  (2 ρ  + 4 proof each)
  • players:  4    muls  to verify share/proof
"""

import secrets, hashlib, time
import curve25519_python as cp

# ──────────────────────────────────────────────────────────────────────
# 1.  safe wrapper + multiplication counter
# ──────────────────────────────────────────────────────────────────────
MULTI_CNT = 0
def ensure_bytes(x): return bytes(x) if isinstance(x, list) else x
def _bump():          globals().__setitem__("MULTI_CNT", MULTI_CNT + 1)
def reset_cnt():      globals().__setitem__("MULTI_CNT", 0)
cnt = lambda: MULTI_CNT

# point / scalar wrappers
P_add = lambda A, B: ensure_bytes(cp.point_addition(A, B))
def P_mul(s, P):  _bump(); return ensure_bytes(cp.point_multiply(s, P))
def S_mul(s):     _bump(); return ensure_bytes(cp.scalar_multiply(s))
S_mul2 = lambda a, b: ensure_bytes(cp.scalar_multiply_scalar(a, b))
S_sub  = lambda a, b: ensure_bytes(cp.scalar_subtraction(a, b))
S_inv  = lambda x:    ensure_bytes(cp.scalar_inverse(x))
rand_scalar = lambda: secrets.token_bytes(32)

# ──────────────────────────────────────────────────────────────────────
# 2.  parameters & fixed bases
# ──────────────────────────────────────────────────────────────────────
def run_case(n, k):
    print(f"\n── case n={n:3d}  k={k:3d} ──")

    B     = ensure_bytes(b'\x58' + b'\x66'*31)   # Ed25519 generator
    IDENT = b'\x01' + b'\x00'*31
    H1, H2 = P_mul(b'\x02'+b'\x00'*31, B), P_mul(b'\x03'+b'\x00'*31, B)

    # polynomial helpers ------------------------------------------------
    def gen_glist(t): return [S_mul(rand_scalar()) for _ in range(t)]

    def hp_value(x, glist):
        acc, x_pow = IDENT, x
        for g_k in glist:
            acc   = P_add(acc, P_mul(x_pow, g_k))   # mul here
            x_pow = S_mul2(x_pow, x)
        return acc

    # dealer helpers ----------------------------------------------------
    dealer_cm  = lambda r,s: P_add(P_mul(r,H1), P_mul(s,H2))
    dealer_rho = lambda h,r,s: P_add(P_mul(r,h),  P_mul(s,B))

    def proof_single(h, rho, cm, r, s):
        k_r, k_s = rand_scalar(), rand_scalar()
        A_cm  = P_add(P_mul(k_r,H1), P_mul(k_s,H2))     # 2 muls
        A_rho = P_add(P_mul(k_r,h),  P_mul(k_s,B))      # 2 muls
        c     = hashlib.sha256(b"".join([cm,rho,A_cm,A_rho])).digest()[:32]
        z_r   = S_sub(k_r, S_mul2(c,r))
        z_s   = S_sub(k_s, S_mul2(c,s))
        return c,z_r,z_s                                # +0 muls

    def verify_single(h,rho,cm,prf):
        c,z_r,z_s = prf
        M_cm  = P_add(P_mul(c,cm), P_add(P_mul(z_r,H1), P_mul(z_s,H2)))
        M_rho = P_add(P_mul(c,rho),P_add(P_mul(z_r,h),  P_mul(z_s,B)))
        c2 = hashlib.sha256(b"".join([cm,rho,M_cm,M_rho])).digest()[:32]
        return c2==c                                    # 4 muls

    # generate public polynomial coefficients ---------------------------
    glist  = gen_glist(k-1)
    x_vals = [int.to_bytes(i+1, 32, 'little') for i in range(n)]

    # player cost: f(x_i) ----------------------------------------------
    reset_cnt(); t0=time.perf_counter()
    h_list=[hp_value(x,glist) for x in x_vals]
    t_players = 1e3*(time.perf_counter()-t0)
    print(f"[player] share   avg {t_players/n:5.2f} ms   {cnt()//n:3d} muls (k-1)")
    
    # dealer: shares + proofs ------------------------------------------
    r,s = rand_scalar(), rand_scalar()
    cm  = dealer_cm(r,s)           # 2 muls (not counted in 6n)  ← done *before* reset
    reset_cnt(); t0=time.perf_counter()
    rho_list=[]; proofs=[]
    for h in h_list:
        rho   = dealer_rho(h,r,s)              # 2 muls
        proof = proof_single(h,rho,cm,r,s)     # 4 muls
        rho_list.append(rho); proofs.append(proof)
    t_dealer=1e3*(time.perf_counter()-t0)
    print(f"[dealer] total        {t_dealer:5.2f} ms   {cnt():4d} muls (expect {6*n})")

    # verification cost per player -------------------------------------
    reset_cnt(); t0=time.perf_counter()
    oks=[verify_single(h,rho,cm,prf) for h,rho,prf in zip(h_list,rho_list,proofs)]
    t_ver = 1e3*(time.perf_counter()-t0)
    print(f"[player] verify  avg {t_ver/n:5.2f} ms   {cnt()//n:3d} muls (expect 4)")
    assert all(oks)

    # reconstruction ----------------------------------------------------
    def reconstruct(threshold_shares):
        res=IDENT
        for j,(xj,yj) in enumerate(threshold_shares):
            Δ=IDENT
            for m,(xm,_) in enumerate(threshold_shares):
                if m==j: continue
                Δ=S_mul2(Δ,S_mul2(xm,S_inv(S_sub(xm,xj))))
            res=P_add(res,P_mul(Δ,yj))
        return res
    rec = reconstruct(list(zip(x_vals,rho_list))[:k])
    print("Reconstruction OK? ", rec==S_mul(s))

# ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    total_t = time.perf_counter()
    run_case(32,22)
    run_case(64,44)
    run_case(128,88)
    run_case(256,176)
    print("\nFinished in {:.2f} ms".format(1e3*(time.perf_counter()-total_t)))
