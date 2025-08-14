#!/usr/bin/env python3
"""
TSS-PV (no-mask) demo on Curve25519 / Ed25519

Adds tiny wrappers  P_add, P_mul, … so that every result coming
from curve25519_python is passed through  ensure_bytes(..) exactly once.
This prevents  TypeError: list → bytes  while keeping the code fast.
"""
import secrets, hashlib, time
import curve25519_python as cp          # low-level EC helpers

# ---------------------------------------------------------------------------#
# Small helper ensures every value is bytes                                   #
# ---------------------------------------------------------------------------#
def ensure_bytes(x):
    return bytes(x) if isinstance(x, list) else x   # no copy if already bytes

# one-line wrappers around cp-calls
P_add  = lambda a, b: ensure_bytes(cp.point_addition(a, b))
P_mul  = lambda s, P: ensure_bytes(cp.point_multiply(s, P))
S_mul  = lambda s:    ensure_bytes(cp.scalar_multiply(s))
S_mul2 = lambda a, b: ensure_bytes(cp.scalar_multiply_scalar(a, b))
S_sub  = lambda a, b: ensure_bytes(cp.scalar_subtraction(a, b))
S_inv  = lambda x:    ensure_bytes(cp.scalar_inverse(x))

# ---------------------------------------------------------------------------#
# Constants                                                                   #
# ---------------------------------------------------------------------------#
n, k  = 32, 22
B     = ensure_bytes(                      # Ed25519 base point (compressed y)
        b'\x58' + b'\x66'*31)
IDENT = b'\x01' + b'\x00'*31

# independent public bases  H1 = 2·B,  H2 = 3·B
H1 = P_mul(b'\x02'+b'\x00'*31, B)
H2 = P_mul(b'\x03'+b'\x00'*31, B)

rand_scalar = lambda: secrets.token_bytes(32)

# ---------------------------------------------------------------------------#
# Polynomial helpers                                                          #
# ---------------------------------------------------------------------------#
def gen_glist(t):                        # random curve points g_k
    return [S_mul(rand_scalar()) for _ in range(t)]

def hp_value(x, glist):                  # H(x) = Σ g_k^{x^k}
    acc, x_pow = IDENT, x
    for g_k in glist:
        acc   = P_add(acc, P_mul(x_pow, g_k))
        x_pow = S_mul2(x_pow, x)
    return acc

def masked_shares(h_list, w_list, r, s):
    return [P_add(P_mul(r, h), P_mul(s, w))
            for h, w in zip(h_list, w_list)]

# ---------------------------------------------------------------------------#
# Dealer commitment and FS proof                                              #
# ---------------------------------------------------------------------------#
def dealer_commit(r, s):
    return P_add(P_mul(r, H1), P_mul(s, H2))

def fs_proof(h_list, w_list, rho_list, cm, r, s):
    k_r, k_s  = rand_scalar(), rand_scalar()
    A_cm      = P_add(P_mul(k_r, H1), P_mul(k_s, H2))
    A_rho     = [P_add(P_mul(k_r, h), P_mul(k_s, w))
                 for h, w in zip(h_list, w_list)]
    c = hashlib.sha256(b"".join([cm] + rho_list + [A_cm] + A_rho)).digest()[:32]
    z_r = S_sub(k_r, S_mul2(c, r))
    z_s = S_sub(k_s, S_mul2(c, s))
    return c, z_r, z_s

def fs_verify(h_list, w_list, rho_list, cm, proof):
    c, z_r, z_s = proof
    M_cm  = P_add(P_mul(c, cm), P_add(P_mul(z_r, H1), P_mul(z_s, H2)))
    M_rho = [P_add(P_mul(c, rho),
                   P_add(P_mul(z_r, h), P_mul(z_s, w)))
             for rho, h, w in zip(rho_list, h_list, w_list)]
    c2 = hashlib.sha256(b"".join([cm] + rho_list + [M_cm] + M_rho)).digest()[:32]
    return c2 == c

# ---------------------------------------------------------------------------#
# Group Lagrange interpolation                                                #
# ---------------------------------------------------------------------------#
def reconstruct(shares):
    res = IDENT
    for j, (xj, yj) in enumerate(shares):
        Δ = IDENT
        for m, (xm, _) in enumerate(shares):
            if m == j: continue
            Δ = S_mul2(Δ,
                       S_mul2(xm,
                              S_inv(S_sub(xm, xj))))
        res = P_add(res, P_mul(Δ, yj))
    return res

# ---------------------------------------------------------------------------#
# Demo                                                                        #
# ---------------------------------------------------------------------------#
def test_protocol_no_mask():
    print("─ TSS-PV demo  (n={}, k={})".format(n, k))
    glist  = gen_glist(k-1)
    x_vals = [int.to_bytes(i+1, 32, 'little') for i in range(n)]
    h_list = [hp_value(x, glist) for x in x_vals]
    w_list = [B]*n

    r, s = rand_scalar(), rand_scalar()
    cm   = dealer_commit(r, s)
    rho  = masked_shares(h_list, w_list, r, s)

    ok   = fs_verify(h_list, w_list, rho, cm,
                     fs_proof(h_list, w_list, rho, cm, r, s))
    rec  = reconstruct(list(zip(x_vals, rho))[:k])
    print("   • ZK proof ok? ", ok)
    print("   • Reconstruction ok? ", rec == S_mul(s))

if __name__ == "__main__":
    t0 = time.perf_counter()
    test_protocol_no_mask()
    print("Finished in {:.2f} ms".format(1e3*(time.perf_counter()-t0)))
