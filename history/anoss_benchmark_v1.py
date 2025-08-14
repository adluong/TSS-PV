#!/usr/bin/env python3
"""Minimized TSS-PV benchmark for Curve25519 with ZKPs"""

import secrets, hashlib, time, tracemalloc, curve25519_python as cp

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
CNT = 0
bump = lambda: globals().__setitem__("CNT", CNT + 1)
reset = lambda: globals().__setitem__("CNT", 0)

# ───── constants (32 B) ─────
ID = b"\x01" + b"\0" * 31  # identity point
B = E(cp.scalar_multiply(b"\x01" + b"\0" * 31))  # base point

# ───── curve wrappers ─────
P = lambda A, B_: E(cp.point_addition(E(A), E(B_)))
Pm = lambda s, P_: bump() or E(cp.point_multiply(E(s), E(P_)))
Sm = lambda s: bump() or E(cp.scalar_multiply(E(s)))
S2 = lambda a, b: E(cp.scalar_multiply_scalar(E(a), E(b)))
Sub = lambda a, b: E(cp.scalar_subtraction(E(a), E(b)))
Inv = lambda z: E(cp.scalar_inverse(E(z)))
rand = lambda: secrets.token_bytes(32)

# ───── protocol functions ─────
def gen_glist(t):
    return [Sm(rand()) for _ in range(t)]

def compute_hp(x, m, g):
    h, xp = ID, x
    for gk in g:
        h = P(h, Pm(xp, gk))
        xp = S2(xp, x)
    return Pm(m, h)

def compute_w(m):
    return Sm(m)

def compute_fp(hlist, wlist, r, s):
    if len(hlist) != len(wlist):
        raise ValueError("hlist and wlist must have same length")
    return [P(Pm(r, h), Pm(s, w)) for h, w in zip(hlist, wlist)]

def reconstruct(shs):
    acc = ID
    for j, (xj, yj) in enumerate(shs):
        delta = b"\x01" + b"\0" * 31
        for m, (xm, _) in enumerate(shs):
            if m != j:
                delta = S2(delta, S2(xm, Inv(Sub(xm, xj))))
        acc = P(acc, Pm(delta, yj))
    return acc

def zkpok_proof(hlist, wlist, flist, r, s):
    y, z = rand(), rand()
    commits = [P(Pm(y, h), Pm(z, w)) for h, w in zip(hlist, wlist)]
    c = hashlib.sha256(b"".join(flist + commits)).digest()[:32]
    u = Sub(y, S2(c, r))
    v = Sub(z, S2(c, s))
    return c, u, v, flist

def zkpok_verify(hlist, wlist, proof):
    c, u, v, flist = proof
    M = [P(Pm(c, f), P(Pm(u, h), Pm(v, w))) for f, h, w in zip(flist, hlist, wlist)]
    return hashlib.sha256(b"".join(flist + M)).digest()[:32] == c

def player_zkpok_proof(w, m):
    k = rand()
    commit = Pm(k, B)
    c = hashlib.sha256(b"".join([w, commit])).digest()[:32]
    z = Sub(k, S2(c, m))
    return c, z

def player_zkpok_verify(w, proof):
    c, z = proof
    M = P(Pm(c, w), Pm(z, B))
    return hashlib.sha256(b"".join([w, M])).digest()[:32] == c

# ───── benchmark driver ─────
def bench(n=32, k=22):
    t = k - 1
    print(f"\n== n={n} k={k} t={t} ==")
    g = gen_glist(t)
    xs = [int.to_bytes(i + 1, 32, "little") for i in range(n)]

    reset(); t0 = time.perf_counter()
    ms = [rand() for _ in range(n)]
    m_inv = [Inv(m) for m in ms]
    hlist = [compute_hp(x, m, g) for x, m in zip(xs, ms)]
    wlist = [compute_w(m) for m in ms]
    pi_proofs = [player_zkpok_proof(w, m) for w, m in zip(wlist, ms)]
    print(f"[hp+w+pi] {1e3*(time.perf_counter()-t0)/n:7.1f} ms  {CNT//n:3d} mul/ply"); heap()

    r, s = rand(), rand()
    reset(); t0 = time.perf_counter()
    flist = compute_fp(hlist, wlist, r, s)
    proof = zkpok_proof(hlist, wlist, flist, r, s)
    pi_ok = all(player_zkpok_verify(w, pi_proof) for w, pi_proof in zip(wlist, pi_proofs))
    print(f"[fp+D+pi] {1e3*(time.perf_counter()-t0):7.1f} ms  {CNT:4d} mul  pi_ok={pi_ok}"); heap()

    shs = [(x, Pm(m_inv, f)) for x, m_inv, f in zip(xs, m_inv, flist)]
    reset(); t0 = time.perf_counter()
    S = reconstruct(shs[:k])
    ok = S == Sm(s)
    print(f"[Recon ] {1e3*(time.perf_counter()-t0):7.1f} ms  {CNT:4d} mul  ok={ok}"); heap()

# ───── main ─────
if __name__ == "__main__":
    for n, k in [(64, 43), (128, 86), (256, 172), (512, 344)]:
        bench(n, k)