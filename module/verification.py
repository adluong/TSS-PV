from poly_helpers import hp, check
def ShD(T, cm):
    # Check cm consistency
    if any(t[7] != cm for t in T):  # Use cm field (index 7)
        return 0
    # Verify proofs for px_i and pz_i
    for fx_i, fz_i, px_i, pz_i, π_i, πp_i, ct_i, _ in T:
        if not (check(fx_i, px_i, cm, π_i) and check(fz_i, pz_i, cm, πp_i)):
            return 0
    return 1

def ShS(sh, T, g):
    x, y = sh
    fx = hp(x, g)
    return int(any((fx == fx_j and y == px_j) or (fx == fz_j and y == pz_j) for fx_j, fz_j, px_j, pz_j, _, _, _, _ in T))