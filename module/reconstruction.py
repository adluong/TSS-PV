import time
from utils import rb_bump, bump, RB_TIME
from poly_helpers import montgomery_batch_invert
from curve_ops import ONE, Pm, P, S2, Inv, Sub, Sm, rand, H1, H2, B

def recon(shs):
    acc = Sm(b"\0" * 32)  # ID
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
                l = S2(l, shs[m][0])  # l_j = ∏_{m≠j} x_m / (x_m - x_j)
        acc = P(acc, Pm(l, yj))
    return acc