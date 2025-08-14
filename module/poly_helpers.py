import hashlib
from curve_ops import ONE, Pm, P, S2, Inv, Sub, Sm, rand, H1, H2, B
from utils import E, bump

def gpoly(t):
    return [Sm(rand()) for _ in range(t)]

def hp(x, g):
    acc, xp = ONE, x
    for gk in g:
        acc = P(acc, Pm(xp, gk))
        xp = S2(xp, x)
    return acc  # (k-1) muls

def proof(h, p, cm_, r, s):
    kr, ks = rand(), rand()
    Acm = P(Pm(kr, H1), Pm(ks, H2))
    Ar = P(Pm(kr, h), Pm(ks, B))
    c = hashlib.sha256(b''.join([cm_, p, Acm, Ar])).digest()[:32]
    zr = Sub(kr, S2(c, r))
    zs = Sub(ks, S2(c, s))
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

