import random
from verification import ShS
from poly_helpers import hp, check
from curve_ops import rho, Sm, cm, B
from utils import reset, rb_bump, CNT
from reconstruction import recon
import time

def Trace(tk, T, g, f, t, R, cm):
    ζ, dsh, shv = tk
    r, s, S = ζ
    n = len(shv)
    I = set(range(n))
    # Select t-f-1 dummy shares
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
    ζ, dsh, shv = vk
    r, s, S = ζ
    # Check zeta and cm consistency
    if any(t[7] != cm for t in T):  # Use cm field (index 7)
        return 0
    # Check sh_i in shv
    if not all(sh in shv for sh in π):
        return 0
    # Verify px_i, pz_i consistency
    for x, px in shv:
        fx = hp(x, g)
        if px != rho(fx, r, s):
            return 0
    for z, pz in dsh:
        fz = hp(z, g)
        if pz != rho(fz, r, s):
            return 0
    # Select t-f-1 dummy shares
    DSH = []
    banned = {x for x, _ in shv}
    for z, pz in random.sample(dsh, t - f - 1):
        if z not in banned:
            DSH.append((z, pz))
            banned.add(z)
    # Check faulty shares
    for idx in I:
        x, y = shv[idx]
        if not any((y == px_j and hp(x, g) == fx_j) or (y == pz_j and hp(x, g) == fz_j) for fx_j, fz_j, px_j, pz_j, _, _, _, _ in T):
            return 0
        if R(DSH + [shv[idx]]) == S:
            return 0
    return 1

def Rbox(k, embeds, T, g, cm):
    def R(shares):
        shares = [shares] if isinstance(shares, tuple) else list(shares)
        all_shares = shares
        valid_shares = [] + embeds
        for x, y in all_shares:
            if not ShS((x, y), T, g):
                return None
            valid_shares.append((x, y))
        uniq = {}
        [uniq.setdefault(x, y) for x, y in valid_shares]
        if len(uniq) < k:
            return None
        before, t0 = CNT, time.perf_counter()
        res = recon(list(uniq.items())[:k])
        rb_bump(CNT - before)
        globals().__setitem__("RB_TIME", globals().get("RB_TIME", 0.0) + (time.perf_counter() - t0))
        return res
    return R