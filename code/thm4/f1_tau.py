"""f1_tau.py -- exact reproduction of the CDT denominator rates tau(b;e) (arXiv:2408.15403v2, Def. 6.0.1, (6.0.4)-(6.0.6)).

Checks the two published instances:
  Theorem A (14 functions, Sec. 13):  tau_flat = 191/49,  tau_sharp = 27/80 at xi = 2,     tau = 16603/3920
  Theorem C (17 functions, Sec. 14.5): tau_flat = 1136/289, tau_sharp = 78419/242760 at xi = 7/4, tau = 1032659/242760
and scans xi on a rational grid to confirm the minimiser.  Exact rationals throughout (fractions.Fraction).
Output: runs/f1_tau.json
"""
from fractions import Fraction as Fr
import json, math, sys, os

def floor_fr(x):
    return x.numerator // x.denominator

def harmonic(k):
    return sum(Fr(1, h) for h in range(1, k + 1)) if k >= 1 else Fr(0)

def I_uvw(u, v, w):
    """Definition 6.0.1: I_u^v(w) for 0 <= max(u,1) <= v, w <= v. Exact."""
    u, v, w = Fr(u), Fr(v), Fr(w)
    total = Fr(0)
    # term 1: int_{min(u,1)}^{1} max(t-w,0) dt
    lo = min(u, 1)
    if lo < 1:
        a = max(lo, w)   # max(t-w,0) > 0 for t > w
        if a < 1:
            total += Fr(1, 2) * ((1 - w) ** 2 - (a - w) ** 2)
    M = max(Fr(1), w)
    lo2 = max(u, 1)
    # breakpoints of floor((t-1)/M) and floor((t+max(0,w-1))/M) on [lo2, v]
    pts = {lo2, v}
    k = 0
    while True:
        t1 = 1 + k * M
        t2 = k * M - max(Fr(0), w - 1)
        if t1 > v and t2 > v:
            break
        for t in (t1, t2):
            if lo2 < t < v:
                pts.add(t)
        k += 1
    pts = sorted(pts)
    for a, b in zip(pts[:-1], pts[1:]):
        mid = (a + b) / 2
        k2 = floor_fr((mid - 1) / M)
        k3 = floor_fr((mid + max(Fr(0), w - 1)) / M)
        assert k3 >= 1
        total += harmonic(k2) * (b - a)
        # int_a^b max(t/k3 - w, 0) dt ; positive for t > k3*w
        c = max(a, k3 * w)
        if c < b:
            total += (b * b - c * c) / (2 * k3) - w * (b - c)
    return total

def tau_flat(b_rows):
    m = len(b_rows)
    sig = sorted(sum(r) for r in b_rows)
    return Fr(sum((2 * i - 1) * s for i, s in enumerate(sig, 1)), m * m)

def tau_sharp_at(e, xi):
    m = len(e)
    return Fr(2, m * m) * (xi * sum(e) + max(e) * I_uvw(xi, m, xi))

def tau_sharp_scan(e, step=Fr(1, 200)):
    m = len(e)
    best = None
    grid = []
    xi = Fr(0)
    while xi <= m:
        val = tau_sharp_at(e, xi)
        grid.append((xi, val))
        if best is None or val < best[1]:
            best = (xi, val)
        xi += step
    return best, grid

def instance(name, nb, cols, e, xi_paper, expect):
    # b array: column j has zeros in the first u_j rows then 2's (CDT (6.0.3)); cols = list of u_j
    rows = [[0 if i < u else 2 for u in cols] for i in range(nb)]
    tf = tau_flat(rows)
    ts = tau_sharp_at(e, Fr(xi_paper))
    best, grid = tau_sharp_scan(e)
    out = dict(name=name, m=nb, tau_flat=str(tf), tau_sharp_at_paper_xi=str(ts), xi_paper=str(Fr(xi_paper)),
               tau=str(tf + ts), tau_float=float(tf + ts), scan_best_xi=str(best[0]), scan_best=str(best[1]),
               expect=expect, ok=(str(tf) == expect['tau_flat'] and str(ts) == expect['tau_sharp']))
    # minimising interval (within 1e-12 of best) on the grid
    flat = [str(x) for x, v in grid if v == best[1]]
    out['minimiser_grid_points'] = flat[:3] + (['...'] if len(flat) > 6 else []) + flat[-3:] if len(flat) > 6 else flat
    print(f"{name}: tau_flat={tf} ({float(tf):.6f}) tau_sharp(xi={Fr(xi_paper)})={ts} ({float(ts):.6f}) "
          f"tau={tf+ts} ({float(tf+ts):.6f}); scan min {best[1]} at xi={best[0]}; expected {expect}; OK={out['ok']}")
    return out

if __name__ == '__main__':
    res = {}
    # Theorem A: 14 functions, (13.0.2): col1 u=1, col2 u=3; e=(0,0,1;0,0,0,0,0,0;1,1,1,1,1)
    res['thmA_14'] = instance('ThmA_14', 14, [1, 3], [0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1], 2,
                              dict(tau_flat='191/49', tau_sharp='27/80', tau='16603/3920'))
    # Theorem C: 17 functions, Sec 14.5: same column shape; e=(0,0,1;0*6;1*8)
    res['thmC_17'] = instance('ThmC_17', 17, [1, 3], [0, 0, 1] + [0] * 6 + [1] * 8, Fr(7, 4),
                              dict(tau_flat='1136/289', tau_sharp='78419/242760', tau='1032659/242760'))
    # sanity: I_2^14(2) should give (12 + I)/98 = 27/80  ->  I = 27*98/80 - 12
    I = I_uvw(2, 14, 2)
    print('I_2^14(2) =', I, ' check 27/80 ->', Fr(27, 80) * 98 - 12)
    I17 = I_uvw(Fr(7, 4), 17, Fr(7, 4))
    print('I_{7/4}^17(7/4) =', I17, ' check ->', Fr(78419, 242760) * 289 / 2 - 9 * Fr(7, 4))
    os.makedirs('runs', exist_ok=True)
    json.dump(res, open('runs/f1_tau.json', 'w'), indent=1)
