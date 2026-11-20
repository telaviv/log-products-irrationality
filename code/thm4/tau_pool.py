"""tau_pool.py -- exact tau(b;e) (CDT Def. 6.0.1, (6.0.4)-(6.0.6)) for the pooled list of 7 + 10K functions:
pure B1..B7 (sigma = 0,2,2,4,4,4,4; e = 0,0,1,0,0,1,1) shared, plus per pair G, G', G'', G''' (sigma 4, e 0) and six
integrals (sigma 4, e 1).  tau_sharp is evaluated at an explicit rational xi (any xi gives a valid upper bound for the
minimum in (6.0.6), hence a valid tau); xi is the best point of a grid scan.  Exact Fractions.  Output runs/tau_pool.json."""
import sys, json
from fractions import Fraction as Fr
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # logprod copy: f1_tau.py is in this folder
import f1_tau as FT

def lists(K):
    sig = [0, 2, 2, 4, 4, 4, 4] + [4] * (10 * K)
    e = [0, 0, 1, 0, 0, 1, 1] + ([0] * 4 + [1] * 6) * K
    return sig, e

out = {}
for K in (1, 2, 3, 5):
    sig, e = lists(K)
    m = len(sig)
    tf = FT.tau_flat([[s] for s in sig])
    (xi, ts), _ = FT.tau_sharp_scan(e, Fr(1, 100))
    # refine around xi on a finer grid
    best = (xi, ts)
    x = xi - Fr(1, 100)
    while x <= xi + Fr(1, 100):
        if x > 0:
            v = FT.tau_sharp_at(e, x)
            if v < best[1]:
                best = (x, v)
        x += Fr(1, 2000)
    xi, ts = best
    tau = tf + ts
    out[K] = dict(m=m, sum_e=sum(e), tau_flat=str(tf), xi=str(xi), tau_sharp=str(ts), tau=str(tau), tau_float=float(tau))
    print(K, out[K])
assert out[1]['tau'] == '1032659/242760', out[1]['tau']   # CDT Theorem C value (xi = 7/4)
json.dump(out, open('runs/tau_pool.json', 'w'), indent=1)
