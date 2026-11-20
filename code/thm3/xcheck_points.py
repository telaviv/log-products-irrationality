"""xcheck_points.py map.json N -- independent cross-check of the range certificate: for the specific pair (n, m) = (N, N-1)
(|m/n - 1| = 1/N), compute the Gamma_0(2)-orbits of its singular values (lib_obs, mp 80 digits) and check winding numbers of
psi(T) about every obstruction with |q| < 0.95 (must be 0) and about q_alpha (must be 1).  [Numerical]"""
import sys, json, numpy as np
from lib_map import *; from lib_obs import *
c = np.array(json.load(open(sys.argv[1]))['c']); N = int(sys.argv[2])
O, rec = get_obs(N); gen, qa0, qan = split(O)
E = np.array(gen + qan)
wn, mx = winding(c, E, 2 ** 17)
w0, _ = winding(c, [qa0, 0j], 2 ** 17)
d = min_dist(c, E[np.abs(E) < mx + 0.05], 2 ** 17)
print(dict(N=N, n_obs=len(E), nonzero_windings=int(np.sum(wn != 0)), winding_alpha_and_0=w0.tolist(), max_psi=mx, min_dist=float(d.min())))
