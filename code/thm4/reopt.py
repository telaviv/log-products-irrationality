"""reopt.py eps m tau start.json tag -- re-optimise a range map for the POOLED objective L - T/m - tau (m = 7 + 10K) instead of
CDT's m = 17, with the ray endpoint at Y1 = 1/eps - 12 (6 units of safety below Lemma S's 1/eps - 6; the range pilot's optimiser
used 1/eps, which let the 10^-5 map cross the real axis at 100000.06).  Uses range_pilot's lib_range/lib_opt unchanged.
Output runs/reopt_<tag>.json.  [Numerical]; certification is cert_pool.py."""
import sys, os, json, math, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "thm3"))  # logprod copy: was an absolute path to range_pilot
import numpy as np
from fractions import Fraction
from lib_map import *; from lib_range import *; from lib_opt import lbfgs
eps = float(sys.argv[1]); m = int(sys.argv[2]); tau = float(Fraction(sys.argv[3])); src = sys.argv[4]; tag = sys.argv[5]
c = np.array(json.load(open(src))['c'], dtype=float)
logY1 = math.log(1 / eps - 12)
def prob(mu):
    P = RangeProblem(eps, mu=mu, m=m, tau=tau)
    P.logY1 = logY1
    return P
P = prob(1e-3)
rho = 1.0
while not P.admissible(c + np.eye(len(c))[0] * math.log(rho))[0]:
    rho *= 0.999
    if rho < 0.3: raise SystemExit('no admissible rescaling')
c = c + np.eye(len(c))[0] * math.log(rho)
t0 = time.time()
log = open(f'runs/reopt_{tag}.log', 'w')
print(f'start rho={rho:.4f} L={LOG256 + c[0]:.5f}', flush=True)
for mu in (1e-3, 3e-4, 1e-4, 3e-5, 1e-5):
    c = lbfgs(prob(mu), c, iters=200, log=log, tag=f'eps={eps:.1e} m={m} mu={mu:.0e}')
L = LOG256 + c[0]
Ts = {NN: T_and_grad(c, NN, False)[0] for NN in (1024, 4096)}
res = dict(eps=eps, m=m, tau=tau, logY1_used=logY1, L=L, T=Ts, M={k: L - v / m - tau for k, v in Ts.items()},
           admissible=prob(1e-3).admissible(c), c=list(c), secs=time.time() - t0, src=src)
print(json.dumps({k: v for k, v in res.items() if k != 'c'}), flush=True)
json.dump(res, open(f'runs/reopt_{tag}.json', 'w'), indent=1)
