"""lib_obs.py -- obstruction sets in the q-disc for the CDT Theorem C machine (Gamma_0(2) orbits of the singular values),
built with cdt4_budget.obstructions (read-only import). Cached to runs/obs_<eps>.json.  [Numerical, mp 80 digits for preimages]"""
import sys, os, json, math, cmath
B = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..', '..', 'schanuel', 'five_exp', 'cdt4'))
sys.path.insert(0, B)
import cdt4_budget as CB
from mpmath import mp

def ab_of(eps_n):
    """t = m/n = 1 - 1/N with N = eps_n (task normalisation (n, m) = (N, N-1)), CDT aC = 2N-1, bC = 2(N-1)-1."""
    N = eps_n
    return 2 * N - 1, 2 * (N - 1) - 1

def get_obs(N, qmax=0.95, cache=True):
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'runs', f'obs_N{N}_q{qmax}.json')
    if cache and os.path.exists(fn):
        d = json.load(open(fn))
        return {k: [complex(*p) for p in v] for k, v in d['q'].items()}, d
    aC, bC = ab_of(N)
    mp.dps = 80
    CB.QMAX = qmax
    O = {}
    S = CB.F3.sing_y(aC, bC)
    rec = dict(N=N, aC=aC, bC=bC, qmax=qmax, y={}, check={})
    for k, v in S.items():
        t0 = CB.one_preimage(v['y'])
        qs, _ = CB.orbit_q(t0, qmax)
        O[k] = qs
        rec['y'][k] = str(v['y'])
        rec['check'][k] = max(abs(complex(CB.F3.h_mp(q)) - complex(v['y'])) / max(1.0, abs(complex(v['y']))) for q in qs)
    rec['q'] = {k: [(q.real, q.imag) for q in v] for k, v in O.items()}
    json.dump(rec, open(fn, 'w'), indent=1)
    return O, rec

def split(O):
    qa = O['a-b-']
    genuine = [q for k in ('a+b-', 'a-b+', 'a+b+') for q in O[k]]
    return genuine, qa[0], qa[1:]
