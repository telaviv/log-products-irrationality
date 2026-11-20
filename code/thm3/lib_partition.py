"""lib_partition.py -- CDT Theorem 6.0.2 with a nontrivial partition and dilated maps phi_k(z) = phi(r_k z), 0 < r_0 < ... < r_l = 1:
numerator = E[max(X, Y)], X, Y iid from the mixture sum_k ((theta_{k+1}-theta_k)/m) * law(log|phi_k| on T)
            + (1/m) sum_{k=1}^{l} theta_k^2 log(r_k / r_{k-1}),  denominator log|phi'(0)| - tau.  [Numerical]"""
import numpy as np, math
from scipy.optimize import minimize
from lib_map import *

def samples(c, r, N):
    th = 2 * np.pi * (np.arange(N) + 0.5) / N
    return np.real(logh_and_dlog(Map(c).psi(r * np.exp(1j * th)))[0])

def emax_mixture(vals, wts):
    v = np.concatenate(vals); w = np.concatenate([np.full(len(x), wt / len(x)) for x, wt in zip(vals, wts)])
    o = np.argsort(v); v = v[o]; w = w[o]
    W = np.cumsum(w); Wm = W - w
    return float(np.sum(v * (W * W - Wm * Wm)))

def numerator(c, rs, thetas, m=17, N=8192, cache=None):
    """rs: r_0 < ... < r_l = 1; thetas: theta_1 < ... < theta_l in (0, m)."""
    th = np.concatenate([[0.0], thetas, [float(m)]])
    vals = [samples(c, r, N) if cache is None else cache(r) for r in rs]
    wts = np.diff(th) / m
    pen = sum(th[k] ** 2 * math.log(rs[k] / rs[k - 1]) for k in range(1, len(rs))) / m
    return emax_mixture(vals, wts) + pen

def best_l1(c, m=17, N=8192):
    best = (numerator(c, [1.0], [], m, N), None)
    for r in (0.5, 0.7, 0.8, 0.85, 0.9, 0.95):
        for t in np.linspace(0.5, 12, 24):
            v = numerator(c, [r, 1.0], [t], m, N)
            if v < best[0]:
                best = (v, (r, t))
    return best

def optimise_partition(c, l, m=17, N=8192, x0=None):
    """minimise the numerator over r_0 < .. < r_{l-1} < 1 and theta_1 < .. < theta_l (parametrised monotonically)."""
    def unpack(x):
        u = np.exp(-np.cumsum(np.exp(x[:l]))[::-1])     # r_0 < ... < r_{l-1} < 1
        rs = list(u) + [1.0]
        t = np.cumsum(np.exp(x[l:]))
        return rs, t
    def f(x):
        rs, t = unpack(x)
        if t[-1] >= m: return 1e9
        return numerator(c, rs, t, m, N)
    if x0 is None:
        x0 = np.concatenate([np.full(l, math.log(0.7 / l)), np.full(l, math.log(2.0 / l))])
    r = minimize(f, x0, method='Nelder-Mead', options=dict(maxiter=600, xatol=1e-4, fatol=1e-6))
    rs, t = unpack(r.x)
    return r.fun, rs, list(t), r.x
