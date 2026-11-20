"""lib_refine.py -- Bost-Charles characteristic T(r, phi) = iint log|phi(z) - phi(r w)| (CDT Def. 7.1.2) and the convexity
refinements: Thm 7.1.6 / (8.0.4) numerator T(1) - (1/m) sum (dT)^2/dlog r (slopes <= m); Thm 7.6.4 / (8.0.5) with the
s_h weights solving s_h = (kappa_{h+1} - kappa_h)/m.  [Numerical]"""
import numpy as np, math
from lib_map import *

def T_r(c, r, N):
    z, w = grids(N)
    m = Map(c)
    fz = h_of(m.psi(z)); fw = h_of(m.psi(r * w))
    tot = 0.0
    ch = max(1, 2_000_000 // N)
    for i in range(0, N, ch):
        df = fz[i:i + ch, None] - fw[None, :]
        tot += float(np.sum(np.log(np.abs(df))))
        if r == 1.0:
            tot -= float(np.sum(np.log(np.abs(z[i:i + ch, None] - w[None, :]))))
    return tot / (N * N)

def profile(c, N=2048, rs=None):
    if rs is None:
        rs = np.concatenate([np.geomspace(0.05, 0.5, 6), np.linspace(0.55, 0.95, 9), [0.97, 0.985, 1.0]])
    return np.array(rs), np.array([T_r(c, float(r), N) for r in rs])

def refined_716(rs, Ts, m=17):
    """numerator of (7.1.8): T(1) - (1/m) sum_k (T_k - T_{k-1})^2 / (log r_k - log r_{k-1}); requires all slopes <= m."""
    lr = np.log(rs)
    dT = np.diff(Ts); dl = np.diff(lr)
    kap = dT / dl
    return Ts[-1] - np.sum(dT ** 2 / dl) / m, kap

def refined_764(rs, Ts, m=17):
    """numerator of (8.0.5): T(1) - sum_{h<l} s_h (T(1) - T(r_h)), s_h = (kappa_{h+1}-kappa_h)/m, kappa_0 = 0, kappa_{l+1} = m.
    NOTE: in Thm 7.6.4 the kappa_k are themselves defined with the s-weighted line bundle (7.6.1), giving a linear system (7.6.3);
    here we use the first-order (s = (0..0,1)) slopes as an approximation -- used only for orientation, not for claims."""
    lr = np.log(rs)
    kap = np.diff(Ts) / np.diff(lr)
    kk = np.concatenate([[0.0], kap, [float(m)]])
    s = np.diff(kk) / m          # s_0..s_l
    return Ts[-1] - np.sum(s[:-1] * (Ts[-1] - Ts[:-1])), s
