"""lib_range.py -- the RANGE version of the admissibility condition for Theorem 14.0.1 at |m/n - 1| < eps:
every singular value other than alpha lies on the real rays |y| >= Y1 = 1/eps, alpha in [-Y0, Y0], Y0 = eps^4/16 (CDT Sec. 14.4).
phi = h o psi is admissible for ALL such pairs at once iff (i) phi(T) never meets the rays |y| >= Y1 (then every point of the rays
has winding number 0, i.e. no preimage), (ii) |phi| > 2 Y0 on T (then each point of D(0, 2Y0) has exactly one preimage: the
univalent leaf of Cor. 9.0.19 with U = D(0, 2 Y0)), (iii) psi^{-1}(0) = {0} (built in) and |psi| < 1.
Checked here on a fine sample of T by unwrapping arg phi (floats: [Numerical]; the certified version is cert_*.py)."""
import numpy as np, math, time
from lib_map import *
from lib_opt import TAU, MFUN

def y_on_circle(c, Nw, off=0.0):
    th = 2 * np.pi * (np.arange(Nw) + off) / Nw
    z = np.exp(1j * th)
    ps = Map(c).psi(z)
    lh, dl = logh_and_dlog(ps)
    return z, ps, lh, dl

def range_admissible(c, logY1, logY0x2, qlim=0.97, Nw=16384, maxN=2 ** 18, safety=0.0):
    while True:
        th = 2 * np.pi * np.arange(Nw) / Nw
        ps = Map(c).psi(np.exp(1j * th))
        if not np.all(np.isfinite(ps)) or np.max(np.abs(ps)) >= qlim:
            return False, 'qlim'
        lh, _ = logh_and_dlog(ps)
        ell = lh.real
        if np.min(ell) <= logY0x2:
            return False, 'small'
        A = np.unwrap(np.concatenate([lh.imag, lh.imag[:1]]))
        dA = np.diff(A)
        if np.max(np.abs(dA)) > 0.3 or np.max(np.abs(np.diff(np.concatenate([ell, ell[:1]])))) > 0.5:
            need = Nw * max(np.max(np.abs(dA)) / 0.2, np.max(np.abs(np.diff(np.concatenate([ell, ell[:1]])))) / 0.3)
            if need > maxN:
                return False, 'unresolved'
            Nw = int(2 ** math.ceil(math.log2(need)))
            continue
        if safety > 0:
            near = ell >= logY1 - 0.2
            if np.any(np.abs(np.sin(lh.imag[near])) < safety):
                return False, 'graze'
        k0 = np.floor(A[:-1] / np.pi); k1 = np.floor(A[1:] / np.pi)
        idx = np.nonzero(k0 != k1)[0]
        ell1 = np.concatenate([ell, ell[:1]])
        for j in idx:
            for k in range(int(min(k0[j], k1[j])) + 1, int(max(k0[j], k1[j])) + 1):
                t = (k * np.pi - A[j]) / (A[j + 1] - A[j])
                lc = ell1[j] + t * (ell1[j + 1] - ell1[j])
                if lc >= logY1 - 1e-9:
                    return False, 'ray'
        if abs(A[-1] - A[0] - 2 * np.pi) > 1e-6:
            return False, 'winding'
        return True, 'ok'

def range_barrier(c, logY1, N):
    z, ps, lh, dl = y_on_circle(c, N, 0.37)
    K = len(c)
    V = np.vander(z, K, increasing=True)
    ell = lh.real; th = lh.imag
    s = np.sin(th); co = np.cos(th)
    a = np.maximum(0.0, logY1 - ell)
    d2 = a * a + s * s
    B = np.mean(np.log(d2) - np.log1p(d2))
    # d(d2) = -2 a dell + 2 s cos(th) dth ; dell = Re(v dg), dth = Im(v dg), v = dl
    coef = (1.0 / d2 - 1.0 / (1.0 + d2))
    gre = -2 * a * coef      # multiplies Re(v z^n)
    gim = 2 * s * co * coef  # multiplies Im(v z^n)
    W = dl[:, None] * V
    g = (gre @ W.real + gim @ W.imag) / N
    return B, g

class RangeProblem:
    def __init__(self, eps, N=1024, mu=1e-3, m=MFUN, tau=TAU, qlim=0.97, safety=0.0):
        self.safety = safety
        self.logY1 = -math.log(eps); self.logY0x2 = 4 * math.log(eps) - math.log(16) + math.log(2)
        self.N = N; self.mu = mu; self.m = m; self.tau = tau; self.qlim = qlim; self.nev = 0; self.best = None
    def admissible(self, c):
        ok, why = range_admissible(c, self.logY1, self.logY0x2, self.qlim, safety=self.safety)
        return ok, why
    def f(self, c):
        self.nev += 1
        ok, why = self.admissible(c)
        if not ok:
            return 1e3, np.zeros_like(c)
        T, gT = T_and_grad(c, self.N)
        M = LOG256 + c[0] - T / self.m - self.tau
        gM = -gT / self.m; gM[0] += 1.0
        bv, bg = range_barrier(c, self.logY1, self.N)
        if self.best is None or M > self.best[0]:
            self.best = (M, c.copy(), T)
        return -(M + self.mu * bv), -(gM + self.mu * bg)
