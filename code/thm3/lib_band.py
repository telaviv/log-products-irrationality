"""lib_band.py -- the BAND version of admissibility: Theorem 14.0.1 for all m != n in Z minus {-1,0} with
delta := |m/n - 1| in [d1, d2].  Lemma B (RESULT.md): the genuine singular values are real and lie in
  F = [-(1/d1 + 8), -(1/d2 - 6)]  U  [1/d2 - 6, 1/d1 + 8]  U  [(4/d2 - 7)^2, +inf),   and |alpha| <= d2^4/16 =: Y0.
phi = h o psi is admissible for every pair of the band iff (i) phi(T) misses F and |phi| > 2 Y0 on T, (ii) the winding number of
phi(T) about one point of each finite segment is 0 (the far ray is joined to infinity), (iii) psi^{-1}(0) = {0}, |psi| < 1.
Float check here [Numerical]; certified version cert_band.py."""
import numpy as np, math
from lib_map import *
from lib_opt import TAU, MFUN


class Band:
    def __init__(self, d1, d2):
        self.d1, self.d2 = d1, d2
        self.lo = 1 / d2 - 6
        self.hi = 1 / d1 + 8
        self.far = (4 / d2 - 7) ** 2
        self.logY0x2 = 4 * math.log(d2) - math.log(16) + math.log(2)
        self.llo, self.lhi, self.lfar = math.log(self.lo), math.log(self.hi), math.log(self.far)
        self.ystar = math.sqrt(self.lo * self.hi)


def band_admissible(c, B, qlim=0.97, Nw=16384, maxN=2 ** 18, safety=0.0):
    while True:
        th = 2 * np.pi * np.arange(Nw) / Nw
        ps = Map(c).psi(np.exp(1j * th))
        if not np.all(np.isfinite(ps)) or np.max(np.abs(ps)) >= qlim:
            return False, 'qlim'
        lh, _ = logh_and_dlog(ps)
        ell = lh.real
        if np.min(ell) <= B.logY0x2:
            return False, 'small'
        if np.max(ell) >= B.lfar - 1:
            return False, 'far'
        A = np.unwrap(np.concatenate([lh.imag, lh.imag[:1]]))
        dA = np.diff(A)
        ell1 = np.concatenate([ell, ell[:1]])
        dl = np.abs(np.diff(ell1))
        if np.max(np.abs(dA)) > 0.3 or np.max(dl) > 0.5:
            need = Nw * max(np.max(np.abs(dA)) / 0.2, np.max(dl) / 0.3)
            if need > maxN:
                return False, 'unresolved'
            Nw = int(2 ** math.ceil(math.log2(need)))
            continue
        if safety > 0:
            inseg = (ell >= B.llo - 0.05) & (ell <= B.lhi + 0.05)
            if np.any(np.abs(np.sin(lh.imag[inseg])) < safety):
                return False, 'graze'
        k0 = np.floor(A[:-1] / np.pi); k1 = np.floor(A[1:] / np.pi)
        for j in np.nonzero(k0 != k1)[0]:
            for k in range(int(min(k0[j], k1[j])) + 1, int(max(k0[j], k1[j])) + 1):
                t = (k * np.pi - A[j]) / (A[j + 1] - A[j])
                lc = ell1[j] + t * (ell1[j + 1] - ell1[j])
                if B.llo - 1e-9 <= lc <= B.lhi + 1e-9:
                    return False, 'segment'
        if abs(A[-1] - A[0] - 2 * np.pi) > 1e-6:
            return False, 'winding0'
        # winding about +ystar and -ystar
        y = np.exp(lh)
        for ys in (B.ystar, -B.ystar):
            d = np.concatenate([y - ys, y[:1] - ys])
            w = np.sum(np.angle(d[1:] / d[:-1])) / (2 * np.pi)
            if abs(w) > 0.5:
                return False, 'wind_seg'
        return True, 'ok'


def band_barrier(c, B, N):
    z = np.exp(2j * np.pi * (np.arange(N) + 0.37) / N)
    K = len(c)
    ps = Map(c).psi(z)
    lh, dlg = logh_and_dlog(ps)
    V = np.vander(z, K, increasing=True)
    W = dlg[:, None] * V                     # d(log phi)/dc_n = dlog * z^n
    ell = lh.real; th = lh.imag
    a = np.maximum(0.0, np.maximum(B.llo - ell, ell - B.lhi))
    sa = np.where(ell < B.llo, -1.0, np.where(ell > B.lhi, 1.0, 0.0))   # d a / d ell
    tot = 0.0; g = np.zeros(K)
    for kind in (0, 1):
        s = 2 * np.sin(th / 2) if kind == 0 else 2 * np.cos(th / 2)
        ds = np.cos(th / 2) if kind == 0 else -np.sin(th / 2)
        d2 = a * a + s * s
        tot += np.mean(np.log(d2) - np.log1p(d2))
        coef = 1.0 / d2 - 1.0 / (1.0 + d2)
        gre = 2 * a * sa * coef
        gim = 2 * s * ds * coef
        g += (gre @ W.real + gim @ W.imag) / N
    return tot, g


class BandProblem:
    def __init__(self, B, N=1024, mu=1e-3, m=MFUN, tau=TAU, qlim=0.97, safety=0.0):
        self.B = B; self.N = N; self.mu = mu; self.m = m; self.tau = tau; self.qlim = qlim; self.safety = safety
        self.nev = 0; self.best = None

    def admissible(self, c):
        return band_admissible(c, self.B, self.qlim, safety=self.safety)

    def f(self, c):
        self.nev += 1
        ok, why = self.admissible(c)
        if not ok:
            return 1e3, np.zeros_like(c)
        T, gT = T_and_grad(c, self.N)
        M = LOG256 + c[0] - T / self.m - self.tau
        gM = -gT / self.m; gM[0] += 1.0
        bv, bg = band_barrier(c, self.B, self.N)
        if self.best is None or M > self.best[0]:
            self.best = (M, c.copy(), T)
        return -(M + self.mu * bv), -(gM + self.mu * bg)
