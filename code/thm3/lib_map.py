"""lib_map.py -- maps psi(z) = z exp(g(z)), g real polynomial (so psi^{-1}(0) = {0} automatically), phi = h o psi;
L = log 256 + g(0); Bost-Charles integral T (shifted trapezoid grids, as in f3_contour.bc_integral) with analytic gradient;
admissibility by winding numbers.  Floats only: [Numerical]."""
import numpy as np, math, cmath

LOG256 = math.log(256.0)
OFF1, OFF2 = 0.21, 0.63

def logh_and_dlog(q, tol=1e-17):
    """log h(q) (principal-ish, only |.| and the derivative matter) and q h'(q)/h(q) for h = -256 q prod(1+q^n)^24."""
    q = np.asarray(q, dtype=complex)
    r = float(np.max(np.abs(q)))
    if r >= 1:
        raise ValueError('|q|>=1')
    nmax = min(int(math.log(tol) / math.log(max(r, 1e-300))) + 2, 20000) if r > 0 else 1
    acc = np.zeros_like(q); dacc = np.zeros_like(q)
    qn = np.ones_like(q)
    for n in range(1, nmax + 1):
        qn = qn * q
        acc += np.log1p(qn)
        dacc += n * qn / (1 + qn)
    logh = np.log(256.0) + 1j * np.pi + np.log(q) + 24 * acc
    return logh, 1 + 24 * dacc        # q h'/h

def h_of(q):
    lh, _ = logh_and_dlog(q)
    return np.exp(lh)

class Map:
    def __init__(self, c):
        self.c = np.asarray(c, dtype=float)
    def g(self, z):
        return np.polynomial.polynomial.polyval(z, self.c)
    def psi(self, z):
        return z * np.exp(self.g(z))
    @property
    def L(self):
        return LOG256 + self.c[0]

def grids(N):
    th1 = 2 * np.pi * (np.arange(N) + OFF1) / N
    th2 = 2 * np.pi * (np.arange(N) + OFF2) / N
    return np.exp(1j * th1), np.exp(1j * th2)

def T_and_grad(c, N, want_grad=True, chunk=None, kappa=2.0):
    z, w = grids(N)
    K = len(c)
    m = Map(c)
    pz, pw = m.psi(z), m.psi(w)
    lz, dz = logh_and_dlog(pz); lw, dw = logh_and_dlog(pw)
    fz, fw = np.exp(lz), np.exp(lw)
    Dz, Dw = fz * dz, fw * dw            # d phi / d g  (since d psi = psi dg, d phi = h'(psi) psi dg = h * (q h'/h) dg)
    # grid-scale regulariser for the gradient only: eps ~ kappa * |d phi/d theta| * (pi/N)
    gp = np.polynomial.polynomial.polyder(m.c) if len(m.c) > 1 else np.zeros(1)
    vz = np.abs(Dz * (1 + z * np.polynomial.polynomial.polyval(z, gp))) * np.pi / N * kappa
    vw = np.abs(Dw * (1 + w * np.polynomial.polynomial.polyval(w, gp))) * np.pi / N * kappa
    tot = 0.0
    Sz = np.zeros(N, dtype=complex); Sw = np.zeros(N, dtype=complex)
    if chunk is None:
        chunk = max(1, 2_000_000 // N)
    for i in range(0, N, chunk):
        df = fz[i:i + chunk, None] - fw[None, :]
        dzw = z[i:i + chunk, None] - w[None, :]
        tot += float(np.sum(np.log(np.abs(df))) - np.sum(np.log(np.abs(dzw))))
        if want_grad:
            eps2 = (vz[i:i + chunk, None] + vw[None, :]) ** 2
            inv = np.conj(df) / (np.abs(df) ** 2 + eps2)
            Sz[i:i + chunk] = inv.sum(axis=1)
            Sw += inv.sum(axis=0)
    T = tot / (N * N)
    if not want_grad:
        return T, None
    # dT/dc_n = (1/N^2) Re[ sum_j Sz_j Dz_j z_j^n - sum_k Sw_k Dw_k w_k^n ]
    a = Sz * Dz; b = Sw * Dw
    V1 = np.vander(z, K, increasing=True); V2 = np.vander(w, K, increasing=True)
    grad = np.real(a @ V1 - b @ V2) / (N * N)
    return T, grad

def winding(c, pts, Nw=16384):
    th = 2 * np.pi * np.arange(Nw + 1) / Nw
    p = Map(c).psi(np.exp(1j * th))
    pts = np.asarray(pts, dtype=complex)
    out = np.zeros(len(pts), dtype=int)
    for i in range(0, len(pts), 64):
        d = p[None, :] - pts[i:i + 64, None]
        ang = np.angle(d[:, 1:] / d[:, :-1])
        out[i:i + 64] = np.rint(ang.sum(axis=1) / (2 * np.pi)).astype(int)
    return out, float(np.max(np.abs(p)))

def min_dist(c, pts, Nw=16384):
    th = 2 * np.pi * np.arange(Nw) / Nw
    p = Map(c).psi(np.exp(1j * th))
    pts = np.asarray(pts, dtype=complex)
    return np.array([np.min(np.abs(p - e)) for e in pts])

def fit_g(psi_fun, K, M=4096):
    """g = log(psi(z)/z) truncated to degree K-1 (real parts of Fourier coefficients; psi assumed real-symmetric)."""
    th = 2 * np.pi * np.arange(M) / M
    z = np.exp(1j * th)
    v = psi_fun(z) / z
    lg = np.log(np.abs(v)) + 1j * np.unwrap(np.angle(v))
    co = np.fft.fft(lg) / M
    return np.real(co[:K])

def rearr(c, N):
    """rearrangement integral iint max(g(s), g(t)), g = log|phi| on T (CDT Thm 6.0.2 numerator, trivial partition)."""
    th = 2 * np.pi * (np.arange(N) + 0.5) / N
    g = np.sort(np.real(logh_and_dlog(Map(c).psi(np.exp(1j * th)))[0]))
    t = (np.arange(N) + 0.5) / N
    return float(np.sum(2 * t * g) / N)
