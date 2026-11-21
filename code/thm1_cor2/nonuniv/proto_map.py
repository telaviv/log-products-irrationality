"""proto_map.py -- FLOAT prototype (mpmath, not a certificate) of the map phi of Section 4.6 of the paper and of its price.

Map (all algebraic, no elliptic functions):
  J(u) = c + d4 (u + 1/u),  c = (y2+y3)/2, d4 = (y3-y2)/4  (u = -1 -> y2, u = 1 -> y3);
  alpha = u0m (inner lift of 0), beta = u4m, gamma = u4p (lifts of y4), k^2 = (beta-alpha)/(gamma-alpha);
  R(zeta) = [gamma(2k zeta - 1 - k)^2 - beta(2 zeta - 1 - k)^2] / [4(1-k^2) zeta (1 - zeta)]
     (= u(v) on the rectangle (-K,3K)x(0,K') written in zeta = s^2, s = sn((v-K)(1+k)/2, 2sqrt(k)/(1+k)); Landen);
  U = C minus ((-inf, zeta1] u [0, inf)),  zeta1 = center-line root of R = u1m (slit at the inner lift of y1);
  Z(z) = zeta1 B^2 (1+z)^2 / (B^2 (1+z)^2 + (1-z)^2): D -> U conformal, Z(0) = zeta0 (center-line root of R = u0p);
  phi = J o R o Z  (rational of degree 8), phi(0) = 0.
Price: L = log|phi'(0)|;  T - L = (2/pi) int_0^{pi/2} N(B tan psi) dpsi, N = sum of Green values at the center-line
preimages (Section 4.6 of the paper).  Cross-check: T by direct double sampling of log|phi(z) - phi(w)| on the circle.
Usage: python proto_map.py q [q ...]
"""
import sys, json
from mpmath import mp, mpf, polyroots, sqrt, log, pi, tan, quad, mpc, exp, fabs, findroot
mp.dps = 40

def ys(q):
    W = q * q
    co = [1, 16384*W**3 - 20480*W**2 + 5248*W - 64, 135168*W**2 + 159744*W + 1536, 294912*W - 16384, 65536]
    rts = polyroots(co, maxsteps=400, extraprec=400)
    return sorted([1 / fabs(r) for r in rts])

def setup(q):
    y1, y2, y3, y4 = ys(q)
    c = (y2 + y3) / 2; d4 = (y3 - y2) / 4
    def lifts(yy):   # (outer, inner)
        s = (yy - c) / d4; r = sqrt(s * s - 4)
        a, b = (s + r) / 2, (s - r) / 2
        return (a, b) if fabs(a) > fabs(b) else (b, a)
    u0p, u0m = lifts(0); u1p, u1m = lifts(y1); u4p, u4m = lifts(y4)
    al, be, ga = u0m, u4m, u4p
    k = sqrt((be - al) / (ga - al))
    def R(z): return (ga * (2*k*z - 1 - k)**2 - be * (2*z - 1 - k)**2) / (4 * (1 - k*k) * z * (1 - z))
    def Rp(z):
        h = mpf(10) ** -15
        return (R(z + h) - R(z - h)) / (2 * h)
    def Rinv_center(ust):
        # quadratic: gamma(2kz-1-k)^2 - beta(2z-1-k)^2 - 4 u (1-k^2) z (1-z) = 0
        A2 = 4*ga*k*k - 4*be + 4*ust*(1 - k*k)
        A1 = -4*ga*k*(1 + k) + 4*be*(1 + k) - 4*ust*(1 - k*k)
        A0 = (ga - be) * (1 + k)**2
        disc = A1*A1 - 4*A2*A0
        r1 = (-A1 + sqrt(disc)) / (2*A2); r2 = (-A1 - sqrt(disc)) / (2*A2)
        neg = [r for r in (r1, r2) if r < 0]
        assert len(neg) == 1, (ust, r1, r2)
        return neg[0]
    z1 = Rinv_center(u1m); z0 = Rinv_center(u0p)
    B = sqrt(z0 / (z1 - z0))
    J = lambda u: c + d4 * (u + 1/u)
    def Z(z): return z1 * B*B * (1 + z)**2 / (B*B * (1 + z)**2 + (1 - z)**2)
    phi = lambda z: J(R(Z(z)))
    Zp0 = 4 * z1 * B*B / (B*B + 1)**2
    L = log(fabs(d4 * (1 - 1/u0p**2) * Rp(z0) * Zp0))
    S = dict(q=q, y=(y1, y2, y3, y4), c=c, d4=d4, u0p=u0p, u0m=u0m, u1p=u1p, u1m=u1m, u4p=u4p, u4m=u4m, k=k, R=R,
             Rinv_center=Rinv_center, z1=z1, z0=z0, B=B, phi=phi, L=L, J=J, Z=Z)
    return S

def N_of_sigma(S, sg):
    """sum of Green values of U at zeta0 over the interior preimages of phi at the boundary point with parameter sigma>0."""
    z1, B, R, u1m, u1p = S['z1'], S['B'], S['R'], S['u1m'], S['u1p']
    eta = sg * sg
    if eta == 1: return mpf('inf')
    zb = z1 * eta / (eta - 1)
    ub = R(zb)
    cands = []
    for ust in (ub, 1 / ub):
        if ust < u1m:   # center-line preimage inside U
            cands.append(ust)
    tot = mpf(0)
    for ust in cands:
        zp = S['Rinv_center'](ust)
        Sg = sqrt(zp / (z1 - zp))
        tot += log(fabs((Sg + B) / (Sg - B)))
    return tot

def T_minus_L(S):
    B = S['B']
    f = lambda psi: N_of_sigma(S, B * tan(psi))
    # breakpoints: sigma = 1 (top/slit junction), top corner, side point with u_b = u0p, bottom/side junction
    z1 = S['z1']; k = S['k']
    def sig_of_zeta(zb): return sqrt(zb / (zb - z1))
    from mpmath import atan
    pts = [mpf(0)]
    kk2 = (1 + k)**2 / (4 * k)      # 1/k2^2: top corner
    for zb in (mpf(1), kk2):
        pts.append(atan(sig_of_zeta(zb) / B))
    pts.append(atan(1 / B))
    pts.append(pi / 2)
    pts = sorted(pts)
    val = quad(f, pts, maxdegree=8)
    return 2 * val / pi

def T_direct(S, N=1500):
    import numpy as np
    phi = S['phi']
    a = [2 * pi * (i + mpf(1) / 2) / N for i in range(N)]
    b = [2 * pi * (i + mpf('0.8183098861837907')) / N for i in range(N)]
    fa = np.array([complex(phi(exp(mpc(0, t)))) for t in a])
    fb = np.array([complex(phi(exp(mpc(0, t)))) for t in b])
    tot = 0.0
    for i in range(N):
        tot += np.sum(np.log(np.abs(fa[i] - fb)))
    return tot / N**2

if __name__ == '__main__':
    tau = 6 * mpf(195) / 196
    for q in map(int, sys.argv[1:]):
        S = setup(q)
        TL = T_minus_L(S)
        M = S['L'] - (S['L'] + TL) / 14 - tau
        print(json.dumps(dict(q=q, L=float(S['L']), T_minus_L=float(TL), M=float(M), k=float(S['k']), B=float(S['B']),
                              zeta0=float(S['z0']), zeta1=float(S['z1']))), flush=True)
        if '--direct' in sys.argv[0:0] or q in (22, 26):
            Td = T_direct(S, 800)
            print('   direct T - L (N=800 sampling, includes the boundary poles) =', Td - float(S['L']), flush=True)
