"""cert_map.py -- CERTIFIED construction of the map phi of Section 4.6 of the paper, its admissibility table, L = log|phi'(0)| and an upper
bound for T(phi) (CDT (7.0.3)), hence a certified lower bound for M = L - T/14 - tau.  All arithmetic in Arb balls
(python-flint arb, ctx.prec bits); every inequality used is an Arb comparison that is certainly true; the only inputs are
the exact integer y-quartic (PROOF_I1 Sec. 7 / unit 1: c34prim(y) = 65536 y^4 - a1 y^3 + a2 y^2 - a3 y + 1).

Objects (proofs: Section 4.6 of the paper):
  y1 < y2 < y3 < y4: the real roots of the y-quartic (certified by exact sign changes over Q).
  J(u) = c + d4 (u + 1/u), c = (y2+y3)/2, d4 = (y3-y2)/4.  Lifts: u_outer * u_inner = 1.
  alpha = u0m, beta = u4m, gamma = u4p, k = sqrt((beta-alpha)/(gamma-alpha)).
  R(zeta) = Nz/Dz,  Nz = gamma(2k zeta-1-k)^2 - beta(2 zeta-1-k)^2,  Dz = 4(1-k^2) zeta(1-zeta);  in w = 1/zeta:
  R = Nw/Dw, Nw = gamma(2k-(1+k)w)^2 - beta(2-(1+k)w)^2, Dw = 4(1-k^2)(w-1).
  zeta1 = center-line root of R = u1m, zeta0 = center-line root of R = u0p, B = sqrt(zeta0/(zeta1-zeta0)).
  U = C minus ((-inf, zeta1] u [0, inf));  Z(z) = zeta1 B^2 (1+z)^2/(B^2(1+z)^2 + (1-z)^2);  phi = J o R o Z.
Boundary parametrisation: z = e^{i t} <-> sigma = i B (1+z)/(1-z) real <-> psi = atan(sigma/B), uniform measure dpsi/pi
on (-pi/2, pi/2); by the symmetry sigma -> -sigma, T - L = (2/pi) int_0^{pi/2} N(psi) dpsi.
  w(psi) = (1 - cot(psi)^2/B^2)/zeta1 (= 1/zeta_b), u_b = R(zeta_b).
  N(psi) = G(u_b) [side only] + G(1/u_b) [side, top, slit], G(u*) = g_U(R_center^{-1}(u*), zeta0) if u* < u1m, else 0,
  g_U(zeta', zeta0) = log|(S+B)/(S-B)|, S = sqrt(zeta'/(zeta1 - zeta')).
Upper sum: on each cell, sup N <= sum of G_sup over the hull of endpoint values (all maps are monotone on the pieces);
windows around the three log singularities are bounded by the mean value theorem (derivative chain enclosed in Arb).
Usage: python cert_map.py q [prec] [tol]      Output: runs/cert_map_q<q>.json
"""
import sys, json, time, os
from flint import arb, fmpz_poly, fmpq, fmpz, ctx

HERE = os.path.dirname(os.path.abspath(__file__))
q = int(sys.argv[1]); PREC = int(sys.argv[2]) if len(sys.argv) > 2 else 512
TOL = float(sys.argv[3]) if len(sys.argv) > 3 else 2e-6
ctx.prec = PREC
T0 = time.time()
def log(*a): print('[%6.1fs]' % (time.time() - T0), *a, flush=True)
def sure(cond, msg):
    if not cond: raise AssertionError('CERTIFICATE FAILED: ' + msg)
def hull(*xs):
    h = xs[0]
    for x in xs[1:]: h = arb.union(h, x)
    return h
def to_fmpq(x):
    m, e = x.mid().man_exp()
    m = fmpz(m); e = int(e)
    return fmpq(m) * fmpq(2) ** e if e >= 0 else fmpq(m, fmpz(2) ** (-e))

# ---------------------------------------------------------------- 1. certified singular points
W = q * q
a3 = 16384 * W**3 - 20480 * W**2 + 5248 * W - 64
a2 = 135168 * W**2 + 159744 * W + 1536
a1 = 294912 * W - 16384
Yp = fmpz_poly([1, -a3, a2, -a1, 65536])
approx = sorted([r for r, mult in Yp.complex_roots()], key=lambda r: float(r.real.mid()))
ys = []
for r in approx:
    m = to_fmpq(r.real)
    eps = abs(m) * fmpq(1, 10**120) + fmpq(1, 10**300)
    lo, hi = m - eps, m + eps
    vlo, vhi = Yp(lo), Yp(hi)
    sure(vlo * vhi < 0, 'sign change for a root of the y-quartic')
    ys.append(arb.union(arb(lo), arb(hi)))
sure(len(ys) == 4, 'four real roots')
y1, y2, y3, y4 = ys
sure(arb(0) < y1 and y1 < y2 and y2 < y3 and y3 < y4, 'ordering 0 < y1 < y2 < y3 < y4')
log('q=%d y = %s' % (q, [y.str(12) for y in ys]))

c = (y2 + y3) / 2; d4 = (y3 - y2) / 4
def J(u): return c + d4 * (u + 1 / u)
def lifts(yy):
    s = (yy - c) / d4; r = (s * s - 4).sqrt()
    if s < arb(-2):
        outer = (s - r) / 2
    else:
        sure(s > arb(2), 'lift parameter outside [-2,2]'); outer = (s + r) / 2
    return outer, 1 / outer
u0p, u0m = lifts(arb(0)); u1p, u1m = lifts(y1); u4p, u4m = lifts(y4)
al, be, ga = u0m, u4m, u4p
sure(u0p < u1p and u1p < arb(-1) and arb(-1) < u1m and u1m < u0m and u0m < 0, 'u0p < u1p < -1 < u1m < u0m < 0')
sure(arb(0) < u4m and u4m < 1 and arb(1) < u4p, '0 < u4m < 1 < u4p')
k = ((be - al) / (ga - al)).sqrt()
sure(arb(0) < k and k < 1, '0 < k < 1')
kk = 1 - k * k
zc = (1 + k) ** 2 / (4 * k)          # 1/k2^2: top corners

def Nz(z): return ga * (2 * k * z - 1 - k) ** 2 - be * (2 * z - 1 - k) ** 2
def Dz(z): return 4 * kk * z * (1 - z)
def R(z): return Nz(z) / Dz(z)
def dNz(z): return ga * 4 * k * (2 * k * z - 1 - k) - be * 4 * (2 * z - 1 - k)
def dDz(z): return 4 * kk * (1 - 2 * z)
def dR(z): return (dNz(z) * Dz(z) - Nz(z) * dDz(z)) / Dz(z) ** 2
def Nw(w): return ga * (2 * k - (1 + k) * w) ** 2 - be * (2 - (1 + k) * w) ** 2
def Dw(w): return 4 * kk * (w - 1)
def dNw(w): return -2 * (1 + k) * (ga * (2 * k - (1 + k) * w) - be * (2 - (1 + k) * w))
def dRw(w): return (dNw(w) * Dw(w) - Nw(w) * 4 * kk) / Dw(w) ** 2

def quad_coeffs(ust):        # Nz - ust*Dz = A2 z^2 + A1 z + A0
    A2 = 4 * ga * k * k - 4 * be + 4 * ust * kk
    A1 = -4 * ga * k * (1 + k) + 4 * be * (1 + k) - 4 * ust * kk
    A0 = (ga - be) * (1 + k) ** 2
    return A2, A1, A0
def Rinv_center(ust):
    """the unique negative root of R(z) = ust, for ust < alpha (A2 < 0 < A0 => one negative, one positive root)."""
    A2, A1, A0 = quad_coeffs(ust)
    sure(A2 < 0, 'A2 < 0 (ust < alpha)')
    D = A1 * A1 - 4 * A2 * A0
    sq = D.sqrt()
    if A1 > 0: z = 2 * A0 / (-A1 - sq)
    elif A1 < 0: z = (-A1 + sq) / (2 * A2)
    else: z = (-A1 + sq) / (2 * A2)
    return z
def roots_R(ust):
    """both roots of R(z) = ust as (complex) pairs: returns ('real', r1, r2) or ('complex', re, im)."""
    A2, A1, A0 = quad_coeffs(ust)
    D = A1 * A1 - 4 * A2 * A0
    if D > 0:
        sq = D.sqrt(); return ('real', (-A1 - sq) / (2 * A2), (-A1 + sq) / (2 * A2))
    if D < 0:
        return ('complex', -A1 / (2 * A2), (-D).sqrt() / (2 * abs(A2)))
    return ('undecided', D)

z1 = Rinv_center(u1m); z0 = Rinv_center(u0p); zm1 = Rinv_center(arb(-1)); z1p = Rinv_center(u1p)
sure(z1 < zm1 and zm1 < z1p and z1p < z0 and z0 < 0, 'zeta1 < zeta(-1) < zeta(1p) < zeta0 < 0')
B = (z0 / (z1 - z0)).sqrt()
log('k=%s zeta0=%s zeta1=%s B=%s' % (k.str(10), z0.str(12), z1.str(8), B.str(8)))

# ---------------------------------------------------------------- 2. critical points of R and the preimage table
# numerator of R': dNz*Dz - Nz*dDz is a quadratic (cubic terms cancel); its roots are the two real critical points
cr = []
# extract coefficients by evaluation at 0, 1, -1 (exact quadratic in z; check cubic coefficient vanishes via a 4th point)
def crit_num(z): return dNz(z) * Dz(z) - Nz(z) * dDz(z)
v0, v1, vm1, v2 = crit_num(arb(0)), crit_num(arb(1)), crit_num(arb(-1)), crit_num(arb(2))
C0 = v0; C2 = (v1 + vm1) / 2 - v0; C1 = (v1 - vm1) / 2
sure(abs(C2 * 4 + C1 * 2 + C0 - v2) < arb(10) ** (-60) * (1 + abs(v2)), 'critical numerator is quadratic')
Dc = C1 * C1 - 4 * C2 * C0
sure(Dc > 0, 'two real critical points')
zcr = sorted([(-C1 - Dc.sqrt()) / (2 * C2), (-C1 + Dc.sqrt()) / (2 * C2)], key=lambda x: float(x.mid()))
z_gam, z_bet = zcr
sure(arb(0) < z_gam and z_gam < 1, 'critical point in (0,1) (bottom)')
sure(zc < z_bet, 'critical point beyond the top corner')
# (referee K-2) exact identities, checked symbolically: R - gamma = (gamma-beta)(2z-1-k)^2/Dz, R - beta = (gamma-beta)(2kz-1-k)^2/Dz,
# so the critical points are (1+k)/2 and (1+k)/(2k) and gamma, beta are taken only there; R(inf) = R(zc) = alpha by the choice of k.
sure(abs(R(z_gam) - ga) < arb(10) ** (-40), 'critical value gamma'); sure(abs(R(z_bet) - be) < arb(10) ** (-40), 'critical value beta')
def inU(z):        # real z certainly in U  <=> z1 < z < 0
    return (z1 < z) and (z < 0)
def outU(z):       # real z certainly outside U
    return (z < z1) or (z > 0) or (z == z1)
table = {}
for name, val in [('u0p (y=0)', u0p), ('u1p (y=y1)', u1p), ('u1m (y=y1)', u1m), ('-1 (y=y2)', arb(-1)), ('+1 (y=y3)', arb(1)),
                  ('u4p (y=y4)', u4p), ('u4m (y=y4)', u4m), ('0 (y=inf)', arb(0))]:
    rr = roots_R(val)
    if rr[0] == 'real':
        r1, r2 = rr[1], rr[2]
        st = []
        for r in (r1, r2):
            if name.startswith('u1m') and r.overlaps(z1) and abs(R(r) - u1m) < arb(10) ** -40:
                st.append('outside U'); continue    # the slit endpoint zeta1 itself (defined as this root)
            st.append('in U' if inU(r) else ('outside U' if outU(r) else 'UNDECIDED'))
        table[name] = dict(kind='real', roots=[r1.str(25), r2.str(25)], status=st)
    elif rr[0] == 'complex':
        table[name] = dict(kind='complex pair (in U)', re=rr[1].str(25), im=rr[2].str(25))
    else:
        table[name] = dict(kind='UNDECIDED')
# u0m = alpha: R(inf) = alpha, R(zc) = alpha (check), u = inf: zeta in {0, 1}
sure(abs(R(zc) - al) < arb(10) ** (-40), 'R(top corner) = alpha')
table['u0m=alpha (y=0)'] = dict(kind='real', roots=['inf', zc.str(25)], status=['outside U', 'outside U'])
table['inf (y=inf)'] = dict(kind='real', roots=['0', '1'], status=['outside U', 'outside U'])
# expected pattern
exp_in = {'u0p (y=0)': 1, 'u1p (y=y1)': 1, 'u1m (y=y1)': 0, '-1 (y=y2)': 1, 'u4p (y=y4)': 0, 'u4m (y=y4)': 0, '0 (y=inf)': 0}
for name, nin in exp_in.items():
    t = table[name]
    sure(t['kind'] == 'real' and 'UNDECIDED' not in t['status'] or name in ('u4p (y=y4)', 'u4m (y=y4)'), 'table decided: ' + name)
# u4p, u4m are the critical values gamma, beta: double roots z_gam, z_bet (outside U); +1 must be a complex pair
table['u4p (y=y4)'] = dict(kind='double root', roots=[z_gam.str(25)], status=['outside U'])
table['u4m (y=y4)'] = dict(kind='double root', roots=[z_bet.str(25)], status=['outside U'])
sure(table['+1 (y=y3)']['kind'].startswith('complex'), '+1 has a complex-conjugate pair of preimages')
for name in ('u0p (y=0)', 'u1p (y=y1)', '-1 (y=y2)'):
    sure(table[name]['status'].count('in U') == 1 and table[name]['status'].count('outside U') == 1, name)
for name in ('u1m (y=y1)', '0 (y=inf)'):
    sure(table[name]['status'].count('outside U') == 2, name)
# the in-U roots are the center-line points found above
sure(abs(R(z0) - u0p) < arb(10) ** -40 and abs(R(z1p) - u1p) < arb(10) ** -40 and abs(R(zm1) + 1) < arb(10) ** -40, 'center-line roots')
sure(abs(dR(zm1)) > 0, "R'(zeta(-1)) != 0")
log('preimage table certified')

# ---------------------------------------------------------------- 3. L
Zp0 = 4 * z1 * B * B / (B * B + 1) ** 2
phip0 = d4 * (1 - 1 / (u0p * u0p)) * dR(z0) * Zp0
Lval = abs(phip0).log()
log('L = %s' % Lval.str(20))

# ---------------------------------------------------------------- 4. T - L
def w_of_psi(psi):
    ct = psi.cot() if not psi.contains(arb.pi() / 2) else None
    if ct is None: raise ValueError
    return (1 - ct * ct / (B * B)) / z1
def ub_of_w(w): return Nw(w) / Dw(w)
def g_of_zeta(zp):
    S = (zp / (z1 - zp)).sqrt()
    return ((S + B) / abs(S - B)).log()
INF = None
def G_point(ust):   # ust certainly < u1m and != u0p
    return g_of_zeta(Rinv_center(ust))
def G_sup(lo, hi):
    """upper bound of G on the u*-interval [lo, hi]; lo = None means -infinity.  G = 0 on [u1m, inf), increasing on
    (-inf, u0p), decreasing on (u0p, u1m).  Returns arb upper bound or arb('inf')."""
    if hi < u0p:
        return G_point(hi).upper()
    if lo is not None and lo > u0p:
        if lo > u1m or lo == u1m: return arb(0)
        if lo < u1m: return G_point(lo).upper()
        lo2 = arb(lo.lower())
        if lo2 > u0p and lo2 < u1m: return G_point(lo2).upper()
        return arb('inf')
    return arb('inf')
def cand_ranges(ubL, ubR, side):
    """u*-ranges of the candidates over a cell with endpoint values ubL, ubR (u_b monotone on the cell);
    ubL/ubR = 'NEGINF' allowed at the side start.  Returns list of (lo, hi) with lo possibly None (= -inf)."""
    out = []
    vals = [v for v in (ubL, ubR) if v != 'NEGINF']
    neginf = len(vals) < 2
    lo = None if neginf else (ubL if ubL < ubR else ubR)
    hi = vals[0] if neginf else (ubR if ubL < ubR else ubL)
    if not neginf and not (ubL < ubR or ubR < ubL):
        lo = hull(ubL, ubR); hi = lo; lo = arb(lo.lower()); hi = arb(hi.upper())
    if side:                          # candidate A: u* = u_b
        out.append((lo, hi))
    # candidate B: u* = 1/u_b
    if hi < 0:
        out.append((1 / hi, (1 / lo) if lo is not None else arb(0)))   # [1/hi, 1/lo]
    elif lo is not None and lo < 0 and hi > 0:
        out.append((None, 1 / lo))
    elif lo is not None and lo > 0:
        pass
    else:
        raise ValueError('cand B range undecided')
    return out
def N_sup(ubL, ubR, side):
    s = arb(0)
    for lo, hi in cand_ranges(ubL, ubR, side):
        if lo is not None and (lo > u1m or lo == u1m):
            continue
        s += G_sup(lo, hi)
    return s
def N_est(ub, side):   # midpoint estimate (float-level, for refinement decisions only)
    s = 0.0
    for ust in ([ub] if side else []) + [1 / ub]:
        if ust < u1m and not ust.overlaps(u0p):
            try: s += float(G_point(ust).mid())
            except Exception: pass
    return s

PI = arb.pi()
def psi_of_zeta(zb):     # boundary point with zeta_b real (sigma>0 half)
    sg = (zb / (zb - z1)).sqrt()
    return (sg / B).atan()
psi1 = psi_of_zeta(arb(1))                       # bottom/side junction
zs = roots_R(u0p); sure(zs[0] == 'real', 'side root of R=u0p')
zside = zs[1] if (zs[1] > 1) else zs[2]
sure(arb(1) < zside and zside < zc, 'side root of R = u0p lies on the side')
psis = psi_of_zeta(zside)                        # side point with u_b = u0p
psic = psi_of_zeta(zc)                           # top corner
r0 = roots_R(arb(0)); sure(r0[0] == 'real', 'R=0 real roots')
z0a, z0b = sorted([r0[1], r0[2]], key=lambda x: float(x.mid()))
sure(zc < z0a and z0a < z_bet and z_bet < z0b, 'zeros of R on the top, either side of the turning point')
psi0a = psi_of_zeta(z0a); psi0b = psi_of_zeta(z0b)
psiE = (1 / B).atan()                            # top/slit junction (sigma = 1)
sure(psi1 < psis and psis < psic and psic < psi0a and psi0a < psi0b and psi0b < psiE and psiE < PI / 2, 'breakpoint order')
log('breakpoints psi: 1:%s s:%s c:%s 0a:%s 0b:%s E:pi/2-%s' % (psi1.str(8), psis.str(8), psic.str(8), psi0a.str(8),
    psi0b.str(8), (PI / 2 - psiE).str(6)))

def state(psi):
    w = w_of_psi(psi)
    return ub_of_w(w)

total = arb(0); ncell = [0]; EXC = [0.0]
def integrate_piece(a, b, side, label, fa=None, fb=None):
    """upper sum over [a, b] (a, b exact arb points), adaptive bisection; fa/fb: u_b at the ends (or 'NEGINF')."""
    if fa is None: fa = state(a)
    if fb is None: fb = state(b)
    stack = [(a, b, fa, fb, 0)]
    acc = arb(0)
    while stack:
        x, y, ux, uy, dep = stack.pop()
        width = y - x
        sup = N_sup(ux, uy, side)
        mid = (x + y) / 2
        um = state(mid)
        est = N_est(um, side)
        excess = (float(sup.upper()) - est) * float(width.mid()) if sup.is_finite() else float('inf')
        if excess > TOL and dep < 80:
            stack.append((x, mid, ux, um, dep + 1)); stack.append((mid, y, um, uy, dep + 1))
            continue
        sure(sup.is_finite(), 'finite cell bound in piece ' + label)
        acc += width * sup; ncell[0] += 1; EXC[0] += max(excess, 0.0)
    log('piece %-8s [%s, %s]: upper %s (cells so far %d)' % (label, a.str(6), b.str(6), acc.upper().str(8), ncell[0]))
    return acc.upper()

def deriv_lower(psiL, psiR, cand, side):
    """lower bound of |dS/dpsi| on [psiL, psiR] for the singular candidate (all intermediate variables monotone)."""
    def chain_vals(p):
        w = w_of_psi(p); ub = ub_of_w(w); ust = ub if cand == 'A' else 1 / ub
        zp = Rinv_center(ust); S = (zp / (z1 - zp)).sqrt()
        return dict(p=p, w=w, ub=ub, ust=ust, zp=zp, S=S)
    vL, vR = chain_vals(psiL), chain_vals(psiR)
    H = {key: hull(vL[key], vR[key]) for key in ('p', 'w', 'ub', 'ust', 'zp', 'S')}
    p = H['p']
    # dw/dpsi: w = (1 - cot^2/B^2)/z1 -> dw/dpsi = (2 cot csc^2 / B^2)/z1
    ct = p.cot(); cs2 = 1 + ct * ct
    dw = 2 * ct * cs2 / (B * B) / z1
    dub = dRw(H['w'])
    dust = arb(1) if cand == 'A' else -1 / (H['ub'] * H['ub'])
    dzp = 1 / dR(H['zp'])
    dS = z1 / (2 * H['S'] * (z1 - H['zp']) ** 2)
    prod = dw * dub * dust * dzp * dS
    sure(not prod.contains(0), 'derivative chain bounded away from 0 in a window')
    return abs(prod).lower(), H['S'].upper()

def window(center, delta, cand, side, label):
    """bound of int N over [center-delta, center+delta] (center = singular point ball, radius << delta)."""
    L_, R_ = arb(center.mid()) - delta, arb(center.mid()) + delta
    sure(center.lower() > L_ and center.upper() < R_, 'singular point inside window ' + label)
    ell, Smax = deriv_lower(L_, R_, cand, side)
    width = R_ - L_
    # singular candidate: g <= log(Smax + B) + log(1/(ell |psi - psi_p|)); int over window of log(1/|x|) <= width(1 + log(2/width))
    bsing = width * ((Smax + B).log() + (1 / arb(ell)).log() + 1 + (2 / width).log())
    # the other candidate on the window: sup by G_sup over endpoint hull
    uL, uR = state(L_), state(R_)
    other = arb(0)
    for (lo, hi), cnd in zip(cand_ranges(uL, uR, side), (['A', 'B'] if side else ['B'])):
        if cnd == cand: continue
        if lo is not None and (lo > u1m or lo == u1m): continue
        other += G_sup(lo, hi)
    sure(other.is_finite(), 'window other-candidate bound finite ' + label)
    tot = (bsing + width * other).upper()
    log('window %-4s half-width %s: bound %s (ell=%s)' % (label, delta.str(3), tot.str(8), arb(ell).str(5)))
    return tot, L_, R_

parts = {}
dS = arb('1e-7'); dE = B * arb('1e-4')
ws, sL, sR = window(psis, dS, 'A', True, 's')
wc, cL, cR = window(psic, dS, 'B', True, 'c')       # side and top-a share the monotone stretch through the corner
wE, eL, eR = window(psiE, dE, 'B', True, 'E')
pa = arb(psi1.lower())      # the first side cell also covers the bottom/side junction (bottom part: N = 0)
sure(Nz(arb(0)) > 0 and Nz(arb(1)) > 0, 'R -> +inf at both ends of (0,1): with the single critical value gamma, R >= gamma on the bottom (N = 0 there)')
parts['side1'] = integrate_piece(pa, sL, True, 'side1', fa='NEGINF')
parts['win_s'] = ws
parts['side2'] = integrate_piece(sR, cL, True, 'side2')
parts['win_c'] = wc
# top-a: from cR to psi0a (u_b increases to 0-); cell rule uses only the left value near psi0a: end at psi0a.lower,
# then the tiny cell [psi0a.lower, psi0a.upper] has candidate-B range (-inf, 1/u_b(psi0a.lower)]
EPS0 = arb(10) ** -60
a0 = arb(psi0a.mid()) - EPS0; b0 = arb(psi0a.mid()) + EPS0
sure(a0 < psi0a and psi0a < b0, 'tiny cell around the first zero')
parts['topa'] = integrate_piece(cR, a0, True, 'topa')
u_a0 = state(a0); sure(u_a0 < 0, 'u_b < 0 just before the zero')
parts['topa_tiny'] = ((b0 - a0) * G_sup(None, 1 / u_a0)).upper()
# [psi0a.upper, psi0b.lower]: u_b in (0, beta]: N = 0 (no candidate)
a1_ = arb(psi0b.mid()) + EPS0; b1_ = arb(psi0b.mid()) - EPS0
sure(b1_ < psi0b and psi0b < a1_, 'tiny cell around the second zero')
# (referee K-4) the pieces and windows tile [pa, pi/2] in order, with no overlap (an overlap would subtract area)
sure(pa < sL and sR < cL and cR < a0 and b0 < b1_ and a1_ < eL and eR < PI / 2, 'pieces and windows ordered, no overlap')
u_b1 = state(a1_); sure(u_b1 < 0, 'u_b < 0 just after the second zero')
sure(state(b0) > 0 and state(b1_) > 0, 'u_b > 0 between the zeros (no candidate there)')
parts['topb_tiny'] = ((a1_ - b1_) * G_sup(None, 1 / u_b1)).upper()
parts['topb'] = integrate_piece(a1_, eL, True, 'topb')
parts['win_E'] = wE
# slit: [eR, pi/2]; at pi/2, w = 1/zeta1 exactly
def state_end(): return ub_of_w(1 / z1)
parts['slit'] = integrate_piece(eR, PI / 2, True, 'slit', fb=state_end())
S_int = arb(0)
for key, v in parts.items(): S_int += arb(v)
TL = (2 * S_int / PI).upper()
tau = arb(6) * 195 / 196
T_up = Lval.upper() + TL
M = Lval - T_up / 14 - tau
res = dict(q=q, prec=PREC, tol=TOL, y=[y.str(30) for y in ys], k=k.str(20), zeta0=z0.str(20), zeta1=z1.str(20), B=B.str(20),
           L=Lval.str(30), T_minus_L_upper=TL.str(15), T_upper=T_up.str(15), tau=tau.str(20), M_lower=arb(M.lower()).str(15),
           M=M.str(15), M_positive=bool(M > 0), cells=ncell[0], excess_estimate=EXC[0], parts={kk_: arb(v).str(8) for kk_, v in parts.items()},
           preimage_table=table, crit=[z_gam.str(12), z_bet.str(12)], seconds=round(time.time() - T0, 1))
log('estimated excess of the upper sum (diagnostic): %.3e' % (2 * EXC[0] / 3.14159))
log('T - L <= %s ; L = %s ; M >= %s ; positive: %s' % (TL.str(10), Lval.str(15), arb(M.lower()).str(10), M > 0))
json.dump(res, open(os.path.join(HERE, 'runs', 'cert_map_q%d.json' % q), 'w'), indent=1)
