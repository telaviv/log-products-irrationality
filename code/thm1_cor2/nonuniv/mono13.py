"""mono13.py -- CERTIFIED continuation of g = y^MU * H_q with the minimal operator L13 (exact recurrence of order RHO = 25,
n-degree R = 13: runs/op_q<q>_25_13.json) and a certified lower bound for dim span{continuations of H_q}.

Operator for g: recurrence polys p_j(n - MU); L = sum_j y^(RHO-j) p_j(theta - j - MU) = sum_i A_i(y) D^i, A_i = y^i B_i(y).
L g = E(y), E polynomial computed from g's Taylor coefficients (checked: higher coefficients vanish).
Singular points: 0, the 25 roots of the leading coefficient (y1..y4 and 21 apparent points), infinity.

Rigour (all in Arb balls; no midpoint truncation):
* initial data at b = 2^-21: g's series with |R3_n(q)| <= (2n+1) 2^(20n) q^(-2n) (4 log 2)/(q-1)  [RF1-F1: R3_n(q) =
  int Q_n dmu3/(q-x) = q^(-2n) int Q_n(x) x^(2n) dmu3/(q-x), using the orthogonality int Q_n x^k dmu3 = 0 for k < 2n (Lysov's
  type II conditions, multi-index (2n,2n,2n); PROVED as PROOF_I1 Lemma 4, eqs. (O12), (O3)); then |Q_n| <= (2n+1)2^(20n) on
  [-1,1], |x^(2n)/(q-x)| <= 1/(q-1), |mu3| = 4 log 2.  The earlier text "no orthogonality used" was wrong: without q^(-2n)
  the kap < 0.5 assertion fails at b = 2^-21.  Orthogonality-free variant: b = 2^-22 and no q^(-2n) (referee rerun,
  runs/referee_rerun/, env NOORTH=1 there), which certifies the same d = 14 at q = 26 and 32.]
* Taylor step c -> c+h: coefficients by the recurrence in balls; tail by the majorant |y_j| <= A s^-j (all j), valid when
  F(n0) <= 1, F(N) = M_e s^R/(N+1)_R + (1-s/sigma)^-1 sum_{i<R} M_i s^(R-i) (N+1)_i/(N+1)_R, with
  M_i >= max_{|y-c|=sigma} |A_i/A_R|, M_e >= max |E/A_R| enclosed by arc subdivision (each arc inside a ball; acb evaluation).
* rank: det of the 14 x 14 matrix of continuation vectors (13 Taylor coefficients + affine coordinate) excludes 0.
Usage: python mono13.py q [prec] [MU] [tag] [trunc_bits] [diag]   (certificate runs: q 2200 5 _cert 1700 diag, runcert.sh)
Version history: mono13_v1.py = the script that produced the filed q = 26..32 certificates (q = 32 rerun 2026-11-19, RF1-F10);
this file = v1 + RF1-F1 docstring + RF1-F3 arb comparisons in circle_max + RF1-F4 arb-string radii / diagnostic label
(determinant-neutral: checked by a q = 26 rerun, runs/mono13_q26_rf3*).
"""
import sys, json, os, time, math, pickle
from fractions import Fraction
from flint import fmpz_poly, fmpq, acb, arb, acb_mat, acb_poly, ctx
HERE = os.path.dirname(os.path.abspath(__file__))
def arbmax(xs):   # RF1-F4: max of nonnegative arb upper bounds, as an arb (prints radii < 1e-308 correctly)
    m = arb(0)
    for x in xs:
        u = arb(x.upper())
        if u > m: m = u
    return m
q = int(sys.argv[1]); PREC = int(sys.argv[2]) if len(sys.argv) > 2 else 600
MU = int(sys.argv[3]) if len(sys.argv) > 3 else 5
TAG = sys.argv[4] if len(sys.argv) > 4 else ''
BEXP = 21
FR, SIG, SS = 0.2, 0.6, 0.45      # |h| <= FR*delta, sigma = SIG*delta, s = SS*delta
NARC = 512
TB = int(sys.argv[5]) if len(sys.argv) > 5 else 260
DIAG = len(sys.argv) > 6 and 'diag' in sys.argv[6]
EXTRA = len(sys.argv) > 6 and 'extra' in sys.argv[6]   # target truncation bits per step
ctx.prec = PREC
T0 = time.time()
def log(*a): print('[%7.1fs]' % (time.time() - T0), *a, flush=True)

dop = json.load(open(os.path.join(HERE, 'runs', 'op_q%d_25_13.json' % q)))
RHO, R = dop['rho'], dop['dg']
assert dop['failures'] == 0
P = [[int(x) for x in pp] for pp in dop['polys']]
PM = [fmpz_poly(p)(fmpz_poly([-MU, 1])) for p in P]          # p_j(n - MU)
def stirling2(K):
    S = [[0] * (K + 1) for _ in range(K + 1)]; S[0][0] = 1
    for k in range(1, K + 1):
        for i in range(1, k + 1): S[k][i] = i * S[k - 1][i] + S[k - 1][i - 1]
    return S
S2 = stirling2(R)
Bc = [[0] * (RHO + 1) for _ in range(R + 1)]
for j in range(RHO + 1):
    sh = PM[j](fmpz_poly([-j, 1])); a = [int(sh[k]) for k in range(R + 1)]
    for i in range(R + 1): Bc[i][RHO - j] += sum(a[k] * S2[k][i] for k in range(i, R + 1))
A_int = [[0] * i + Bc[i] for i in range(R + 1)]
Bpoly = [acb_poly([acb(x) for x in Bc[i]]) for i in range(R + 1)]
lead = fmpz_poly(Bc[R])
roots = [r for r, m in lead.complex_roots()]
assert len(roots) == RHO
rootsf = [complex(float(r.real.mid()), float(r.imag.mid())) for r in roots]
log('leading coefficient roots (abs <= 5000):', sorted([('%.4g%+.4gi' % (z.real, z.imag)) for z in rootsf if abs(z) < 5000], key=len))

# exact Q, P3 to N0 by the order-4 recurrence; R3_n; g coefficients
d4 = json.load(open(os.path.join(HERE, 'runs', 'rec_q%d.json' % q))); P4 = [[int(x) for x in pp] for pp in d4['polys']]
sys.argv = ['rec.py', str(q)]
from rec import exact_values
rate = arb(2) ** 20 * arb(2) ** (-BEXP) / arb(q) ** 2
N0 = int((PREC + 60 + R * 12) / (-math.log2(float(rate.mid())))) + 40
vals = exact_values(3)
Qs = [Fraction(v[0]) for v in vals]; P3s = [Fraction(v[3]) for v in vals]
for n in range(0, N0 - 3):
    pj = [sum(P4[j][e] * n ** e for e in range(35)) for j in range(5)]
    Qs.append(-sum(pj[j] * Qs[n + j] for j in range(4)) / pj[4]); P3s.append(-sum(pj[j] * P3s[n + j] for j in range(4)) / pj[4])
ctx.prec = PREC + 60 * N0 + 200
theta = (arb(q + 1) / q).log() * (arb(q - 1) / q).log()
Rn = [theta * arb(fmpq(x.numerator, x.denominator)) - arb(fmpq(y.numerator, y.denominator)) for x, y in zip(Qs, P3s)]
ctx.prec = PREC
Rn = [arb(x) for x in Rn]
cg = [arb(0)] * MU + Rn                       # g = y^MU H
# E = L g (polynomial): compute coefficients of degree < RHO + MU + 2 and check the next ones vanish
def pval(poly, x):
    return sum(int(poly[k]) * x ** k for k in range(poly.degree() + 1))
Ecoef = [arb(0)] * (RHO + MU + 40)
for j in range(RHO + 1):
    for m in range(min(len(cg), RHO + MU + 40)):
        dd = m + RHO - j
        if 0 <= dd < len(Ecoef):
            Ecoef[dd] += pval(PM[j], m - j) * cg[m]
degE = RHO + MU - 1
resid = max(float(abs(Ecoef[k]).upper()) for k in range(degE + 1, len(Ecoef)))
emax = max(float(abs(Ecoef[k]).upper()) for k in range(degE + 1))
assert all(Ecoef[k].contains(0) for k in range(degE + 1, len(Ecoef))) and resid < emax * 2.0 ** (-PREC // 3), (resid, emax)
Ec = Ecoef[:degE + 1]
log('E: degree <= %d, max |higher coeff| (float; 0.00e+00 = below 1e-308) = %.2e (balls contain 0), max |E coeff| = %.2e' % (degE, resid, emax))
Epoly = acb_poly([acb(e) for e in Ec])

# initial data at b = 2^-BEXP (Taylor coefficients of g, rigorous tail)
b = arb(2) ** (-BEXP)
init = []
for k in range(R):
    s = arb(0)
    for n in range(k, len(cg)):
        s += math.comb(n, k) * cg[n] * b ** (n - k)
    Ntail = len(cg)       # tail over g index n >= Ntail  <->  H index n - MU
    nh = Ntail - MU
    tn = math.comb(Ntail, k) * (2 * nh + 1) * (arb(2) ** 20 / arb(q) ** 2) ** nh * 4 * arb(2).log() / (q - 1) * b ** (Ntail - k)
    kap = rate * arb(Ntail + 1) / (Ntail + 1 - k) * arb(2 * nh + 3) / (2 * nh + 1)
    assert float(kap.upper()) < 0.5
    s += arb(0, (tn / (1 - kap)).upper())
    init.append(acb(s))
log('N0 = %d; initial data at 2^-%d: g_0 = %s, rad max %s' % (N0, BEXP, init[0].str(8), arbmax([x.rad() for x in init]).str(3)))   # RF1-F4

FF = [fmpz_poly([1])]
for i in range(1, R + 1): FF.append(FF[-1] * fmpz_poly([-(i - 1), 1]))
FFa = [acb_poly([acb(int(f[t])) for t in range(f.degree() + 1)]) for f in FF]
def taylor_shift(coeffs, c):
    s_ = acb_poly([acb(x) for x in coeffs])(acb_poly([c, acb(1)]))
    return [s_[m] for m in range(s_.degree() + 1)]
sing_acb = [acb(0)] + roots
def delta_of(c):
    return min([abs(c - z) for z in sing_acb], key=lambda x: float(x.mid()))

# factored numerators: |B_i(y)| = |lc_i| prod |y - r| with certified root balls (avoids interval-Horner blow-up)
Broots, Blc = [], []
for i in range(R):
    fp = fmpz_poly(Bc[i])
    Blc.append(abs(arb(int(fp[fp.degree()]))))
    rr_ = []
    for r, m in fp.complex_roots(): rr_ += [r] * m
    assert len(rr_) == fp.degree()
    Broots.append(rr_)
def circle_bounds(c, sigma):
    M, Me, mn = circle_max(c, sigma, want_min=True)
    return M, Me, mn
def circle_max(c, sigma, want_min=False):
    """rigorous upper bounds M_i (i<R) for |y^(i-R) B_i(y)/B_R(y)| and M_e for |E(y)/(y^R B_R(y))| on |y - c| = sigma.
    Numerators: polynomials Taylor-shifted to c at working precision, evaluated on arc balls t; denominator: factored form
    |B_R(y)| = |lc| prod |y - r_k| with the root balls (rigorous lower bound)."""
    she = acb_poly(taylor_shift(Ec, c))
    lp = ctx.prec; ctx.prec = 256
    M = [arb(0)] * R; Me = arb(0)
    pi2 = arb.pi() * 2
    rad = (sigma * pi2 / NARC * arb('1.01') / 2).upper()
    lcR = abs(arb(Bc[R][RHO]))
    mnAR = None
    for k in range(NARC):
        ang = pi2 * (2 * k + 1) / (2 * NARC)
        tm = acb(sigma * ang.cos(), sigma * ang.sin())
        t = acb(arb(tm.real.mid(), (rad + tm.real.rad()).upper()), arb(tm.imag.mid(), (rad + tm.imag.rad()).upper()))
        y = c + t
        bR = lcR
        for r in roots: bR = bR * abs(y - r).lower()
        bR = bR.lower()
        ya = abs(y).lower()
        assert bR > 0 and ya > 0, 'arc touches a singular point'     # RF1-F3: arb comparisons (no float under/overflow)
        ar = (bR * ya ** R).lower()
        if mnAR is None or ar < mnAR: mnAR = arb(ar)                 # RF1-F3 (exact points: comparison decided)
        for i in range(R):
            num = Blc[i]
            for r in Broots[i]: num = num * abs(y - r).upper()
            v = (num / bR * ya ** (i - R)).upper()
            if v > M[i]: M[i] = arb(v)                               # RF1-F3
        v = (abs(she(t)) / bR * ya ** (-R)).upper()
        if v > Me: Me = arb(v)                                       # RF1-F3
    ctx.prec = lp
    if want_min: return M, Me, mnAR
    return M, Me

NSTEP = [0]; MAXN0 = [0]; SIGE, SE = 0.35, 0.25
def step(X, c, c2, Nt=None):
    """Rigorous Taylor step (residual-majorant form).  yhat: midpoint recurrence from mid(X).  e = y - yhat solves
    L e = -r, r = L yhat - aE (computed in balls), e_j (j<R) bounded by rad(X_j).  On |t| = sigma' = SIGE*delta:
    |a_i| <= M_i sigma'^-m, |r/A_R| <= M_r sigma'^-N.  Scalar majorant b_n (prefix-sum form) for n < N2, then
    |e_n| <= A_e s'^-n (s' = SE*delta) once F_hom(n0) <= 1/2 and A_e >= 2 M_r s'^R/(n0+1)_R.
    The step matrix T is enclosed this way from the exact identity (zero input radii); the input X (balls) is then mapped
    by the ball product T_encl * X, which propagates input radii through |T| rather than through the majorant."""
    Xin = X
    X = acb_mat(R + 1, R + 1)
    for i_ in range(R + 1): X[i_, i_] = acb(1)
    K = X.ncols()
    delta = arb(delta_of(c).lower())
    sig = delta * SIGE; sp = delta * SE
    S = arb(2.0 ** math.floor(math.log2(float(delta.mid()))))
    S2_ = arb(2.0 ** math.floor(math.log2(float(delta_of(c2).mid()))))
    h = c2 - c
    rho = abs(h) / sp
    assert float(rho.upper()) < 0.85, float(rho.upper())
    if Nt is None: Nt = int(TB / (-math.log2(SIGE))) + R + 20
    alpha = [[a_ * S ** (m - i) for m, a_ in enumerate(taylor_shift(A_int[i], c))] for i in range(R + 1)]
    ecs = [e_ * S ** m for m, e_ in enumerate(taylor_shift(Ec, c))]
    NE = Nt + RHO + 1
    V = {}
    for sft in range(-R, RHO + 1):
        acc = acb_poly([0]); nz = False
        for i in range(R + 1):
            m = i + sft
            if 0 <= m < len(alpha[i]):
                acc += FFa[i] * alpha[i][m]; nz = True
        if nz: V[sft] = acc.evaluate([acb(n - sft) for n in range(NE)])
    arow = acb_mat([[X[R, col] for col in range(K)]]).mid()
    rows = [acb_mat([[X[k, col] for col in range(K)]]).mid() for k in range(R)]
    radin = [max(float(X[k, col].rad().upper()) * 1.0000001 for k in range(R)) for col in range(K)]
    # (input radii are per entry; we use the max over k of rad/S^k-scaled... see below)
    for n in range(0, Nt - R):
        acc = None
        for sft, vv in V.items():
            if sft == -R: continue
            idx = n - sft
            if 0 <= idx < len(rows):
                t = rows[idx] * vv[n]
                acc = t if acc is None else acc + t
        acc = (arow * ecs[n] - acc) if n < len(ecs) else acc * acb(-1)
        rows.append((acc * (1 / V[-R][n])).mid())
    # residual r_n (scaled tau-units), n = 0 .. NE-1, in balls; bound sum |r_n| (sig/S)^n per column
    sigS = sig / S
    rsum = [arb(0)] * K
    pw = arb(1)
    for n in range(NE):
        acc = None
        for sft, vv in V.items():
            idx = n - sft
            if 0 <= idx < Nt:
                t = rows[idx] * vv[n]
                acc = t if acc is None else acc + t
        if n < len(ecs): acc = (acc - arow * ecs[n]) if acc is not None else arow * (-ecs[n])
        if acc is None:
            pw *= sigS; continue
        for col in range(K):
            rsum[col] += abs(acc[0, col]) * pw
        pw *= sigS
    M, Me, minAR = circle_bounds(c, sig)
    # minAR is a lower bound of |A_R(y)| (unscaled) on the circle; r(t) in tau units equals sum r_n tau^n with
    # A_R scaled by S^-R (alpha scaling): scaled A_R(c+S tau) S^-R  ->  M_r = rsum / (minAR * S^-R)
    Mr = [(rsum[col] / minAR).upper() for col in range(K)]
    # homogeneous majorant constants on sig (unscaled): F_hom(N) with s' = sp
    def Fh(N):
        acc_ = arb(0)
        for i in range(R):
            den = arb(1)
            for t in range(i + 1, R + 1): den *= (N + t)
            acc_ += M[i] * sp ** (R - i) / den
        return acc_ / (1 - sp / sig)
    n0 = 0
    while float(Fh(n0).upper()) > 0.5:
        n0 += 4
        assert n0 < 20000, ('n0', float(delta.mid()))
    MAXN0[0] = max(MAXN0[0], n0)
    N2 = max(n0 + R, 120)
    # scalar majorant in unscaled units (arb, 64 bits; all quantities nonnegative upper bounds): B_n = b_n sig^n
    lp = ctx.prec; ctx.prec = 64
    sigu = arb(sig.upper()); spl = arb(sp.lower()); hu = arb(abs(h).upper()); sigSu = arb(sigS.upper())
    Mu = [arb(x.upper()) for x in M]
    facR = []
    for N in range(N2):
        fr = arb(1)
        for t in range(1, R + 1): fr *= (N + t)
        facR.append(fr)
    rho_ = hu / spl
    out = acb_mat(R + 1, K)
    rr = S2_ / S
    errs = [[None] * K for _ in range(R)]
    for col in range(K):
        Bv = [arb(X[j, col].rad().upper()) * sigSu ** j for j in range(R)]
        pref = [arb(0)]
        for x in Bv: pref.append(pref[-1] + x)
        Mrc = arb(Mr[col])
        for N in range(N2):
            tot = Mrc * sigu ** R
            fi = arb(1)
            for i in range(R):
                if i > 0: fi *= (N + i)
                tot += Mu[i] * sigu ** (R - i) * fi * pref[N + i + 1]
            Bv.append(arb((tot / facR[N]).upper())); pref.append(pref[-1] + Bv[-1])
        Ae = arb(0)
        q_ = spl / sigu
        for j in range(n0 + R):
            v_ = (Bv[j] * q_ ** j).upper()
            if v_ > Ae: Ae = arb(v_)
        rr0 = arb(1)
        for t in range(1, R + 1): rr0 *= (n0 + t)
        v_ = (2 * Mrc * spl ** R / rr0).upper()
        if v_ > Ae: Ae = arb(v_)
        for k in range(R):
            err = arb(0)
            for n in range(k, N2):
                err += math.comb(n, k) * Bv[n] * sigu ** (-n) * hu ** (n - k)
            kap = rho_ * (N2 + 1) / (N2 + 1 - k)
            assert float(kap.upper()) < 0.97
            err += Ae * spl ** (-k) * math.comb(N2, k) * rho_ ** (N2 - k) / (1 - kap)
            errs[k][col] = (err * S2_ ** k).upper()
    ctx.prec = lp
    for k in range(R):
        for col in range(K):
            e_ = errs[k][col]
            out[k, col] = acb(arb(0, e_), arb(0, e_))
    # midpoint evaluation of yhat at h (balls track rounding)
    hS = h / S
    hp_ = [acb(1)]
    for n in range(1, Nt): hp_.append(hp_[-1] * hS)
    Ymat = acb_mat([[rows[n][0, col] for n in range(Nt)] for col in range(K)])
    Vh = acb_mat(Nt, R)
    for n in range(Nt):
        for k in range(min(n, R - 1) + 1): Vh[n, k] = math.comb(n, k) * hp_[n - k]
    Z = Ymat * Vh
    for k in range(R):
        f = rr ** k
        for col in range(K):
            out[k, col] = out[k, col] + Z[col, k] * f
    for col in range(K): out[R, col] = X[R, col]
    out = out * Xin
    NSTEP[0] += 1
    if DIAG: print('step %d c=%s delta=%.3g n0=%d maxerr=%s |M12 s|=%.3g' % (NSTEP[0], c.str(5), float(delta.mid()), n0, arbmax([out[k, col].rad() for k in range(R) for col in range(out.ncols())]).str(3), float((M[R-1]*sp).mid())), flush=True)
    return out

# ---------- paths (exact dyadic points), avoiding all singular points ----------
singf = [0j] + rootsf
def fdist(z): return min(abs(z - s_) for s_ in singf)
def dy(x, scale):
    e = max(0, int(-math.floor(math.log2(scale))) + 30)
    return Fraction(round(x * 2 ** e), 2 ** e)
def pt(z, scale=None): return (dy(z.real, scale or fdist(z)), dy(z.imag, scale or fdist(z)))
def to_acb(p): return acb(arb(fmpq(p[0].numerator, p[0].denominator)), arb(fmpq(p[1].numerator, p[1].denominator)))
def zc(p): return complex(float(p[0]), float(p[1]))
STEPF = FR * 0.95
def seg_points(p0, p1):
    z1 = zc(p1); pts = [p0]; z = zc(p0)
    while True:
        dd = fdist(z); L = abs(z1 - z)
        if L <= STEPF * dd: pts.append(p1); return pts
        z = z + (z1 - z) / L * STEPF * dd
        pts.append(pt(z)); z = zc(pts[-1])
def arc_points(p0, center, radius, ang0, ang1, pend):
    import cmath
    pts = [p0]; a = ang0; z = zc(p0)
    while True:
        dd = fdist(z)
        da = 2 * math.asin(min(1.0, STEPF * dd / (2 * radius))) * 0.97
        if abs(ang1 - a) <= da: pts.append(pend); return pts
        a = a + da * (1 if ang1 > a else -1)
        zz = center + radius * cmath.exp(1j * a)
        pts.append(pt(zz)); z = zc(pts[-1])

if __name__ == '__main__':
    import cmath
    yt = sorted([z.real for z in rootsf if abs(z.imag) < 1e-30 and z.real > 100], key=float)
    ap_real = [z.real for z in rootsf if abs(z.imag) < 1e-30]
    log('real singular points:', ['%.6g' % x for x in sorted(ap_real)])
    lam_ = sorted([r.real for r, m in fmpz_poly([P4[j][-1] for j in range(5)]).complex_roots()], key=lambda x: -float(x.mid()))
    ytrue = [1 / float(l.mid()) for l in lam_]          # y1 < y2 < y3 < y4
    def is_true(x): return any(abs(x - t) < 1e-6 * t for t in ytrue)
    y1, y2, y3, y4 = ytrue
    BASE = Fraction(round(0.42 * y2))
    pos = sorted([x for x in ap_real if x > 0])
    apos = [x for x in pos if not is_true(x)]
    r2 = (y3 - y2) / 3
    assert not any(y2 - r2 <= x <= y3 + r2 for x in apos), 'apparent point near y2, y3'
    assert min(abs(complex(float(BASE), 0) - s_) for s_ in singf) > 0.2 * float(BASE)
    bpt = (Fraction(1, 2 ** BEXP), Fraction(0)); base = (BASE, Fraction(0))
    def axis_path(p0, xend, aps):
        """real-axis path from p0 to xend, detouring above each apparent point in aps (sorted)"""
        path = [p0]
        for x in aps:
            rr_ = min(min(abs(complex(x, 0) - s_) for s_ in singf if abs(s_ - complex(x, 0)) > 0), x, abs(x - zc(p0).real)) / 2
            pa = pt(complex(x - rr_, 0)); pb = pt(complex(x + rr_, 0))
            path += seg_points(path[-1], pa)[1:]
            path += arc_points(pa, x, x - zc(pa).real, math.pi, 0.0, pb)[1:]
        path += seg_points(path[-1], xend)[1:]
        return path
    small_ap = [x for x in apos if float(bpt[0]) < x < float(BASE)]
    mid_ap = [x for x in apos if float(BASE) < x < y2 - r2]
    log('q=%d true y_k %s, base %s, apparent on (b,base) %s, on (base, y2) %s' % (q, ['%.6g' % t for t in ytrue], BASE, small_ap, mid_ap))
    path0 = axis_path(bpt, base, small_ap)
    P2a = pt(complex(y2 - r2, 0)); P2b = pt(complex(y2 + r2, 0)); P3a = pt(complex(y3 - r2, 0))
    S0path = axis_path(base, P2a, mid_ap)
    pieces = {
        'CA': arc_points(base, 0.0, float(BASE), 0.0, 2 * math.pi, base),
        'S0': S0path,
        'C2': arc_points(P2a, y2, r2, math.pi, 3 * math.pi, P2a),
        'U2': arc_points(P2a, y2, r2, math.pi, 0.0, P2b),
        'S23': seg_points(P2b, P3a),
        'C3': arc_points(P3a, y3, r2, math.pi, 3 * math.pi, P3a),
    }
    if EXTRA:
        a39 = [x for x in apos if x < float(bpt[0])][-1]
        r3b = r2; r4 = (y4 - y3) / 3
        P3b = pt(complex(y3 + r3b, 0)); P4a = pt(complex(y4 - r4, 0))
        pieces['U3'] = arc_points(P3a, y3, r3b, math.pi, 0.0, P3b)
        pieces['S34'] = seg_points(P3b, P4a)
        pieces['C4'] = arc_points(P4a, y4, r4, math.pi, 3 * math.pi, P4a)
        # loop around y1 alone, from b: real axis down to 3.9e-8 + r, detour above it, down to y1 + y1/2, circle
        ra_ = (a39 - y1) / 3
        pa_ = pt(complex(a39 + ra_, 0)); pb_ = pt(complex(a39 - ra_, 0)); P1a = pt(complex(1.5 * y1, 0))
        pieces['P0'] = path0
        pieces['S1'] = seg_points(bpt, pa_)[:-1] + arc_points(pa_, a39, ra_, 0.0, math.pi, pb_)[:-1] + seg_points(pb_, P1a)
        pieces['C1'] = arc_points(P1a, y1, 0.5 * y1, 0.0, 2 * math.pi, P1a)
    log('path0 points', len(path0), ' pieces', {k: len(v) for k, v in pieces.items()})
    # transport g from b to base
    v = acb_mat(R + 1, 1)
    Sb = arb(2.0 ** math.floor(math.log2(float(delta_of(to_acb(bpt)).mid()))))
    for k in range(R): v[k, 0] = init[k] * Sb ** k
    v[R, 0] = acb(1)
    for t in range(len(path0) - 1):
        v = step(v, to_acb(path0[t]), to_acb(path0[t + 1]))
    log('g at base: rad max %s, |u_0| = %s, max n0 so far %d' % (arbmax([v[k, 0].rad() for k in range(R)]).str(3), abs(v[0, 0]).str(8), MAXN0[0]))
    I = acb_mat(R + 1, R + 1)
    for i in range(R + 1): I[i, i] = acb(1)
    Tm = {}
    for nm, pl in pieces.items():
        X = I
        for t in range(len(pl) - 1):
            X = step(X, to_acb(pl[t]), to_acb(pl[t + 1]))
        Tm[nm] = X
        log('piece %s: steps %d, max rad %s, max n0 %d' % (nm, NSTEP[0], arbmax([X[i, j].rad() for i in range(R + 1) for j in range(R + 1)]).str(3), MAXN0[0]))
    S0i = Tm['S0'].inv(); U2i = Tm['U2'].inv(); S23i = Tm['S23'].inv()
    G = {'A': Tm['CA'], '2': S0i * Tm['C2'] * Tm['S0'], '3': S0i * U2i * S23i * Tm['C3'] * Tm['S23'] * Tm['U2'] * Tm['S0']}
    if EXTRA:
        U3i = Tm['U3'].inv(); S34i = Tm['S34'].inv(); P0 = Tm['P0']; P0i = P0.inv(); S1i = Tm['S1'].inv()
        G['4'] = S0i * U2i * S23i * U3i * S34i * Tm['C4'] * Tm['S34'] * Tm['U3'] * Tm['S23'] * Tm['U2'] * Tm['S0']
        G['1'] = P0 * S1i * Tm['C1'] * Tm['S1'] * P0i
    pickle.dump({'G': {k: [[(str(M_[i, j].real), str(M_[i, j].imag)) for j in range(R + 1)] for i in range(R + 1)] for k, M_ in G.items()},
                 'v': [(str(v[i, 0].real), str(v[i, 0].imag)) for i in range(R + 1)], 'q': q, 'mu': MU, 'prec': PREC},
                open(os.path.join(HERE, 'runs', 'mono13_q%d%s.pkl' % (q, TAG)), 'wb'))
    dv = G['A'] * v - v
    log('check gA v = v: max |diff| upper %s' % arbmax([abs(dv[i, 0]) for i in range(R + 1)]).str(3))
    # Krylov words: '', 2, 3, A2, A3, AA2, AA3, ... (14 vectors)
    words = ['', '2', '3']
    while len(words) < R + 1:
        words.append('A' + words[len(words) - 2])
    vecs = {}
    def apply(w):
        if w in vecs: return vecs[w]
        x = v if w == '' else G[w[0]] * apply(w[1:])
        vecs[w] = x; return x
    Mx = acb_mat(R + 1, R + 1)
    for j, w in enumerate(words):
        x = apply(w)
        for i in range(R + 1): Mx[i, j] = x[i, 0]
    det = Mx.det()
    log('words', words)
    # RF1-F4: this is a DIAGNOSTIC with the default word list (which contains the dependent word A^5 3); the certificate
    # is cert_det.py's greedy word choice, filed as runs/mono13_q<q><tag>_det.json and runs/cert_det_q<q>.txt.
    log('DIAGNOSTIC (default words, not the certificate; see cert_det.py): det =', det.str(10), ' contains 0:', det.contains(0))
    json.dump({'q': q, 'mu': MU, 'prec': PREC, 'words': words, 'det': det.str(20), 'det_excludes_0': not det.contains(0),
               'note': 'diagnostic (default words); the certificate is cert_det.py -> *_det.json',
               'steps': NSTEP[0], 'max_n0': MAXN0[0]}, open(os.path.join(HERE, 'runs', 'mono13_q%d%s.json' % (q, TAG)), 'w'))
