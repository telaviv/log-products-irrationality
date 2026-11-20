"""exponents.py -- local exponents of the order-34 operator L = sum_j y^(4-j) p_j(theta - j) (p_j = the symbolic order-4
recurrence of PROOF_I1 Prop. 3, runs/i1_L4_symbolic.json, w = z^2) at its finite nonzero singular points y_k (roots of
c34(y) = sum_j y^(4-j) lc_n(p_j)), symbolically in w.

At a simple root y_k of c34 the theta^34-coefficient, L in d/dy form has a_34 = y^34 c34(y) with a simple zero, so y_k is a
regular singular point with exponents 0, 1, ..., 32 and beta_k = 33 - c33(y_k) / (y_k c34'(y_k)), where c33 is the
theta^33-coefficient: c33 = sum_j y^(4-j) (p_{j,33} - 34 j lc(p_j)).  (theta^N = y^N D^N + C(N,2) y^(N-1) D^(N-1) + ...; the
C(34,2) c34 term vanishes at y_k.)

Exact output: the polynomial Res_y(c34(y), (33 - b) y c34'(y) - c33(y)) in b, w (sympy), its factorisation, and whether
beta is constant in w.  Also: the nonnegative integer roots of p_0(n, q^2), p_4(n, q^2) for q = 2..60, and the Casoratian
C(0)(z) = det[v_0, v_1, v_2, v_3], v_n = (Q_n, P1_n, P2_n, P3_n)(z), as a polynomial in z (rec.py exact values).
Writes exponents.json.  Usage: python exponents.py   (from inside fourterm/)"""
import json, os, sys, time
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__)); UP = os.path.dirname(HERE)
T0 = time.time()
def log(*a): print('[%6.1fs]' % (time.time() - T0), *a, flush=True)
d = json.load(open(os.path.join(UP, 'runs', 'i1_L4_symbolic.json')))
w, n, y, b = sp.symbols('w n y b')
P = []
for pj in d['p']:
    P.append(sum(sp.Integer(int(v)) * n ** int(k.split(',')[0]) * w ** int(k.split(',')[1]) for k, v in pj.items()))
Pn = [sp.Poly(p, n) for p in P]
lc = [Pn[j].coeff_monomial(n ** 34) for j in range(5)]
c33j = [Pn[j].coeff_monomial(n ** 33) for j in range(5)]
c34 = sp.expand(sum(y ** (4 - j) * lc[j] for j in range(5)))
c33 = sp.expand(sum(y ** (4 - j) * (c33j[j] - 34 * j * lc[j]) for j in range(5)))
log('c34 degree in y:', sp.degree(c34, y), ' lc_j w-degrees:', [sp.degree(l, w) for l in lc])
out = {}
# beta as a function on the roots of c34
expr = sp.expand((33 - b) * y * sp.diff(c34, y) - c33)
res = sp.resultant(c34, expr, y)
res = sp.factor(res)
log('Res_y factorised:', str(res)[:2000])
out['resultant_factored'] = str(res)
# direct check: is c33 == (33 - beta0) y c34' mod c34 for a constant beta0?
rem_form = sp.rem(sp.Poly(c33, y), sp.Poly(c34, y))
rem_der = sp.rem(sp.Poly(sp.expand(y * sp.diff(c34, y)), y), sp.Poly(c34, y))
ratio = None
cs1 = rem_form.all_coeffs(); cs2 = rem_der.all_coeffs()
if len(cs1) == len(cs2):
    rats = set(sp.simplify(a / bb) for a, bb in zip(cs1, cs2) if bb != 0)
    log('c33 / (y c34\') mod c34, coefficient ratios:', rats)
    if len(rats) == 1:
        ratio = rats.pop(); log('beta = 33 - ratio =', sp.nsimplify(33 - ratio))
        out['beta_constant'] = str(33 - ratio)
# root-free certificate: pseudo-remainder of 2 c33 - 65 y c34' by c34 in Z[w][y] is identically 0
prem = sp.prem(sp.Poly(sp.expand(2 * c33 - 65 * y * sp.diff(c34, y)), y), sp.Poly(c34, y))
log('prem(2 c33 - 65 y c34\', c34) == 0:', prem.is_zero)
out['prem_zero'] = bool(prem.is_zero)
lc0 = sp.Poly(lc[0], w); lc4 = sp.Poly(lc[4], w)
log('lc_n(p_0)(w) =', sp.factor(lc[0]), ';  lc_n(p_4)(w) =', sp.factor(lc[4]))
out['lc_p0_factored'] = str(sp.factor(lc[0])); out['lc_p4_factored'] = str(sp.factor(lc[4]))
# p_0, p_4 nonnegative integer roots, q = 2..60
roots = {}
for q in range(2, 61):
    rr = {}
    for j in (0, 4):
        pq = sp.Poly(P[j].subs(w, q * q), n)
        rr['p%d' % j] = [int(r) for r in pq.ground_roots() if r.is_integer and r >= 0]
    roots[q] = rr
bad = {q: r for q, r in roots.items() if r['p0'] or r['p4']}
log('q in 2..60 with a nonnegative integer root of p_0 or p_4:', bad)
out['p0_p4_nonneg_int_roots_q2_60'] = bad
# Casoratian C(0)(z)
from math import comb
z = sp.symbols('z')
def qcoefs(nn):   # identical to rec.py qcoefs: Q_n(z) = sum cs[i] z^i
    cs = [0] * (6 * nn + 1)
    for l in range(3 * nn + 1):
        s = 0
        for k in range(max(0, l - nn), min(2 * nn, l) + 1):
            s += comb(nn, l - k) * comb(2 * nn + k, 2 * k) * comb(nn + k, k) * comb(4 * nn + 2 * k, 2 * nn + k)
        cs[2 * l] = (-1) ** (nn + l) * comb(nn + l, l) * s
    return cs
K = 30
a1 = [sp.Integer(0)] + [sp.Rational((-1) ** (k + 1), k) for k in range(1, K + 1)]     # log(1+u), u = 1/z
a2 = [sp.Integer(0)] + [sp.Rational(-1, k) for k in range(1, K + 1)]                  # log(1-u)
a3 = [sp.Integer(0)] * (K + 1)
for i in range(1, K + 1):
    for j in range(1, K + 1 - i): a3[i + j] += a1[i] * a2[j]
def polypart(cs, a):   # polynomial part of (sum cs[i] z^i) * (sum_k a[k] z^-k)
    return sp.expand(sum(cs[i] * a[k] * z ** (i - k) for i in range(len(cs)) for k in range(1, i + 1) if cs[i]))
rows = []
for nn in range(4):
    cs = qcoefs(nn)
    rows.append([sp.expand(sum(cs[i] * z ** i for i in range(len(cs))))] + [polypart(cs, a) for a in (a1, a2, a3)])
C0 = sp.factor(sp.Matrix(rows).det())
log('Casoratian C(0)(z) =', C0)
out['casoratian_C0'] = str(C0)
C0p = sp.Poly(sp.numer(sp.together(C0)), z)
out['casoratian_C0_integer_roots'] = [int(r) for r in C0p.ground_roots() if r.is_integer]
log('integer roots of C(0):', out['casoratian_C0_integer_roots'])
json.dump(out, open(os.path.join(HERE, 'exponents.json'), 'w'), indent=1)
log('done')
