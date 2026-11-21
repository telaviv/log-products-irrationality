"""unit1_exponents.py -- Unit 1 (hardest piece): local exponents of the order-34 operator L at its finite nonzero singular
points y_1..y_4, symbolically in w = q^2, from PROOF_I1's symbolic order-4 recurrence (runs/i1_L4_symbolic.json of pilot F1).

L = sum_{j=0}^{4} y^(4-j) p_j(theta - j, w),  theta = y d/dy,  L H = E (E polynomial)  [PROOF_I1 Thm I1, Sec. 5].
Write L = sum_k c_k(y) theta^k (k <= 34).  In d/dy form, a_34 = y^34 c_34, a_33 = y^33 (C(34,2) c_34 + c_33).
At a simple root y* != 0 of c_34 the indicial polynomial is (k)_33 [ (k-33) a_34'(y*) + a_33(y*) ], so the exponents are
0, 1, ..., 32 and rho = 33 - c_33(y*) / (y* c_34'(y*)).

Exact steps (sympy over Q(w)):
 (a) c_34(y) = sum_j y^(4-j) lc_n(p_j)  (lc_n = coefficient of n^34),  c_33(y) = sum_j y^(4-j) [n^33] p_j(n - j, w);
 (b) factor c_34 in Q[w][y]: content (a polynomial in w) times a primitive part; show the primitive part is (up to a
     constant) the y-quartic y^4 quart(-1/y) of FILTER / PROOF_I1 Sec. 7 and that its discriminant has no root w >= 4;
 (c) the polynomial whose roots are rho(y_i):  Res_y( c34prim(y), c_33(y) - (33 - r) y c_34'(y) )  as a polynomial in r over
     Q(w), factored.  If it is (const)(r - r_1)...(r - r_4) with constants r_i in Q, each root y_i has a w-independent
     exponent; the assignment r <-> y_i is then made by an Arb evaluation at each q.
Output: runs/unit1_exponents.json and stdout.
"""
import json, os, sys, time
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
F1 = os.path.join(HERE, '..')   # = code/thm1_cor2 (was the working folder of the q >= 26 proof)
d = json.load(open(os.path.join(F1, 'runs', 'i1_L4_symbolic.json')))
n, w, y, r = sp.symbols('n w y r')
P = []
for pj in d['p']:
    e = 0
    for key, val in pj.items():
        a, b = map(int, key.split(','))
        e += sp.Integer(val) * n**a * w**b
    P.append(sp.Poly(e, n, w))
t0 = time.time()
c34 = 0; c33 = 0
for j, pj in enumerate(P):
    shifted = sp.Poly(sp.expand(pj.as_expr().subs(n, n - j)), n)
    lc = shifted.coeff_monomial(n**34)
    sub = shifted.coeff_monomial(n**33)
    c34 += y**(4 - j) * lc
    c33 += y**(4 - j) * sub
c34 = sp.Poly(sp.expand(c34), y); c33 = sp.Poly(sp.expand(c33), y)
print('deg_y c34 =', c34.degree(), ' deg_y c33 =', c33.degree(), ' [%.1fs]' % (time.time() - t0))
cont, prim = c34.primitive() if False else (None, None)
# content in w
coeffs = [sp.factor(c) for c in c34.all_coeffs()]
g = coeffs[0]
for c in coeffs[1:]:
    g = sp.gcd(g, c)
g = sp.factor(g)
print('content of c34 (in w):', g)
c34p = sp.Poly(sp.expand(sp.cancel(c34.as_expr() / g)), y)
W = w
quart_lam = lambda lam: lam**4 + (16384*W**3 - 20480*W**2 + 5248*W - 64)*lam**3 + (135168*W**2 + 159744*W + 1536)*lam**2 + (294912*W - 16384)*lam + 65536
yq = sp.Poly(sp.expand(y**4 * quart_lam(-1 / y)), y)       # roots y = -1/lam
ratio = sp.cancel(c34p.as_expr() / yq.as_expr())
print('c34prim / (y^4 quart(-1/y)) =', ratio)
disc = sp.factor(sp.discriminant(c34p.as_expr(), y))
print('disc_y(c34prim) =', disc)
# (c) resultant in y
expr = c33.as_expr() - (33 - r) * y * sp.diff(c34.as_expr(), y)
res = sp.resultant(c34p.as_expr(), expr, y)
resf = sp.factor(res)
print('Res_y = ', resf, ' [%.1fs]' % (time.time() - t0))
# also check: g(w) has no root at w = q^2 >= 4 and y = 0 is not a root of c34prim
print('c34prim(0) =', sp.factor(c34p.eval(0)))
out = dict(content=str(g), ratio_to_quartic=str(ratio), disc=str(disc), resultant=str(resf),
           c34prim=str(c34p.as_expr()), c33=str(c33.as_expr()))
json.dump(out, open(os.path.join(HERE, 'runs', 'unit1_exponents.json'), 'w'), indent=1)
