# Independent referee check of RESULT_4T sec. 2-3: local exponent at y_k via explicit theta->D conversion (numeric roots,
# 120 digits), exact symbolic identity in a second formulation, p0 roots (flint factor), Casoratian at z=q (sympy series).
import json, sys
import sympy as sp, mpmath as mp
from flint import fmpz_poly
import os; F = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..")  # logprod copy: was the absolute path of pilot_F1
d = json.load(open(os.path.join(F, "runs", "i1_L4_symbolic.json")))
n, w, y = sp.symbols('n w y')
P = [sum(sp.Integer(int(v)) * n**int(k.split(',')[0]) * w**int(k.split(',')[1]) for k, v in pj.items()) for pj in d['p']]
# theta-operator L = sum_j y^(4-j) p_j(theta - j) = sum_k c_k(y) theta^k ; c_k exactly via expansion of p_j(t - j)
t = sp.symbols('t')
ck = {}
for j in range(5):
    pt = sp.Poly(sp.expand(P[j].subs(n, t - j)), t)
    for (k,), c in pt.terms():
        ck[k] = ck.get(k, 0) + y**(4 - j) * c
c34 = sp.expand(ck[34]); c33 = sp.expand(ck[33])
# D-form: theta^k = sum_i S(k,i) y^i D^i ; a_34 = y^34 c34, a_33 = y^33 (c33 + S(34,33) c34)
from sympy.functions.combinatorial.numbers import stirling
a34 = sp.expand(y**34 * c34); a33 = sp.expand(y**33 * (c33 + stirling(34, 33) * c34))
# second formulation: 2 a33 - 65 a34'  vanishes at roots of c34  <=> prem by c34 == 0
expr = sp.expand(2 * a33 - 65 * sp.diff(a34, y))
pr = sp.prem(sp.Poly(expr, y), sp.Poly(c34, y))
print('prem(2 a33 - 65 a34\', c34) == 0 (symbolic in w):', pr.is_zero)
# discriminant of c34 in y (distinct roots)
disc = sp.factor(sp.discriminant(sp.Poly(c34, y)))
print('disc_y(c34) =', str(disc)[:300])
mp.mp.dps = 120
for q in (26, 33, 40, 53):
    c34q = sp.Poly(c34.subs(w, q * q), y); a34q = sp.Poly(a34.subs(w, q*q), y); a33q = sp.Poly(a33.subs(w, q*q), y)
    rts = mp.polyroots([int(c) for c in c34q.all_coeffs()], maxsteps=400, extraprec=400)
    da34 = a34q.diff(y)
    betas = []
    for r in rts:
        v = mp.polyval([int(c) for c in a33q.all_coeffs()], r) / mp.polyval([int(c) for c in da34.all_coeffs()], r)
        betas.append(33 - v)
    print('q', q, 'y_k', [mp.nstr(mp.re(r), 8) for r in rts], 'beta_k - 1/2:', [mp.nstr(abs(b - mp.mpf(1)/2), 3) for b in betas])
# p0(.,q^2), p4(.,q^2): nonnegative integer roots via flint factorisation, q = 26..53
bad = []
for q in range(26, 54):
    for jj in (0, 4):
        pq = sp.Poly(P[jj].subs(w, q * q), n)
        fp = fmpz_poly([int(c) for c in reversed(pq.all_coeffs())])
        for fac, m in fp.factor()[1]:
            if fac.degree() == 1:
                a0, a1 = int(fac[0]), int(fac[1])
                if a0 % a1 == 0 and -a0 // a1 >= 0: bad.append((q, jj, -a0 // a1))
print('nonneg integer roots of p0/p4 at q=26..53:', bad)
# Casoratian at z = q directly: v_n = (Q_n, polypart(Q_n f_j)) using sympy series of f_j
z, u = sp.symbols('z u')
from math import comb
def Qn(nn):
    return sum((-1)**(nn + l) * comb(nn + l, l) * sum(comb(nn, l - k) * comb(2*nn + k, 2*k) * comb(nn + k, k) * comb(4*nn + 2*k, 2*nn + k)
               for k in range(max(0, l - nn), min(2*nn, l) + 1)) * z**(2*l) for l in range(3*nn + 1))
fs = [sp.log(1 + u), sp.log(1 - u), sp.log(1 + u) * sp.log(1 - u)]
ser = [sp.series(f, u, 0, 22).removeO() for f in fs]
rows = []
for nn in range(4):
    Q = sp.expand(Qn(nn))
    row = [Q]
    for s in ser:
        prod = sp.expand(Q * s.subs(u, 1 / z))
        row.append(sum(c * z**e for (e,), c in sp.Poly(prod * z**30, z).terms() if e >= 30 for _ in [0]) if False else
                   sp.expand(sum(term for term in sp.Add.make_args(prod) if sp.degree(term, z) >= 0 and not term.has(1/z) and sp.Poly(term*z**40, z).degree() >= 40)))
    rows.append(row)
C0 = sp.factor(sp.Matrix(rows).det())
print('Casoratian C(0)(z) =', C0)
vals = {q: C0.subs(z, q) for q in range(26, 54)}
print('C(0)(q) != 0 for q=26..53:', all(v != 0 for v in vals.values()))
