"""i1_extras.py -- facts about the symbolic L4 = sum p_i(n,w) S^i of runs/i1_L4_symbolic.json (w = z^2):
factorisations of p_0, p_4; the characteristic polynomial sum_j lc_n(p_j)(w) lambda^j versus the pilot's lambda-quartic
(geom.py), now symbolically in w = W.  Usage: python i1_extras.py"""
import json, os
import sympy as sp
HERE = os.path.dirname(os.path.abspath(__file__))
n, w, lam = sp.symbols('n w lam')
L4 = json.load(open(os.path.join(HERE, 'runs', 'i1_L4_symbolic.json')))
p = [sum(sp.Rational(v) * n ** int(key.split(',')[0]) * w ** int(key.split(',')[1]) for key, v in pi.items()) for pi in L4['p']]
for i in (0, 4):
    f = sp.factor_list(p[i])
    print('p_%d factors:' % i, [(sp.Poly(g, n, w).total_degree(), sp.degree(g, n), sp.degree(g, w), m) for g, m in f[1]])
    print('   linear-in-n factors:', [g for g, m in f[1] if sp.degree(g, n) == 1 and sp.degree(g, w) == 0])
lc = [sp.Poly(pi, n).LC() for pi in p]
char = sp.expand(sum(lc[j] * lam ** j for j in range(5)))
W = w
quart = lam**4 + (16384*W**3 - 20480*W**2 + 5248*W - 64)*lam**3 + (135168*W**2 + 159744*W + 1536)*lam**2 + (294912*W - 16384)*lam + 65536
r = sp.factor(sp.cancel(char / quart.subs(lam, -lam)))
print('char(lam) / quartic(-lam) =', r)
