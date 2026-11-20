"""geom.py -- characteristic polynomial of the exact recurrence vs the lambda-quartic (R6), integer roots of p_4, contact check.
Usage: python geom.py q"""
import sys, json, os
import sympy as sp
from flint import fmpz_poly, arb, ctx
HERE = os.path.dirname(os.path.abspath(__file__))
q = int(sys.argv[1]); W = q*q
d = json.load(open(os.path.join(HERE, 'runs', 'rec_q%d.json' % q)))
P = [[int(x) for x in pp] for pp in d['polys']]
n, lam = sp.symbols('n lam')
pj = [sp.Poly(list(reversed(pp)), n) for pp in P]
char = sum(sp.Integer(P[j][-1]) * lam**j for j in range(5))
quart = lam**4 + (16384*W**3 - 20480*W**2 + 5248*W - 64)*lam**3 + (135168*W**2 + 159744*W + 1536)*lam**2 + (294912*W - 16384)*lam + 65536
for s in (1, -1):
    r = sp.cancel(char.subs(lam, s*lam) / quart)
    print('char(%+d*lam)/quartic =' % s, r)
print('lc p_j:', [P[j][-1] for j in range(5)])
for j in range(5):
    ir = [rt for rt in sp.Poly(pj[j]).ground_roots() if rt.is_integer]
    print('p_%d integer roots:' % j, ir, ' rational roots:', [rt for rt in sp.Poly(pj[j]).ground_roots()])
ctx.prec = 200
rs = fmpz_poly([P[j][-1] for j in range(5)]).complex_roots()
print('char roots (lambda):', [r.str(12) for r, m in rs])
