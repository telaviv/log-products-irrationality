"""budget_cert.py -- certified margins M = (1 - 1/m) L - tau, L = log(4 y2), y2 = 1/lambda2, lambda2 the second root of the exact
recurrence's characteristic polynomial (= the lambda-quartic, geom.py), tau = 6 (m^2-1)/m^2 (rows (0,6,...,6), e = 0; T = L for
the Koebe map phi(w) = 4 y2 w/(1+w)^2).  Arb.  Usage: python budget_cert.py q1 q2 ..."""
import sys, json
from flint import fmpz_poly, arb, ctx
ctx.prec = 200
for q in map(int, sys.argv[1:]):
    d = json.load(open('runs/rec_q%d.json' % q)); P = [[int(x) for x in pp] for pp in d['polys']]
    ch = fmpz_poly([P[j][-1] for j in range(5)])
    rs = [r for r, m in ch.complex_roots()]
    assert all(r.imag.contains(0) or abs(float(r.imag.mid())) < 1e-50 for r in rs)
    lam = sorted([r.real for r in rs], key=lambda x: -float(x.mid()))
    assert all(float(l.lower()) > 0 for l in lam)
    l2 = lam[1]
    L = arb(4).log() - l2.log()
    row = 'q=%d  lambda2=%s  log y2=%s  L=%s ' % (q, l2.str(10), (-l2.log()).str(10), L.str(10))
    for m in (8, 13, 14):
        tau = arb(6) * (m * m - 1) / (m * m)
        M = L * (1 - arb(1) / m) - tau
        row += ' M(m=%d)=%s' % (m, M.str(6))
    print(row)
