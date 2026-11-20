# independent budget: lambda-quartic (FILTER 1.3 / geom.py), Arb roots, check vs rec json char poly, M = (13/14)L - 6*195/196
import json, sys
from flint import fmpz_poly, arb, ctx, fmpq
ctx.prec = 300
P_DIR = sys.argv[1]
for q in range(20, 34):
    W = q*q
    quart = fmpz_poly([65536, 294912*W - 16384, 135168*W**2 + 159744*W + 1536, 16384*W**3 - 20480*W**2 + 5248*W - 64, 1])
    rs = [r for r, m in quart.complex_roots()]
    assert all(r.imag.contains(0) for r in rs) and len(rs) == 4
    lam = sorted([-r.real for r in rs], key=lambda x: -float(x.mid()))   # Lysov normalisation: positive
    assert all(float(x.lower()) > 0 for x in lam)
    y = [1/x for x in lam]
    L = (4*y[1]).log()
    tau = arb(6)*195/196
    M14 = L*arb(13)/14 - tau
    M13 = L*arb(12)/13 - arb(6)*168/169
    extra = ''
    try:
        d = json.load(open('%s/rec_q%d.json' % (P_DIR, q))); P = [[int(x) for x in pp] for pp in d['polys']]
        ch = fmpz_poly([P[j][-1] for j in range(5)])
        # ch(lam) should be c * quart(-lam)
        qm = fmpz_poly([quart[k]*(-1)**k for k in range(5)])
        c = fmpq(int(ch[4]), int(qm[4]))
        ok = all(fmpq(int(ch[k])) == c*int(qm[k]) for k in range(5))
        extra = ' char==c*quart(-lam): %s, deg34 lc all nonzero: %s' % (ok, all(P[j][-1] != 0 for j in range(5)))
    except FileNotFoundError: pass
    print('q=%d y1=%s y2=%s y3/y2=%s L=%s M14=%s M13=%s  y1<2^-20:%s%s' % (q, y[0].str(4), y[1].str(10), (y[2]/y[1]).str(4), L.str(12), M14.str(8), M13.str(5), float(y[0].upper()) < 2**-20, extra))
