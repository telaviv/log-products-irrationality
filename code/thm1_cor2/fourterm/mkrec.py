"""mkrec.py -- runs/rec_q<q>.json for given q, by specialising the PROVED symbolic order-4 operator (PROOF_I1 Prop. 3,
../runs/i1_L4_symbolic.json) at w = q^2 (made primitive).  Same format as ../rec.py output (p_j coefficient lists in n).
Sanity: the recurrence is checked exactly on the four sequences for n <= 40 (not needed: Prop. 3/5 prove it).
Usage: python mkrec.py q1 q2 ..."""
import json, os, sys, math
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
d = json.load(open(os.path.join(HERE, '..', 'runs', 'i1_L4_symbolic.json')))
for q in map(int, sys.argv[1:]):
    W = q * q
    polys = []
    for pj in d['p']:
        c = [0] * 35
        for k, v in pj.items():
            e, f = map(int, k.split(','))
            c[e] += int(v) * W ** f
        polys.append(c)
    g = 0
    for pp in polys:
        for x in pp: g = math.gcd(g, x)
    polys = [[x // g for x in pp] for pp in polys]
    sys.argv = ['rec.py', str(q)]
    import importlib, rec
    importlib.reload(rec)
    vals = rec.exact_values(44)
    bad = 0
    for s in range(4):
        for n in range(41):
            tot = sum(sum(polys[j][e] * n ** e for e in range(35)) * Fraction(vals[n + j][s]) for j in range(5))
            bad += tot != 0
    json.dump(dict(q=q, order=4, deg=34, polys=[[str(x) for x in pp] for pp in polys], verified_to=40, failures=bad,
                   source='specialisation of i1_L4_symbolic.json at w=q^2, primitive'),
              open(os.path.join(HERE, 'runs', 'rec_q%d.json' % q), 'w'))
    print('q', q, 'content removed', g, 'sanity failures', bad, flush=True)
