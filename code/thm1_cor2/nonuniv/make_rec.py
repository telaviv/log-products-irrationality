"""make_rec.py -- runs/rec_q<q>.json by SPECIALISING the proved symbolic order-4 recurrence of PROOF_I1 (pilot F1,
runs/i1_L4_symbolic.json) at w = q^2 (Theorem I1, symbolic in q).  Same format as F1's rec.py output (polys[j][e] = coefficient
of n^e in p_j(n, q^2)), made primitive (a nonzero constant factor does not change the recurrence).  Sanity: for q = 26 the
result must equal F1's runs/rec_q26.json exactly (PROOF_I1 Sec. 3 says the specialisation is identical, factor 1).
Usage: python make_rec.py q [q ...]"""
import json, os, sys
from math import gcd
HERE = os.path.dirname(os.path.abspath(__file__))
F1 = os.path.join(HERE, '..')   # = code/thm1_cor2 (was the working folder of the q >= 26 proof)
d = json.load(open(os.path.join(F1, 'runs', 'i1_L4_symbolic.json')))
for q in map(int, sys.argv[1:]):
    W = q * q
    polys = []
    for pj in d['p']:
        co = [0] * 35
        for key, val in pj.items():
            a, b = map(int, key.split(','))
            co[a] += int(val) * W ** b
        polys.append(co)
    g = 0
    for co in polys:
        for c in co: g = gcd(g, c)
    polys = [[c // g for c in co] for co in polys]
    out = dict(q=q, order=4, deg=34, polys=[[str(c) for c in co] for co in polys], source='specialisation of PROOF_I1 L4',
               content_removed=str(g))
    json.dump(out, open(os.path.join(HERE, 'runs', 'rec_q%d.json' % q), 'w'))
    if q == 26:
        ref = json.load(open(os.path.join(F1, 'runs', 'rec_q26.json')))
        same = [[int(x) for x in pp] for pp in ref['polys']] == polys
        print('q=26 identical to F1 rec_q26.json:', same)
    print('q=%d written, removed content %d' % (q, g))
