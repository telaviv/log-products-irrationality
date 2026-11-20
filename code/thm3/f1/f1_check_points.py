"""Extra exact checks of det(a,b) (a-b)^86 = kappa * a^6 b^6 (a+b)^24 (...)^4 (...)^2 Phi4 Phi6 Phi14^3 Phi26 Phi79^3
at 40 random points (mixed parity and sign) and at 8 odd pairs in the Theorem C' range.  Output runs/check_points.json."""
import random, json
from flint import fmpq
from f1_phi6_grid import r
from lib_f1 import A_, B_
P = eval(open('runs/phi6.txt').read().replace('^', '**'), {'a': A_, 'b': B_})
T = [(m[0], m[1], int(c)) for m, c in zip(P.monoms(), P.coeffs())]
phi6 = lambda a, b: sum(c * fmpq(a) ** i * fmpq(b) ** j for i, j, c in T)
kappa = fmpq(645700815, 37778931862957161709568)
random.seed(20261119)
pts = [(random.randint(-10 ** 4, 10 ** 4), random.randint(-10 ** 4, 10 ** 4)) for _ in range(40)]
pts = [(x, y) for x, y in pts if x not in (0, 1, -1) and y not in (0, 1, -1) and x != y and x != -y]
odd = [(2 * m + 1, 2 * n + 1) for m, n in [(666664, 666666), (999997, 1000000), (1999997, 1999999), (-666667, -666665),
                                           (333332, 333333), (2999996, 2999999), (5 * 10 ** 6, 5 * 10 ** 6 + 14), (10 ** 9, 10 ** 9 + 2999)]]
res = []
for x, y in pts + odd:
    v = r(x, y)
    res.append(dict(a=x, b=y, ok=bool(v == kappa * phi6(x, y)), det_nonzero=bool(v != 0)))
print('all agree:', all(t['ok'] for t in res), ' points:', len(res), ' det != 0 at the odd C\' pairs:',
      all(t['det_nonzero'] for t in res[-len(odd):]))
json.dump(res, open('runs/check_points.json', 'w'), indent=1)
