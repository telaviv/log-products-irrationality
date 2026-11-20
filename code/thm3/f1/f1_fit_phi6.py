"""Recover Phi6(a, b) from the degree-6 factors f_{b0}(a) of det(a, b0) (runs/det_b*.json).

Unknowns: c_{ij} (0 <= i, j <= DB) and one scale lambda_{b0} per b0;  sum_j c_{ij} b0^j = lambda_{b0} [a^i] f_{b0}.
The integer kernel (expected 1-dimensional) gives Phi6 up to sign; written to runs/phi6.txt.
"""
import json, glob, sys
from flint import fmpz, fmpz_mat, fmpq
from lib_f1 import CTX, A_, B_

DB = 6
data = []
for fn in sorted(glob.glob('runs/det_b*.json')):
    d = json.load(open(fn))
    f6 = [f for f in d['factors'] if f['deg'] == 6]
    if len(f6) != 1:
        print(fn, 'number of degree-6 factors', len(f6)); continue
    data.append((d['b0'], [int(c) for c in f6[0]['coeffs']]))
print('b0 used', [b for b, _ in data])
nc = 7 * (DB + 1)
nl = len(data)
rows = []
for k, (b0, cf) in enumerate(data):
    for i in range(7):
        row = [0] * (nc + nl)
        for j in range(DB + 1):
            row[i * (DB + 1) + j] = b0 ** j
        row[nc + k] = -cf[i]
        rows.append(row)
if '--sym' in sys.argv:  # impose c_ij = c_ji (CDT: Phi_d(a,b) = Phi_d(b,a))
    for i in range(7):
        for j in range(i + 1, 7):
            row = [0] * (nc + nl); row[i * 7 + j] = 1; row[j * 7 + i] = -1; rows.append(row)
M = fmpz_mat(rows)
X, nullity = M.nullspace()
print('nullity', nullity)
assert nullity == 1
v = [X[r, 0] for r in range(nc + nl)]
g = fmpz(0)
for x in v:
    g = g.gcd(x)
v = [x // g for x in v]
P = CTX.from_dict({})
for i in range(7):
    for j in range(DB + 1):
        c = v[i * (DB + 1) + j]
        if c != 0:
            P += c * A_ ** i * B_ ** j
if P.coeffs()[0] < 0:
    P = -P
print('Phi6 =', P)
Ps = CTX.from_dict({(eb, ea, ey): c for (ea, eb, ey), c in zip(P.monoms(), P.coeffs())})
print('symmetric:', Ps == P)
open('runs/phi6.txt', 'w').write(str(P) + '\n')
