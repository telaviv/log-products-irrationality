"""det(a, b0) of the Lemma 14.3.1 system as an exact polynomial in a, for fixed integer b0.

D(a) := det * (a - b0)^E * 2^F is interpolated from exact values at integer a (Newton), checked on extra points,
and factored over Z.  Usage: python f1_det_univ.py b0 [npts]   ->  runs/det_b{b0}.json
"""
import sys, json, os, time
from flint import fmpq, fmpz, fmpz_poly, fmpq_poly
from lib_det import system

b0 = int(sys.argv[1])
NPTS = int(sys.argv[2]) if len(sys.argv) > 2 else 420
E = 86
NCHECK = 6
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'runs', 'det_b%d.json' % b0)
vals = {}
if os.path.exists(out + '.part'):
    vals = {int(k): fmpq(*map(int, v.split('/'))) if '/' in v else fmpq(int(v)) for k, v in json.load(open(out + '.part')).items()}

pts = []
k = 1
while len(pts) < NPTS + NCHECK:
    for a in (b0 + k, b0 - k):
        if a not in (0, 1, -1, b0, -b0) and len(pts) < NPTS + NCHECK:
            pts.append(a)
    k += 1
t0 = time.time()
for i, a in enumerate(pts):
    if a in vals:
        continue
    M, _, _ = system(a, b0)
    vals[a] = M.det() * fmpq(a - b0) ** E
    if i % 20 == 0:
        json.dump({str(k): str(v) for k, v in vals.items()}, open(out + '.part', 'w'))
        print(i, a, '%.0fs' % (time.time() - t0), flush=True)
json.dump({str(k): str(v) for k, v in vals.items()}, open(out + '.part', 'w'))
# common denominator must be a power of 2
den = fmpz(1)
for a in pts:
    d = vals[a].denom()
    den = den * d // den.gcd(d)
F = 0
while den % 2 == 0:
    den //= 2
    F += 1
assert den == 1, ('denominator not a power of 2', den)
xs = pts[:NPTS]
ys = [vals[a] * fmpq(2) ** F for a in xs]
# Newton interpolation
n = len(xs)
c = list(ys)
for j in range(1, n):
    for i in range(n - 1, j - 1, -1):
        c[i] = (c[i] - c[i - 1]) / (xs[i] - xs[i - j])
P = fmpq_poly([c[-1]])
X = fmpq_poly([0, 1])
for i in range(n - 2, -1, -1):
    P = P * (X - xs[i]) + c[i]
ok = all(P(fmpq(a)) == vals[a] * fmpq(2) ** F for a in pts[NPTS:])
print('degree', P.degree(), 'checks', ok, 'F', F, flush=True)
assert ok and P.degree() < NPTS - 10
Pz = P.numer()  # P = Pz / denom(P)
cont, fl = Pz.factor()
res = {'b0': b0, 'denomP': str(P.denom()), 'E': E, 'F': F, 'degree': P.degree(), 'content': str(cont),
       'factors': [{'deg': f.degree(), 'mult': e, 'poly': str(f) if f.degree() <= 30 else None,
                    'coeffs': [str(x) for x in f.coeffs()]} for f, e in fl]}
json.dump(res, open(out, 'w'), indent=1)
for f, e in fl:
    print('  deg %d mult %d %s' % (f.degree(), e, str(f) if f.degree() <= 8 else ''))
