"""findop.py -- exact minimal-order operator (recurrence order RHO, n-degree DG) common to Q, P1, P2, P3 at q, from the exact
order-4 recurrence: residues mod 62-bit primes, nullity-1 nullspace, CRT + rational reconstruction until stable, then exact
verification over Q on the exact sequences for n <= NV (all four sequences).  Output runs/op_q<q>_<RHO>_<DG>.json.
Usage: python findop.py q RHO DG NTERMS NV"""
import sys, json, math, time
from fractions import Fraction
from flint import nmod_mat
from sympy import prevprime
args = sys.argv; sys.argv = ['rec.py', args[1]]
q, RHO, DG, NT, NV = map(int, args[1:6])
from rec import exact_values, ratrec
d = json.load(open('runs/rec_q%d.json' % q)); P = [[int(x) for x in pp] for pp in d['polys']]
t0 = time.time()
vals = exact_values(3)
def run_mod(p):
    def red(x):
        x = Fraction(x); return x.numerator % p * pow(x.denominator % p, -1, p) % p
    S = [[red(v[s]) for v in vals] for s in range(4)]
    for n in range(0, NT - 3):
        pj = [sum(P[j][e] * pow(n, e, p) for e in range(35)) % p for j in range(5)]
        inv = pow(pj[4], -1, p)
        for X in S: X.append((-sum(pj[j] * X[n + j] for j in range(4)) * inv) % p)
    unk = (RHO + 1) * (DG + 1); rows = []
    for X in (S[0], S[3]):
        for n in range(NT - RHO):
            pw = [1]
            for e in range(DG): pw.append(pw[-1] * n % p)
            rows.append([pw[e] * X[n + j] % p for j in range(RHO + 1) for e in range(DG + 1)])
    M = nmod_mat(rows, p)
    Xn, nul = M.nullspace()
    if nul != 1: return None, nul
    return [int(Xn[i, 0]) for i in range(unk)], 1
p = 2 ** 62; acc = None; M_ = 1; prev = None; piv = None; np_ = 0
while True:
    p = prevprime(p)
    v, nul = run_mod(p)
    if v is None: print('prime', p, 'nullity', nul); continue
    if piv is None: piv = max(i for i in range(len(v)) if v[i])
    inv = pow(v[piv], -1, p); v = [x * inv % p for x in v]
    if acc is None: acc, M_ = v, p
    else:
        acc = [a + M_ * (((b - a) * pow(M_, -1, p)) % p) for a, b in zip(acc, v)]; M_ *= p
    np_ += 1
    rec = [ratrec(a, M_) for a in acc]
    if all(r is not None for r in rec):
        if rec == prev: break
        prev = rec
    else: prev = None
    print('primes', np_, '(%.0fs)' % (time.time() - t0), flush=True)
den = 1
for r in rec: den = den * r.denominator // math.gcd(den, r.denominator)
ints = [int(r * den) for r in rec]; g = 0
for x in ints: g = math.gcd(g, x)
ints = [x // g for x in ints]
polys = [[ints[j * (DG + 1) + e] for e in range(DG + 1)] for j in range(RHO + 1)]
ex = exact_values(NV + RHO)
bad = 0
for s in range(4):
    for n in range(NV + 1):
        tot = sum(sum(polys[j][e] * n ** e for e in range(DG + 1)) * Fraction(ex[n + j][s]) for j in range(RHO + 1))
        if tot != 0: bad += 1
print('primes used', np_, 'max digits', len(str(max(abs(x) for x in ints))), 'exact check n<=%d, 4 seqs: failures %d' % (NV, bad))
json.dump(dict(q=q, rho=RHO, dg=DG, polys=[[str(x) for x in pp] for pp in polys], verified_to=NV, failures=bad),
          open('runs/op_q%d_%d_%d.json' % (q, RHO, DG), 'w'))
