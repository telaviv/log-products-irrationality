"""Phi6(a, b) != 0 near the diagonal (exact; fmpq).  Also: consistency of runs/phi6.txt with the univariate
factorisations runs/det_b*.json, a congruence search, and a brute-force search for integer zeros.

Claims: (A) Phi6(a, b) = 16 mod 32 for all odd a, b (section 4b), so Phi6 has no zero with a, b odd.
(B) for all integers a != b with |a/b - 1| < delta (delta = 1/1000, CDT Lemma 14.2.6), Phi6(a, b) != 0.
Proof (checked here): put b = a(1 + s), |s| <= delta' := delta/(1 - delta) (same sign, |b - a| <= delta|b|), so
|a| >= 1/delta' (b - a is a nonzero integer).  Phi6(a, a(1+s)) = sum_k a^k P_k(s), P_k(s) = Phi6^{(k)}(1, 1+s) (homogeneous
parts, k = 8, 6, 4, 2, 0).  P_8(0) = Phi6^{(8)}(1,1) = -2541.  Bound |P_8(s) - P_8(0)| <= sum_{j>=1} |p_{8,j}| delta'^j and
|P_k(s)| <= sum_j |p_{k,j}| delta'^j for k < 8; then |Phi6| >= |a|^8 (2541 - err_8 - sum_{k<8} B_k |a|^(k-8)) > 0 for
|a| >= 1/delta'.
"""
import json, glob
from flint import fmpq, fmpz, fmpz_poly, fmpq_poly
from lib_f1 import CTX, A_, B_, Y_

P = eval(open('runs/phi6.txt').read().replace('^', '**'), {'a': A_, 'b': B_})
terms = [(ea, eb, int(c)) for (ea, eb, ey), c in zip(P.monoms(), P.coeffs())]


def phi6(a, b):
    return sum(c * a ** i * b ** j for i, j, c in terms)


# 1. consistency with univariate determinant factors
print('1. Phi6(a, b0) vs the degree-6 factor of det(a, b0):')
for fn in sorted(glob.glob('runs/det_b*.json')):
    d = json.load(open(fn))
    b0 = d['b0']
    f6 = [f for f in d['factors'] if f['deg'] == 6]
    u = fmpz_poly([sum(c * b0 ** j for i, j, c in terms if i == k) for k in range(7)])
    prim = u // u.content() if u.content() != 0 else u
    if prim.leading_coefficient() < 0:
        prim = -prim
    print('   b0 = %3d: degree-6 factors %d, equal to primitive part of Phi6(a, b0): %s'
          % (b0, len(f6), any(fmpz_poly([int(c) for c in f['coeffs']]) == prim for f in f6)))

# 2. homogeneous parts along the diagonal
s = fmpq_poly([0, 1])
Pk = {}
for i, j, c in terms:
    k = i + j
    Pk[k] = Pk.get(k, fmpq_poly([0])) + c * (1 + s) ** j
for k in sorted(Pk, reverse=True):
    print('2. P_%d(s) = %s' % (k, Pk[k]))
print('   Phi6(a, a) =', ' + '.join('(%s) a^%d' % (Pk[k](0), k) for k in sorted(Pk, reverse=True)))


def lower_bound(delta):
    dp = delta / (1 - delta)
    amin = 1 / dp
    err8 = sum(abs(Pk[8][j]) * dp ** j for j in range(1, Pk[8].degree() + 1))
    tot = abs(Pk[8](0)) - err8
    for k in Pk:
        if k < 8:
            Bk = sum(abs(Pk[k][j]) * dp ** j for j in range(Pk[k].degree() + 1))
            tot -= Bk * amin ** (k - 8)
    return tot, amin


for delta in (fmpq(3, 10 ** 6), fmpq(1, 10 ** 3), fmpq(1, 100), fmpq(1, 20)):
    lb, amin = lower_bound(delta)
    print('3. delta = %s: |Phi6(a,b)| / |a|^8 >= %s (= %.4f) for |a| >= %.1f -> %s'
          % (delta, lb, float(lb.p) / float(lb.q), float(amin.p) / float(amin.q), 'NONZERO' if lb > 0 else 'inconclusive'))

# 4. congruence obstruction search (odd a, b)
print('4. moduli m with Phi6(a,b) != 0 mod m for all odd a, b:',
      [m for m in range(2, 200) if all(phi6(x, y) % m != 0 for x in range(m) for y in range(m)
                                       if (m % 2 == 1 or (x % 2 == 1 and y % 2 == 1)))][:5])

# 5. integer zeros in a box
Z = [(x, y) for x in range(-300, 301) for y in range(-300, 301) if phi6(x, y) == 0]
print('5. integer zeros with |a|, |b| <= 300:', Z)

# 4b. the 2-adic argument used in RESULT.md section 5: Phi6(a, b) = 16 mod 32 for all odd a, b
#     (integer coefficients, so the residue depends only on a, b mod 32; exhaustive over the 16 x 16 odd classes)
res32 = sorted(set(phi6(x, y) % 32 for x in range(1, 32, 2) for y in range(1, 32, 2)))
print('4b. {Phi6(a,b) mod 32 : a, b odd} =', res32, '-> Phi6(a, b) != 0 for all odd a, b:', res32 == [16])
