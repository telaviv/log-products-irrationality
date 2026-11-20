"""i1_rec25.py -- I1' from I1: the order-25 operator of runs/op_q<q>_25_13.json is a left multiple of the order-4
operator of runs/rec_q<q>.json in Q(n)[S]:  mu(n) * Rec25 = Kt * L4  exactly (zero right remainder), with mu(n) a product
of shifts p_4(n+t), t = 0..21, times a nonzero integer, and Kt with polynomial coefficients.

Fraction-free right division: for b = 25..4: A <- p_4(n+b-4) A - A_b(n) S^{b-4} L4; integer contents removed.
Checks: (i) remainder (coefficients of S^0..S^3) is identically 0; (ii) p_4(n, q^2) has no root at an integer n >= 0
(so mu(n) != 0 for all n >= 0).  Then for any sequence X with L4 X = 0 for all n >= 0: mu(n) (Rec25 X)_n = 0, hence
Rec25 X = 0 for all n >= 0.
Output: runs/i1_rec25_q<q>.json.  Usage: python i1_rec25.py [q ...]   (default 26..32)
"""
import json, os, sys, time
from math import gcd
from flint import fmpz_poly, fmpz

HERE = os.path.dirname(os.path.abspath(__file__))


def shift(p, t):
    """p(n) -> p(n+t)"""
    return p(fmpz_poly([t, 1]))


def content(polys):
    g = fmpz(0)
    for p in polys:
        for c in p.coeffs():
            g = g.gcd(c) if g != 0 else abs(c)
    return g


def nonneg_int_roots(p):
    roots = []
    for fac, m in p.factor()[1]:
        if fac.degree() == 1:
            a, b = fac.coeffs()           # b n + a
            if a % b == 0 and -a // b >= 0:
                roots.append(int(-a // b))
    return roots


def run(q):
    t0 = time.time()
    d4 = json.load(open(os.path.join(HERE, 'runs', 'rec_q%d.json' % q)))
    d25 = json.load(open(os.path.join(HERE, 'runs', 'op_q%d_25_13.json' % q)))
    L = [fmpz_poly([int(x) for x in pp]) for pp in d4['polys']]
    A = [fmpz_poly([int(x) for x in pp]) for pp in d25['polys']]
    r4, r25 = len(L) - 1, len(A) - 1
    roots4 = nonneg_int_roots(L[r4])
    for b in range(r25, r4 - 1, -1):
        c = A[b]
        if c == 0:
            continue
        s = b - r4
        lead = shift(L[r4], s)
        A = [lead * x for x in A]
        for j in range(r4 + 1):
            A[j + s] = A[j + s] - c * shift(L[j], s)
        assert A[b] == 0
        g = content(A)
        if g > 1:
            A = [x / g for x in A]   # exact
    rem_zero = all(A[j] == 0 for j in range(r4))
    hi_zero = all(A[j] == 0 for j in range(r4, r25 + 1))
    res = dict(q=q, remainder_zero=rem_zero, all_zero=hi_zero and rem_zero,
               p4_nonneg_integer_roots=roots4, seconds=round(time.time() - t0, 1))
    json.dump(res, open(os.path.join(HERE, 'runs', 'i1_rec25_q%d.json' % q), 'w'))
    print(res, flush=True)
    return res


if __name__ == '__main__':
    qs = [int(x) for x in sys.argv[1:]] or list(range(26, 33))
    for q in qs:
        run(q)
