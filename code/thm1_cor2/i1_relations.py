"""i1_relations.py -- find (by exact linear algebra over Z) the two small relations of Lysov's coefficients
c_{n,l} (Q_n(z) = sum_l c_{n,l} z^{2l}, the pilot's rec.py qcoefs) that drive the proof of I1:

  (R1)  a_0(n,l) c_{n,l} + a_1(n,l) c_{n,l+1} + a_2(n,l) c_{n,l+2} = 0                    (all n >= 0, l in Z)
  (R2)  alpha(n) c_{n+1,l} = sum_{a=0}^{3} sum_{i=0}^{3} beta_{a,i}(n) (l-i)^a c_{n,l-i}    (all n >= 0, l in Z)

In generating-function form (w = z^2, theta = w d/dw):
  (R1) <=> D_n Q_n = 0,  D_n = sum_m w^{2-m} a_m(n, theta - m)          (order-4 ODE in w)
  (R2) <=> alpha(n) Q_{n+1} = sum_a beta_a(n, w) theta^a Q_n,  beta_a(n,w) = sum_i beta_{a,i}(n) w^i.

The relations are only *found* here (guessing); they are *proved* by the certificates of i1_certs.py.
Output: runs/i1_relations.json (integer coefficient lists, sympy-readable strings).
Usage: python i1_relations.py
"""
import json, os
from math import comb, gcd
from flint import fmpz_mat
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
n, l, w = sp.symbols('n l w')


def F(N, L, K):
    """Lysov's summand: c_{N,L} = sum_K F(N,L,K) (natural boundaries)."""
    if L < 0 or K < 0 or K > 2 * N or L - K < 0 or L - K > N:
        return 0
    return (-1) ** (N + L) * comb(N + L, L) * comb(N, L - K) * comb(2 * N + K, 2 * K) * comb(N + K, K) * comb(4 * N + 2 * K, 2 * N + K)


def c(N, L):
    return sum(F(N, L, K) for K in range(0, 2 * N + 1))


def intnull(rows):
    X, nul = fmpz_mat(rows).nullspace()
    assert nul == 1, nul
    v = [int(X[i, 0]) for i in range(X.nrows())]
    g = 0
    for x in v:
        g = gcd(g, x)
    return [x // g for x in v]


if __name__ == '__main__':
    NM = 22
    C = {(N, L): c(N, L) for N in range(NM + 2) for L in range(-6, 3 * N + 12)}
    # (R1): order 2 in l, degree <= 4 in n and in l
    cols = [(m, i, j) for m in range(3) for i in range(5) for j in range(5)]
    rows = [[N ** i * L ** j * C[(N, L + m)] for (m, i, j) in cols] for N in range(NM) for L in range(-4, 3 * N + 4)]
    v = intnull(rows)
    a = [sp.expand(sum(v[cols.index((m, i, j))] * n ** i * l ** j for i in range(5) for j in range(5))) for m in range(3)]
    # (R2): alpha(n) c_{n+1,l} = sum beta_{a,i}(n) (l-i)^a c_{n,l-i}, n-degree <= 7, w-degree 0..3, theta-order 3
    A, WH, DN = 3, 3, 7
    cols = [('a', i, e) for i in range(WH + 1) for e in range(DN + 1)] + \
           [('b', aa, i, e) for aa in range(A + 1) for i in range(WH + 1) for e in range(DN + 1)]
    rows = []
    for N in range(NM - 1):
        for L in range(-2, 3 * N + 3 + WH + 2):
            row = []
            for col in cols:
                if col[0] == 'a':
                    _, i, e = col
                    row.append(N ** e * C.get((N + 1, L - i), 0))
                else:
                    _, aa, i, e = col
                    row.append(-N ** e * (L - i) ** aa * C.get((N, L - i), 0))
            rows.append(row)
    v = intnull(rows)
    alpha_w = sp.expand(sum(v[cols.index(('a', i, e))] * n ** e * w ** i for i in range(WH + 1) for e in range(DN + 1)))
    assert sp.degree(alpha_w, w) == 0, 'alpha depends on w'
    beta = [sp.expand(sum(v[cols.index(('b', aa, i, e))] * n ** e * w ** i for i in range(WH + 1) for e in range(DN + 1)))
            for aa in range(A + 1)]
    out = dict(a=[str(x) for x in a], alpha=str(alpha_w), beta=[str(x) for x in beta],
               note='R1: sum_m a_m(n,l) c_{n,l+m} = 0;  R2: alpha(n) c_{n+1,l} = sum_{a,i} [w^i]beta_a(n,w) (l-i)^a c_{n,l-i}')
    for x in a:
        print('a:', sp.factor(x))
    print('alpha:', sp.factor(alpha_w))
    for x in beta:
        print('beta:', sp.factor(x))
    # sanity: both relations hold on all computed data (incl. l < 0 and l > 3n)
    for N in range(NM - 1):
        for L in range(-6, 3 * N + 6):
            s1 = sum(int(a[m].subs({n: N, l: L})) * C.get((N, L + m), 0) for m in range(3))
            assert s1 == 0
            s2 = int(alpha_w.subs(n, N)) * C.get((N + 1, L), 0) - sum(
                int(sp.Poly(beta[aa], w).coeff_monomial(w ** i).subs(n, N)) * (L - i) ** aa * C.get((N, L - i), 0)
                for aa in range(A + 1) for i in range(WH + 1))
            assert s2 == 0
    print('data check n <=', NM - 2, ': OK')
    json.dump(out, open(os.path.join(HERE, 'runs', 'i1_relations.json'), 'w'), indent=1)
