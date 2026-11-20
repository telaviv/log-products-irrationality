"""i1_contact.py -- the pieces of PROOF_I1.md section 3-4 that are checked by computer:

(a) termwise Rodrigues identity (symbolic in n, k): for 0 <= k <= 2n
      [w^k] D^n[ w^{n+1/2} D^{2n}[ w^{2n-1/2} (w-1)^{2n} ] ]  =  (2n)! n! / 16^n * (-1)^k B_k(n),
    B_k(n) = C(2n+k,2k) C(n+k,k) C(4n+2k,2n+k)   (D = d/dw).
(b) exact check, n <= 8, of  Q_n(w) = 16^n/(n!^2 (2n)!) D^n[w^n (w-1)^n D^n[w^{n+1/2} D^{2n}[w^{2n-1/2}(w-1)^{2n}]]]
    against rec.py's coefficients (w = z^2; half-integer powers handled with w = t^2).
(c) moment formula  [z^{-2m}] f3 = -2 int_0^1 x^{2m-1} log(1+1/x) dx  (m <= 30; exact: both sides rational)
    and an exact contact check: Q_n f_j - P_{j,n} = O(z^{-(2n+1)}) for n <= NC, j = 1, 2, 3 (sanity only; proof in PROOF_I1).
(d) the finitely many cases n = 0, 1, 2 of  L4 P_j = 0  as polynomial identities in z (the contact argument covers n >= 3);
    also n <= 8 for good measure.
Output: runs/i1_contact.json.   Usage: python i1_contact.py [NC]
"""
import json, os, sys, time
from fractions import Fraction as Fr
from math import comb, factorial
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
NC = int(sys.argv[1]) if len(sys.argv) > 1 else 24


def qcoefs(n):   # identical to rec.py
    cs = [0] * (6 * n + 1)
    for l in range(3 * n + 1):
        s = 0
        for k in range(max(0, l - n), min(2 * n, l) + 1):
            s += comb(n, l - k) * comb(2 * n + k, 2 * k) * comb(n + k, k) * comb(4 * n + 2 * k, 2 * n + k)
        cs[2 * l] = (-1) ** (n + l) * comb(n + l, l) * s
    return cs


def logcoef(K):   # identical to rec.py: coefficients of z^{-t} of f1 = log(1+1/z), f2 = log(1-1/z), f3 = f1 f2
    a = [Fr(0)] + [Fr((-1) ** (t + 1), t) for t in range(1, K + 1)]
    b = [Fr(0)] + [Fr(-1, t) for t in range(1, K + 1)]
    t3 = [Fr(0)] * (K + 1)
    for i in range(1, K + 1):
        for j in range(1, K + 1 - i):
            t3[i + j] += a[i] * b[j]
    return a, b, t3


if __name__ == '__main__':
    t0 = time.time()
    res = {}
    # (a)
    n, k = sp.symbols('n k', integer=True, nonnegative=True)
    lhs = sp.binomial(2 * n, k) * (-1) ** k * sp.gamma(2 * n + k + sp.Rational(1, 2)) / sp.gamma(k + sp.Rational(1, 2)) \
        * sp.factorial(n + k) / sp.factorial(k)
    Bk = sp.binomial(2 * n + k, 2 * k) * sp.binomial(n + k, k) * sp.binomial(4 * n + 2 * k, 2 * n + k)
    ratio = sp.simplify(sp.gammasimp(sp.combsimp(lhs / ((-1) ** k * Bk * sp.factorial(2 * n) * sp.factorial(n) / 16 ** n))))
    print('(a) termwise ratio (should be 1):', ratio)
    res['a_termwise_ratio'] = str(ratio)
    # numeric confirmation of (a) for n <= 12 (exact rationals; Gamma(j+1/2)/Gamma(i+1/2) = prod (i+1/2)...(j-1/2))
    oka = True
    for N in range(0, 13):
        for K in range(0, 2 * N + 1):
            g = Fr(1)
            for s in range(K, 2 * N + K):
                g *= Fr(2 * s + 1, 2)
            v = comb(2 * N, K) * (-1) ** K * g * Fr(factorial(N + K), factorial(K))
            B = comb(2 * N + K, 2 * K) * comb(N + K, K) * comb(4 * N + 2 * K, 2 * N + K)
            oka &= (v == Fr((-1) ** K * B * factorial(2 * N) * factorial(N), 16 ** N))
    print('(a) exact check n <= 12:', oka)
    res['a_exact_n_le_12'] = oka
    # (b)
    t = sp.symbols('t')
    def Dw(f, m):   # (d/dw)^m with w = t^2
        for _ in range(m):
            f = sp.expand(sp.diff(f, t) / (2 * t))
        return f
    okb = True
    for N in range(0, 9):
        inner = Dw(t ** (4 * N - 1) * (t ** 2 - 1) ** (2 * N), 2 * N)
        V = Dw(sp.expand(t ** (2 * N + 1) * inner), N) * sp.Rational(16 ** N, factorial(2 * N) * factorial(N))
        Q = Dw(sp.expand(t ** (2 * N) * (t ** 2 - 1) ** N * V), N) / factorial(N)
        cs = qcoefs(N)
        target = sum(cs[i] * t ** i for i in range(len(cs)))   # z = t
        okb &= sp.expand(Q - target) == 0
    print('(b) Rodrigues form = rec.py Q_n, n <= 8:', okb)
    res['b_rodrigues_n_le_8'] = okb
    # (c) moments of f3 and exact contact
    K = 6 * (NC + 8) + 10
    a1, a2, a3 = logcoef(K)
    okm = True
    for m in range(1, 31):
        A2m = sum(Fr((-1) ** (j + 1), j) for j in range(1, 2 * m + 1))
        integral = A2m / (2 * m) + Fr(1, 4 * m * m)    # int_0^1 x^{2m-1} log(1+1/x) dx (log 2 terms cancel)
        okm &= (a3[2 * m] == -2 * integral)
    print('(c) f3 moment formula m <= 30:', okm)
    res['c_f3_moments_m_le_30'] = okm
    okc = True
    for N in range(0, NC + 1):
        cs = qcoefs(N)
        for phi in (a1, a2, a3):
            for s in range(1, 2 * N + 1):     # coefficient of z^{-s} of Q_n f_j
                if sum(cs[i] * phi[i + s] for i in range(len(cs)) if cs[i]) != 0:
                    okc = False
    print('(c) exact contact O(z^-(2n+1)) for n <= %d, j = 1,2,3:' % NC, okc)
    res['c_contact_exact_n_le'] = NC if okc else -1
    # (d) small-n polynomial identities L4 P_j = 0 in Q[z]
    L4 = json.load(open(os.path.join(HERE, 'runs', 'i1_L4_symbolic.json')))
    p = [{tuple(int(x) for x in key.split(',')): sp.Rational(v) for key, v in pi.items()} for pi in L4['p']]
    def P_poly(N, phi):
        cs = qcoefs(N)
        d = 6 * N
        out = [Fr(0)] * (d + 1)       # out[r] = [z^r] polypart(Q_N f) = sum_{i>r} q_i phi_{i-r}
        for r in range(d):
            out[r] = sum(Fr(cs[i]) * phi[i - r] for i in range(r + 1, d + 1) if cs[i])
        return out
    def Q_poly(N):
        return [Fr(x) for x in qcoefs(N)]
    okd = {}
    for N in range(0, 9):
        good = True
        for name, phi in (('Q', None), ('P1', a1), ('P2', a2), ('P3', a3)):
            tot = {}
            for i in range(5):
                X = Q_poly(N + i) if phi is None else P_poly(N + i, phi)
                # p_i(N, z^2) as dict power-of-z -> coeff
                pz = {}
                for (en, ew), cf in p[i].items():
                    pz[2 * ew] = pz.get(2 * ew, Fr(0)) + Fr(int(cf.p), int(cf.q)) * N ** en
                for e1, c1 in pz.items():
                    if c1 == 0:
                        continue
                    for e2, c2 in enumerate(X):
                        if c2:
                            tot[e1 + e2] = tot.get(e1 + e2, Fr(0)) + c1 * c2
            if any(v != 0 for v in tot.values()):
                good = False
        okd[N] = good
    print('(d) L4 X = 0 in Q[z] for X = Q, P1, P2, P3:', okd)
    res['d_small_n_polynomial_identities'] = okd
    res['seconds'] = round(time.time() - t0, 1)
    json.dump(res, open(os.path.join(HERE, 'runs', 'i1_contact.json'), 'w'), indent=1)
    print('done (%.1fs)' % (time.time() - t0))
