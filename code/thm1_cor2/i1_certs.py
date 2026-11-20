"""i1_certs.py -- creative-telescoping certificates proving (R1), (R2) of runs/i1_relations.json for ALL n >= 0, l in Z.

c_{n,l} = sum_k F(n,l,k),  F = (-1)^{n+l} C(n+l,l) C(n,l-k) C(2n+k,2k) C(n+k,k) C(4n+2k,2n+k)   (0 outside the natural range).
For n, l', k >= 0 (and also for l' < 0 when n + l' >= 0) F has the factorial form
   F(n,l,k) = (-1)^{n+l} (n+l)! (n+k)! (4n+2k)! / ( l! (l-k)! (n-l+k)! (2k)! (2n-k)! k! n! (2n+k)! ),   1/(m!) := 0 for m < 0.

Certificates G_r (r = 1, 2) are written in the same factorial form (pole factors absorbed into shifted factorials,
e.g. 1/((l-k)!(l-k+1)(l-k+2)) = 1/(l-k+2)!), and the telescoping identity
      sum_s coeff_s(n,l) F(n+dn_s, l+dl_s, k) = G_r(n,l,k+1) - G_r(n,l,k)
is checked as a POLYNOMIAL identity after rewriting every term over one common factorial base
(1/(L+c)! = (L+c+1)...(L+b)/(L+b)!  for c <= b, valid at every integer point;  (L+c)! = (L+c)...(L+b+1) (L+b)!  for
c >= b, valid when L+b >= 0).  So the identity holds at every integer point of the stated region, without any
'generic point' caveat.  Boundary handling (k = -1, l < 0, small n + l) is done explicitly below.

Output: runs/i1_certs.json (certificates as strings, check results).  Usage: python i1_certs.py
"""
import json, os, time
from math import comb
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))
n, l, k, w = sp.symbols('n l k w')
rel = json.load(open(os.path.join(HERE, 'runs', 'i1_relations.json')))
A_ = [sp.sympify(x) for x in rel['a']]
ALPHA = sp.sympify(rel['alpha'])
BETA = [sp.Poly(sp.sympify(x), w) for x in rel['beta']]

# linear forms (cn, cl, ck, c0)
NUM_F = [(1, 1, 0, 0), (1, 0, 1, 0), (4, 0, 2, 0)]
DEN_F = [(0, 1, 0, 0), (0, 1, -1, 0), (1, -1, 1, 0), (0, 0, 2, 0), (2, 0, -1, 0), (0, 0, 1, 0), (1, 0, 0, 0), (2, 0, 1, 0)]


def lf(f):
    return f[0] * n + f[1] * l + f[2] * k + f[3]


def shift_form(f, dn, dl, dk):
    return (f[0], f[1], f[2], f[3] + f[0] * dn + f[1] * dl + f[2] * dk)


class Term:
    """sign (-1)^{n+l+par} * poly(n,l,k) * prod NUM! / prod DEN!"""
    def __init__(self, poly, num, den, par=0):
        self.poly, self.num, self.den, self.par = sp.sympify(poly), list(num), list(den), par

    def shifted(self, dn=0, dl=0, dk=0):
        return Term(self.poly.subs({n: n + dn, l: l + dl, k: k + dk}, simultaneous=True),
                    [shift_form(f, dn, dl, dk) for f in self.num], [shift_form(f, dn, dl, dk) for f in self.den],
                    self.par + dn + dl)

    def scaled(self, c):
        return Term(self.poly * c, self.num, self.den, self.par)


F_TERM = Term(1, NUM_F, DEN_F)


def common_base_check(terms):
    """returns (sum of reduced polynomials == 0, base num, base den)"""
    def groups(lst_of_lists):
        keys = None
        for lst in lst_of_lists:
            ks = sorted(f[:3] for f in lst)
            assert keys is None or ks == keys, 'shape mismatch'
            keys = ks
        return keys
    gn = groups([t.num for t in terms]); gd = groups([t.den for t in terms])
    bnum = {key: min(f[3] for t in terms for f in t.num if f[:3] == key) for key in gn}
    bden = {key: max(f[3] for t in terms for f in t.den if f[:3] == key) for key in gd}
    total = 0
    for t in terms:
        mu = sp.Integer(1)
        for f in t.num:
            b = bnum[f[:3]]
            L = f[0] * n + f[1] * l + f[2] * k
            for c in range(b + 1, f[3] + 1):
                mu *= (L + c)
        for f in t.den:
            b = bden[f[:3]]
            L = f[0] * n + f[1] * l + f[2] * k
            for c in range(f[3] + 1, b + 1):
                mu *= (L + c)
        total += (-1) ** (t.par % 2) * t.poly * mu
    return sp.expand(total) == 0, bnum, bden


def num_min_forms(bnum):
    return [sp.expand(key[0] * n + key[1] * l + key[2] * k + c) for key, c in bnum.items()]


# ---------- exact binomial F, c for direct checks
def Fb(N, L, K):
    if L < 0 or K < 0 or K > 2 * N or L - K < 0 or L - K > N:
        return 0
    return (-1) ** (N + L) * comb(N + L, L) * comb(N, L - K) * comb(2 * N + K, 2 * K) * comb(N + K, K) * comb(4 * N + 2 * K, 2 * N + K)


def cb(N, L):
    return sum(Fb(N, L, K) for K in range(0, 2 * N + 1))


def eval_term(t, N, L, K):
    """numeric value of a factorial-form term (requires numerator args >= 0)"""
    from math import factorial
    v = sp.Rational(int(t.poly.subs({n: N, l: L, k: K}))) if t.poly.free_symbols else sp.Rational(t.poly)
    if v == 0:
        return 0
    for f in t.num:
        a = f[0] * N + f[1] * L + f[2] * K + f[3]
        assert a >= 0
        v *= factorial(a)
    for f in t.den:
        a = f[0] * N + f[1] * L + f[2] * K + f[3]
        if a < 0:
            return 0
        v /= factorial(a)
    return (-1) ** ((N + L + t.par) % 2) * v


def r2_certificate():
    """solve R(k+1)u - R v = rho v for R = k^2(2k-1) Y(k) / Den, Y cubic in k (coefficients in Q(n,l))."""
    def fall(x, m): return sp.Mul(sp.Integer(1), *[x - j for j in range(m)])
    def rise(x, m): return sp.Mul(sp.Integer(1), *[x + j for j in range(m)])
    u = (n + k + 1) * (4 * n + 2 * k + 1) * (l - k) * (2 * n - k)
    v = (n - l + k + 1) * (2 * k + 1) * (k + 1) ** 2
    r_n = -(n + l + 1) * (n + k + 1) * rise(4 * n + 2 * k + 1, 4) / (
        (n - l + k + 1) * (2 * n - k + 1) * (2 * n - k + 2) * (n + 1) * (2 * n + k + 1) * (2 * n + k + 2))
    def r_l(i): return sp.Integer(-1) ** i * fall(l, i) * fall(l - k, i) / (fall(n + l, i) * rise(n - l + k + 1, i))
    rho = ALPHA * r_n
    for a in range(4):
        for i in range(4):
            rho -= BETA[a].coeff_monomial(w ** i) * (l - i) ** a * r_l(i)
    Den = (2 * n + 1 - k) * (2 * n + 2 - k) * (n - l + k + 1) * (n - l + k + 2)
    y = sp.symbols('y0:4')
    R = k ** 2 * (2 * k - 1) * sum(y[j] * k ** j for j in range(4)) / Den
    expr = sp.together(R.subs(k, k + 1) * u - R * v - rho * v)
    eqs = sp.Poly(sp.expand(sp.numer(expr)), k).coeffs()
    sol = sp.solve(eqs, y, dict=True)
    assert len(sol) == 1
    return sp.factor(sum(sol[0][y[j]] * k ** j for j in range(4)))


if __name__ == '__main__':
    t0 = time.time()
    out = {}
    # ================= (R1)
    # certificate R1 = (n+1)(l+n+1)(l+n+2)/(l+1) * k^2 (2k-1)(k-l+n) / ((l+1-k)(l+2-k)),  G1 = R1 * F:
    P1 = (n + 1) * (l + n + 1) * (l + n + 2) * k ** 2 * (2 * k - 1) * (k - l + n)
    DEN_G1 = [(0, 1, 0, 1), (0, 1, -1, 2)] + DEN_F[2:]          # l! -> (l+1)!, (l-k)! -> (l-k+2)!
    G1 = Term(P1, NUM_F, DEN_G1)
    terms = [F_TERM.shifted(dl=m).scaled(A_[m]) for m in range(3)] + [G1.shifted(dk=1).scaled(-1), G1]
    ok1, bnum1, bden1 = common_base_check(terms)
    print('R1 telescoping identity (polynomial, common factorial base):', ok1, ' numerator base forms:', num_min_forms(bnum1))
    # region: base numerator args >= 0  <=>  n+l >= 0, n+k >= 0, 4n+2k >= 0: holds for n >= 0, k >= 0, l >= -n.
    # remaining cases with l >= -2 and n + l < 0: (n,l) = (0,-1), (0,-2), (1,-2): direct.  l <= -3: all c = 0.
    direct1 = all(sum(int(A_[m].subs({n: N, l: L})) * cb(N, L + m) for m in range(3)) == 0
                  for (N, L) in [(0, -1), (0, -2), (1, -2)])
    # G1(n,l,0) = 0 (factor k^2) and G1 finitely supported in k (1/(2n-k)!): telescoping sum over k in Z is 0.
    g1_k0 = sp.expand(P1.subs(k, 0)) == 0
    print('R1 direct boundary cases:', direct1, ' G1(n,l,0)=0:', g1_k0)
    # numeric sanity with the binomial definition (independent of the factorial bookkeeping)
    bad = 0
    for N in range(0, 7):
        for L in range(-N, 3 * N + 4):
            for K in range(-1, 2 * N + 3):
                lhs = sum(int(A_[m].subs({n: N, l: L})) * Fb(N, L + m, K) for m in range(3))
                g = lambda KK: eval_term(G1, N, L, KK) if KK >= 0 else 0
                if lhs != g(K + 1) - g(K):
                    bad += 1
    print('R1 numeric sanity (n<=6, all l>=-n, k): failures', bad)
    out['R1'] = dict(certificate='G1 = (-1)^(n+l) (n+1)(l+n+1)(l+n+2) k^2 (2k-1)(k-l+n) (n+l)!(n+k)!(4n+2k)! / '
                                 '((l+1)! (l-k+2)! (n-l+k)! (2k)! (2n-k)! k! n! (2n+k)!)',
                     rational_form='R1 = (n+1)(l+n+1)(l+n+2) k^2 (2k-1)(k-l+n) / ((l+1)(l+1-k)(l+2-k)),  G1 = R1*F',
                     identity_ok=ok1, direct_boundary_ok=direct1, G_at_k0_zero=g1_k0, numeric_failures=bad,
                     region='n>=0, k>=0, l>=-n (factorial bases); l<=-3 trivial; (0,-1),(0,-2),(1,-2) direct')
    # ================= (R2)
    Y = r2_certificate()
    print('R2 certificate Y =', Y, ' (%.1fs)' % (time.time() - t0))
    Ynum, Yden = sp.fraction(sp.factor(Y))
    assert sp.expand(Yden - (l + n)) == 0 or sp.expand(Yden + (l + n)) == 0, Yden
    sgn = 1 if sp.expand(Yden - (l + n)) == 0 else -1
    # G2 = k^2(2k-1) Y / ((2n+1-k)(2n+2-k)(n-l+k+1)(n-l+k+2)) * F;   1/(l+n) absorbed: (n+l)!/(n+l) = (n+l-1)!
    P2 = sgn * k ** 2 * (2 * k - 1) * Ynum
    NUM_G2 = [(1, 1, 0, -1)] + NUM_F[1:]
    DEN_G2 = [DEN_F[0], DEN_F[1], (1, -1, 1, 2), DEN_F[3], (2, 0, -1, 2), DEN_F[5], DEN_F[6], DEN_F[7]]
    G2 = Term(P2, NUM_G2, DEN_G2)
    terms = [F_TERM.shifted(dn=1).scaled(ALPHA)]
    for a in range(4):
        for i in range(4):
            cf = BETA[a].coeff_monomial(w ** i)
            if cf != 0:
                terms.append(F_TERM.shifted(dl=-i).scaled(-cf * (l - i) ** a))
    terms += [G2.shifted(dk=1).scaled(-1), G2]
    ok2, bnum2, bden2 = common_base_check(terms)
    print('R2 telescoping identity (polynomial, common factorial base):', ok2, ' numerator base forms:', num_min_forms(bnum2))
    # region: n+l-3 >= 0, n+k >= 0, 4n+2k >= 0  (n,l,k >= 0, n+l >= 3).  l < 0: every c in R2 vanishes.
    # remaining: l >= 0, n + l <= 2: direct.
    def r2_direct(N, L):
        s = int(ALPHA.subs(n, N)) * cb(N + 1, L)
        for a in range(4):
            for i in range(4):
                s -= int(BETA[a].coeff_monomial(w ** i).subs(n, N)) * (L - i) ** a * cb(N, L - i)
        return s == 0
    direct2 = all(r2_direct(N, L) for N in range(3) for L in range(3) if N + L <= 2)
    g2_k0 = sp.expand(P2.subs(k, 0)) == 0
    print('R2 direct cases (n+l<=2):', direct2, ' G2(n,l,0)=0:', g2_k0)
    bad = 0
    for N in range(0, 7):
        for L in range(0, 3 * N + 7):
            if N + L < 3:
                continue
            for K in range(-1, 2 * N + 4):
                lhs = int(ALPHA.subs(n, N)) * Fb(N + 1, L, K)
                for a in range(4):
                    for i in range(4):
                        lhs -= int(BETA[a].coeff_monomial(w ** i).subs(n, N)) * (L - i) ** a * Fb(N, L - i, K)
                g = lambda KK: eval_term(G2, N, L, KK) if KK >= 0 else 0
                if lhs != g(K + 1) - g(K):
                    bad += 1
    print('R2 numeric sanity (n<=6, n+l>=3, all k): failures', bad)
    out['R2'] = dict(Y=str(Y), certificate='G2 = (-1)^(n+l) * k^2 (2k-1) * [(l+n) Y] * (n+l-1)!(n+k)!(4n+2k)! / '
                     '(l! (l-k)! (n-l+k+2)! (2k)! (2n-k+2)! k! n! (2n+k)!)',
                     rational_form='R2 = k^2 (2k-1) Y / ((2n+1-k)(2n+2-k)(n-l+k+1)(n-l+k+2)),  G2 = R2*F',
                     identity_ok=ok2, direct_cases_ok=direct2, G_at_k0_zero=g2_k0, numeric_failures=bad,
                     region='n,l,k>=0 with n+l>=3 (factorial bases); l<0 trivial; n+l<=2 direct')
    json.dump(out, open(os.path.join(HERE, 'runs', 'i1_certs.json'), 'w'), indent=1)
    print('all:', ok1 and direct1 and g1_k0 and ok2 and direct2 and g2_k0, ' (%.1fs)' % (time.time() - t0))
