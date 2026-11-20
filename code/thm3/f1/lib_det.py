"""The 10x10 system in CDT's proof of Lemma 14.3.1, at numeric (a, b), exact over Q.

From (14.3.3): b2 = (c3/c4) b3 - b3', b1 = (c2/c4) b3 - b2', b0 = (c1/c4) b3 - b1',
               sum_{i=-4}^{1} a_i y^i = b0' - (c0/c4) b3,
with (14.3.5) b3 = (y-4)^2 R4 Q / R10, Q = q0 + ... + q9 y^9.  CDT (14.3.6): the right side is S(y) / (y^4 R10^4).
The 10 linear conditions 'R10 divides S' (coefficients of S mod R10) give the matrix; its determinant is CDT's
determinant (up to normalisation of the rows).
"""
from flint import fmpq, fmpq_poly, fmpq_mat
from lib_f1 import R10, R4, spec_y
from lib_L import solve_L, Y


class RF:
    def __init__(self, n, d=None):
        self.n = n
        self.d = d if d is not None else fmpq_poly([1])
        self._red()

    def _red(self):
        g = self.n.gcd(self.d)
        if g.degree() > 0:
            self.n, self.d = self.n // g, self.d // g
        lc = self.d[self.d.degree()]
        self.n, self.d = self.n / lc, self.d / lc

    def __add__(s, o):
        return RF(s.n * o.d + o.n * s.d, s.d * o.d)

    def __sub__(s, o):
        return RF(s.n * o.d - o.n * s.d, s.d * o.d)

    def __mul__(s, o):
        return RF(s.n * o.n, s.d * o.d)

    def der(s):
        return RF(s.n.derivative() * s.d - s.n * s.d.derivative(), s.d ** 2)


def system(a, b, cs=None):
    a, b = fmpq(a), fmpq(b)
    if cs is None:
        cs = solve_L(a, b)
    c0, c1, c2, c3, c4 = cs
    r10, r4 = spec_y(R10, a, b), spec_y(R4, a, b)
    q = [RF(c, c4) for c in (c0, c1, c2, c3)]
    Ss = []
    D = Y ** 4 * r10 ** 4
    for k in range(10):
        b3 = RF((Y - 4) ** 2 * r4 * Y ** k, r10)
        b2 = q[3] * b3 - b3.der()
        b1 = q[2] * b3 - b2.der()
        b0 = q[1] * b3 - b1.der()
        rhs = b0.der() - q[0] * b3
        # check shape S / (y^4 R10^4)
        S, rem = divmod(rhs.n * D, rhs.d)
        assert rem == 0, 'denominator not y^4 R10^4'
        Ss.append(S)
    rows = []
    for S in Ss:
        rr = S % r10
        rows.append([rr[i] for i in range(10)])
    M = fmpq_mat(rows).transpose()
    return M, Ss, r10
