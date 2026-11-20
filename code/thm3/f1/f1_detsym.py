"""K-3 closure, step B: the Lemma 14.3.1 matrix as explicit polynomials, and an a-priori degree bound for its determinant.

Input: runs/Lsym.txt (L~ = CDT's L, verified in f1_Lsym.py).  All arithmetic exact (fmpq_mpoly).
1. S_k(a, b, y) := y^4 R10^4 * rhs_k (k = 0..9), rhs_k from (14.3.3) with b3 = (y-4)^2 R4 y^k / R10 (14.3.5), computed
   symbolically; the shape S_k in Q[a, b, y], deg_y S_k <= 42 (CDT (14.3.6)) is checked (exact division).
2. In the variables (t, b, y), t := a - b: remainder r_k of S_k modulo R10 in y over Q(a, b).  lc_y(R10) = -32 t^4, so
   r_k = P_k / t^{o_k} with P_k in Q[t, b][y]; o_k is minimised (largest t-power common to the coefficients removed).
   Matrix (as in lib_det.system): M[i][k] = [y^i] r_k = M'[i][k] / t^{o_k}, M' polynomial.  det M = det M' / t^E, E = sum o_k.
3. A-PRIORI BOUND.  det M' is a sum of products of one entry from each column, so
      totdeg(det M') <= B := sum_k max_i totdeg(M'[i][k]).
   Total degree is invariant under (a, b) <-> (t, b).
Output: runs/detsym.pkl (M' entries as strings, o_k), runs/detsym.out.
"""
import time, pickle
from flint import fmpq, fmpq_mpoly_ctx
from lib_f1 import R10, R4

QC = fmpq_mpoly_ctx.get(('a', 'b', 'y'), 'lex')
TC = fmpq_mpoly_ctx.get(('t', 'b', 'y'), 'lex')
a_, b_, y_ = QC.gens()
t_, bt_, yt_ = TC.gens()


def zq(P):
    return QC.from_dict({m: fmpq(int(c)) for m, c in zip(P.monoms(), P.coeffs())})


def load_L():
    cs = []
    for line in open('runs/Lsym.txt'):
        s = line.strip()
        if s:
            cs.append(QC.from_dict(dict(eval_poly(s))))
    return cs


def eval_poly(s):
    P = eval(s.replace('^', '**'), {'a': a_, 'b': b_, 'y': y_})
    return zip(P.monoms(), P.coeffs())


class F:
    """N / (y^al R10^be), N in Q[a, b, y]."""
    R = None

    def __init__(self, N, al, be):
        self.N, self.al, self.be = N, al, be

    def lift(self, al, be):
        return self.N * y_ ** (al - self.al) * F.R ** (be - self.be)

    def __add__(s, o):
        al, be = max(s.al, o.al), max(s.be, o.be)
        return F(s.lift(al, be) + o.lift(al, be), al, be)

    def __neg__(s):
        return F(-s.N, s.al, s.be)

    def __sub__(s, o):
        return s + (-o)

    def der(s):
        R = F.R
        N = s.N.derivative('y') * y_ * R - s.al * s.N * R - s.be * s.N * y_ * R.derivative('y')
        return F(N, s.al + 1, s.be + 1)


def ycoeffs(P, ctx, nv):
    out = {}
    for m, c in zip(P.monoms(), P.coeffs()):
        out.setdefault(m[2], {})[(m[0], m[1], 0)] = c
    return {k: ctx.from_dict(v) for k, v in out.items()}


def tval(P):
    """t-adic valuation of a nonzero poly in TC (first variable t)."""
    return min(m[0] for m in P.monoms())


def shift_t(P, e):
    """P * t^e (e may be negative if divisible)."""
    return TC.from_dict({(m[0] + e, m[1], m[2]): c for m, c in zip(P.monoms(), P.coeffs())})


if __name__ == '__main__':
    log = open('runs/detsym.out', 'w')

    def say(*s):
        print(*s, flush=True)
        print(*s, file=log, flush=True)

    t0 = time.time()
    c0, c1, c2, c3, c4 = load_L()
    r4, r10 = zq(R4), zq(R10)
    F.R = r10
    assert c4 == (y_ - 4) ** 2 * y_ ** 3 * r4 * r10
    # q_i b3 = c_i y^k / (y^3 R10^2)   (c4 = (y-4)^2 y^3 R4 R10, b3 = (y-4)^2 R4 y^k / R10)
    S = []
    for k in range(10):
        b3 = F((y_ - 4) ** 2 * r4 * y_ ** k, 0, 1)
        qb = [F(c * y_ ** k, 3, 2) for c in (c0, c1, c2, c3)]
        b2 = qb[3] - b3.der()
        b1 = qb[2] - b2.der()
        b0 = qb[1] - b1.der()
        rhs = b0.der() - qb[0]
        # S = rhs * y^4 R10^4 must be a polynomial
        num = rhs.N
        if rhs.al > 4:
            q, rem = divmod(num, y_ ** (rhs.al - 4))
            assert rem == 0, 'y-denominator'
            num = q
        else:
            num = num * y_ ** (4 - rhs.al)
        if rhs.be > 4:
            q, rem = divmod(num, r10 ** (rhs.be - 4))
            assert rem == 0, 'R10-denominator'
            num = q
        else:
            num = num * r10 ** (4 - rhs.be)
        dy = num.degrees()[2]
        S.append(num)
        say('S_%d: deg_y %d, totdeg(a,b) of coefficients <= %d  (%.0fs)' % (
            k, dy, max(c.total_degree() for c in ycoeffs(num, QC, 3).values()), time.time() - t0))
        assert dy <= 42
    # change of variables a = t + b
    sub = [t_ + bt_, bt_, yt_]
    R10t = r10.compose(*sub, ctx=TC)
    lc = ycoeffs(R10t, TC, 3)[10]
    assert lc == -32 * t_ ** 4, lc
    Mp, offs = [], []
    for k in range(10):
        P = S[k].compose(*sub, ctx=TC)
        off = 0
        while True:
            co = ycoeffs(P, TC, 3)
            d = max(co)
            if d < 10:
                break
            ld = co[d]
            v = tval(ld)
            if v < 4:
                P = shift_t(P, 4 - v)
                off += 4 - v
                ld = shift_t(ld, 4 - v)
            qt = shift_t(ld, -4) / fmpq(-32)
            P = P - qt * yt_ ** (d - 10) * R10t
        co = ycoeffs(P, TC, 3)
        vmin = min(tval(c) for c in co.values())
        col = []
        for i in range(10):
            c = co.get(i, TC.from_dict({}))
            col.append(shift_t(c, -vmin) if not c.is_zero() else c)
        off -= vmin
        Mp.append(col)
        offs.append(off)
        say('column %d: o_k = %d, totdeg of entries %s  (%.0fs)' % (
            k, off, [c.total_degree() if not c.is_zero() else None for c in col], time.time() - t0))
    B = sum(max(c.total_degree() for c in col if not c.is_zero()) for col in Mp)
    E = sum(offs)
    degt = sum(max(c.degrees()[0] for c in col if not c.is_zero()) for col in Mp)
    degb = sum(max(c.degrees()[1] for c in col if not c.is_zero()) for col in Mp)
    say('E = sum o_k =', E)
    say('A-PRIORI BOUND: totdeg(det M\') <= B =', B, '; deg_t <= %d, deg_b (in t,b) <= %d' % (degt, degb))
    pickle.dump({'cols': [[{m: (int(c.p), int(c.q)) for m, c in zip(e.monoms(), e.coeffs())} for e in col] for col in Mp], 'offs': offs, 'B': B, 'E': E},
                open('runs/detsym.pkl', 'wb'))
    say('total %.0fs' % (time.time() - t0))
