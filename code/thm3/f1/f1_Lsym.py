"""K-3 closure, step A: CDT's operator L = sum c_i(y) D^i as an explicit element of Q[a, b, y][D], and its verification.

1. Interpolate every y-coefficient of c_0..c_3 (solve_L, exact) on a 13 x 13 grid in (a, b) (bidegree <= 12; the
   univariate probe gives a-degree 12).  c_4 = (y-4)^2 y^3 R4 R10 is CDT's normalisation (14.2.1).
2. RIGOROUS CHECK (no interpolation trusted): with G_A = sum g_m y^m, g_m in Q[a, b] computed symbolically, the
   candidate L~ satisfies [y^k] L~(G_A) = 0 identically in Q[a, b] for k = 0..K0-1 (K0 = 63).  At (a, b) = (7, 5) the
   first K0 rows of the ansatz system already have rank 63 (= number of unknowns), so the generic rank of these rows
   is 63, the ansatz solution over Q(a, b) is unique, and any operator of CDT's printed shape annihilating G_A
   (CDT's L, Prop. 14.1.11 / (14.2.1)) equals L~ in Q(a, b)[y][D].  At every (a0, b0) where solve_L has rank 63 its
   output therefore equals L~(a0, b0) (spot-checked below at fresh points).
3. Also checks c3/c4 = d/dy log((y-4)^3 y^5 R4^3 / R10) symbolically (CDT's printed formula).
Output: runs/Lsym.txt (c0..c4 as polynomials in a, b, y), runs/Lsym.out.
"""
import time, sys
from flint import fmpq, fmpq_poly, fmpq_mpoly_ctx, fmpz_mpoly_ctx
from lib_L import solve_L
from lib_f1 import R10, R4, legendre

QC = fmpq_mpoly_ctx.get(('a', 'b', 'y'), 'lex')
a_, b_, y_ = QC.gens()
K0 = 63
DEG = 12


def zq(P):
    """fmpz_mpoly (a, b, y) -> fmpq_mpoly in QC."""
    return QC.from_dict({m: fmpq(int(c)) for m, c in zip(P.monoms(), P.coeffs())})


def newton(xs, ys):
    n = len(xs)
    c = list(ys)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            c[i] = (c[i] - c[i - 1]) / (xs[i] - xs[i - j])
    P = fmpq_poly([c[-1]])
    X = fmpq_poly([0, 1])
    for i in range(n - 2, -1, -1):
        P = P * (X - xs[i]) + c[i]
    return P


def interpolate_L():
    A = [fmpq(v) for v in (7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31)]
    B = [fmpq(v) for v in (2, 3, 4, 5, -3, -5, 6, 8, 10, 12, 14, 16, 18)]
    vals = {(x, z): solve_L(x, z) for x in A for z in B}
    cs = []
    for i in range(4):
        dy = max(vals[(x, z)][i].degree() for x in A for z in B)
        P = QC.from_dict({})
        for j in range(dy + 1):
            # for each b in B: univariate in a; then interpolate coefficients in b
            pa = {z: newton(A, [vals[(x, z)][i][j] for x in A]) for z in B}
            for ea in range(DEG + 1):
                pb = newton(B, [pa[z][ea] for z in B])
                for eb in range(pb.degree() + 1):
                    if pb[eb] != 0:
                        P += QC.from_dict({(ea, eb, j): pb[eb]})
        cs.append(P)
    return cs


def spec(P, a, b):
    out = {}
    for (ea, eb, ey), c in zip(P.monoms(), P.coeffs()):
        out[ey] = out.get(ey, fmpq(0)) + c * fmpq(a) ** ea * fmpq(b) ** eb
    d = max(out) if out else 0
    return fmpq_poly([out.get(i, fmpq(0)) for i in range(d + 1)])


def g_series(M):
    """g_0..g_{M-1} in Q[a, b] (fmpq_mpoly in QC, no y).  [y^m] p_n = 0 for n > 2m."""
    N = 2 * M
    # Legendre polynomials as mpolys in a and in b
    Pa = [QC.from_dict({(0, 0, 0): 1}), a_]
    Pb = [QC.from_dict({(0, 0, 0): 1}), b_]
    for n in range(1, N):
        Pa.append((fmpq(2 * n + 1) * a_ * Pa[n] - fmpq(n) * Pa[n - 1]) / fmpq(n + 1))
        Pb.append((fmpq(2 * n + 1) * b_ * Pb[n] - fmpq(n) * Pb[n - 1]) / fmpq(n + 1))
    # p_n(y) truncated mod y^M, as integer lists
    pm2 = [2] + [0] * (M - 1)
    pm1 = [0, 1] + [0] * (M - 2)
    g = [QC.from_dict({}) for _ in range(M)]
    for n in range(N + 1):
        if n == 0:
            cur = pm2
        elif n == 1:
            cur = pm1
        else:
            d = [pm1[i] - pm2[i] for i in range(M)]
            cur = [0] + d[:M - 1]
            pm2, pm1 = pm1, cur
        prod = None
        for i in range(M):
            if cur[i] != 0:
                if prod is None:
                    prod = Pa[n] * Pb[n]
                g[i] += cur[i] * prod
    return g


def ycoeffs(P):
    """fmpq_mpoly in (a, b, y) -> dict ey -> mpoly in (a, b) (y-free)."""
    out = {}
    for (ea, eb, ey), c in zip(P.monoms(), P.coeffs()):
        out.setdefault(ey, {})[(ea, eb, 0)] = c
    return {k: QC.from_dict(v) for k, v in out.items()}


if __name__ == '__main__':
    log = open('runs/Lsym.out', 'w')

    def say(*s):
        print(*s, flush=True)
        print(*s, file=log, flush=True)

    t0 = time.time()
    cs = interpolate_L()
    c4 = (y_ - 4) ** 2 * y_ ** 3 * zq(R4) * zq(R10)
    cs.append(c4)
    say('interpolated c0..c3 (%.0fs); total degrees in (a,b,y):' % (time.time() - t0),
        [c.total_degree() for c in cs], 'degrees (a,b,y):', [c.degrees() for c in cs])
    # c3 / c4 = d/dy log((y-4)^3 y^5 R4^3 / R10)  <=>  c3 * (y-4) y R4 R10 = c4 * (3 y R4 R10 + 5 (y-4) R4 R10 + 3 y (y-4) R4' R10 - y (y-4) R4 R10')
    r4, r10 = zq(R4), zq(R10)
    lhs = cs[3] * (y_ - 4) * y_ * r4 * r10
    rhs = cs[4] * (3 * y_ * r4 * r10 + 5 * (y_ - 4) * r4 * r10 + 3 * y_ * (y_ - 4) * r4.derivative('y') * r10
                   - y_ * (y_ - 4) * r4 * r10.derivative('y'))
    say('c3/c4 = CDT printed log-derivative (symbolic):', lhs == rhs)
    # spot check against solve_L at fresh points
    ok = True
    for (x, z) in [(101, 99), (-7, 3), (35, 37), (12, 35), (1001, 1003), (-21, -23)]:
        num = solve_L(x, z)
        ok &= all(spec(cs[i], x, z) == num[i] for i in range(5))
    say('agrees with solve_L at 6 fresh points:', ok)
    # rigorous check of [y^k] L~(G_A) = 0 for k < K0, identically in (a, b)
    t1 = time.time()
    g = g_series(K0 + 5)
    say('g_0..g_%d symbolic (%.0fs); max total degree %d' % (K0 + 4, time.time() - t1, g[K0 + 4].total_degree()))
    cy = [ycoeffs(c) for c in cs]
    t2 = time.time()
    allzero = True
    for k in range(K0):
        s = QC.from_dict({})
        for i in range(5):
            for j, cij in cy[i].items():
                m = k - j
                if m < 0:
                    continue
                f = 1
                for t in range(m + 1, m + i + 1):
                    f *= t
                s += cij * (f * g[m + i])
        if not s.is_zero():
            allzero = False
            say('row', k, 'NONZERO')
            break
        if k % 10 == 0:
            say('  row %d zero (%.0fs)' % (k, time.time() - t2))
    say('[y^k] L~(G_A) = 0 identically in Q[a,b] for k < %d:' % K0, allzero)
    _, rank = solve_L(fmpq(7), fmpq(5), K=K0, return_rank=True)
    say('rank of the first %d rows at (a, b) = (7, 5):' % K0, rank, '(unknowns: 63)')
    with open('runs/Lsym.txt', 'w') as f:
        for c in cs:
            f.write(str(c) + '\n')
    say('total %.0fs' % (time.time() - t0))
