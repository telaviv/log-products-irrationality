"""Phi6(a, b) by bivariate interpolation of
    r(a, b) = det(a, b) (a - b)^86 / [a^6 b^6 (a+b)^24 (4a^2-4ab+1)^4 (4b^2-4ab+1)^4 (4a^2-4ab+4b^2-3)^4
                                       (16a^2-40ab+16b^2+9)^2 (80a^2-128ab+80b^2-45)^2 Phi4 Phi14^3 Phi26 Phi79^3]
on a 7x7 grid (bidegree <= 6), checked on 12 extra points; Phi6 = primitive part.  The factor pattern is the one found
by f1_det_univ.py at b0 = 3, 5, 7, 9, 11 (identical there).  Output runs/phi6.txt.
"""
from flint import fmpq, fmpz, fmpq_mat
from lib_det import system
from lib_fac import load, ev
from lib_f1 import CTX, A_, B_

F = load()


def rest(a, b):
    a, b = fmpq(a), fmpq(b)
    v = a ** 6 * b ** 6 * (a + b) ** 24 * (4*a*a-4*a*b+1) ** 4 * (4*b*b-4*a*b+1) ** 4 * (4*a*a-4*a*b+4*b*b-3) ** 4
    v *= (16*a*a-40*a*b+16*b*b+9) ** 2 * (80*a*a-128*a*b+80*b*b-45) ** 2
    v *= ev(F['Phi4'], a, b) * ev(F['Phi14'], a, b) ** 3 * ev(F['Phi26'], a, b) * ev(F['Phi79'], a, b) ** 3
    return v


def r(a, b):
    M, _, _ = system(a, b)
    return M.det() * fmpq(a - b) ** 86 / rest(a, b)


if __name__ == '__main__':
    A = [17, 19, 23, 29, 31, 37, 41]
    B = [2, 3, 5, 7, 11, 13, 43]
    vals = {(x, y): r(x, y) for x in A for y in B}
    # solve for coefficients c_ij of a^i b^j, 0 <= i, j <= 6
    rows, rhs = [], []
    for (x, y), v in vals.items():
        rows.append([fmpq(x) ** i * fmpq(y) ** j for i in range(7) for j in range(7)])
        rhs.append([v])
    C = fmpq_mat(rows).solve(fmpq_mat(rhs))
    c = [C[k, 0] for k in range(49)]
    # check points
    chk = [(5, 9), (9, 5), (101, 99), (-7, 3), (61, 59), (3, -11), (15, 15 + 2), (1001, 1003), (8, 21), (-21, -23),
           (999999, 1000001), (12, 35)]
    ok = True
    for (x, y) in chk:
        pv = sum(c[i * 7 + j] * fmpq(x) ** i * fmpq(y) ** j for i in range(7) for j in range(7))
        ok &= (pv == r(x, y))
    print('extra points agree:', ok)
    den = fmpz(1)
    for q in c:
        den = den * q.q // den.gcd(q.q)
    ci = [int(q * den) for q in c]
    g = fmpz(0)
    for q in ci:
        g = g.gcd(q)
    P = CTX.from_dict({})
    for i in range(7):
        for j in range(7):
            if ci[i * 7 + j]:
                P += (ci[i * 7 + j] // int(g)) * A_ ** i * B_ ** j
    print('kappa (r = kappa * Phi6) =', fmpq(int(g)) / den)
    Ps = CTX.from_dict({(eb, ea, ey): cc for (ea, eb, ey), cc in zip(P.monoms(), P.coeffs())})
    print('symmetric:', Ps == P)
    print('Phi6 =', P)
    cont, fl = P.factor()
    print('irreducible over Q:', len(fl) == 1 and fl[0][1] == 1, 'content', cont)
    open('runs/phi6.txt', 'w').write(str(P) + '\n')
