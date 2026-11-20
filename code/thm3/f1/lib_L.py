"""Numeric (exact rational) reconstruction of CDT's operator L = sum c_i(y) D^i at a given (a, b).

Normalisation (CDT (14.2.1)): c4 = (y-4)^2 y^3 R4(y) R10(y).  Ansatz: deg c0 <= 12, deg c1 <= 16, c2 = y*(deg <= 16),
c3 = y^2 (y-4) * (deg <= 15).  Unknowns are fixed by L(G_A) = 0 mod y^K (overdetermined; consistency and uniqueness
are checked by rank).
"""
from flint import fmpq, fmpq_poly, fmpq_mat
from lib_f1 import R10, R4, spec_y, GA_series

Y = fmpq_poly([0, 1])


def solve_L(a, b, K=110, return_rank=False):
    a, b = fmpq(a), fmpq(b)
    g = GA_series(a, b, K + 6)
    r10, r4 = spec_y(R10, a, b), spec_y(R4, a, b)
    c4 = (Y - 4) ** 2 * Y ** 3 * r4 * r10
    # derivative coefficient lists: [y^m] G^(i) = g_{m+i} (m+i)!/m!
    def dcoef(i, m):
        if m < 0:
            return fmpq(0)
        f = fmpq(1)
        for t in range(m + 1, m + i + 1):
            f *= t
        return g[m + i] * f
    # unknown blocks: (i, list of basis polys)
    blocks = []
    blocks.append((0, [Y ** j for j in range(13)]))
    blocks.append((1, [Y ** j for j in range(17)]))
    blocks.append((2, [Y ** (j + 1) for j in range(17)]))
    blocks.append((3, [Y ** 2 * (Y - 4) * Y ** j for j in range(16)]))
    cols = []
    for i, basis in blocks:
        for P in basis:
            cols.append((i, P))
    nunk = len(cols)
    rows = []
    for k in range(K):
        row = []
        for i, P in cols:
            s = fmpq(0)
            for j in range(P.degree() + 1):
                cj = P[j]
                if cj != 0:
                    s += cj * dcoef(i, k - j)
            row.append(s)
        s = fmpq(0)
        for j in range(c4.degree() + 1):
            if c4[j] != 0:
                s += c4[j] * dcoef(4, k - j)
        row.append(-s)
        rows.append(row)
    M = fmpq_mat(rows)
    R, rank = M.rref()
    # consistency: no pivot in last column
    sol = [fmpq(0)] * nunk
    piv_cols = []
    for r in range(rank):
        for c in range(nunk + 1):
            if R[r, c] != 0:
                piv_cols.append(c)
                break
    assert nunk not in piv_cols, 'inconsistent ansatz'
    assert rank == nunk, ('non-unique', rank, nunk)
    for r in range(rank):
        sol[piv_cols[r]] = R[r, nunk]
    cs = [fmpq_poly([]) for _ in range(4)]
    for (i, P), s in zip(cols, sol):
        cs[i] += s * P
    cs.append(c4)
    if return_rank:
        return cs, rank
    return cs
