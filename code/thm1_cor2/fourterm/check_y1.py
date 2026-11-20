"""check_y1.py -- the claim of Section 3.4 of the paper: for every integer 2 <= q <= 60 the roots of chi_{q^2}(-lambda)
(the characteristic polynomial of the order-4 recurrence, Proposition 3.4) are real, positive and distinct, and
y_1 = 1/lambda_1 < 2^-20 q^2, so that H_j is holomorphic at y_1 (Corollary 3.2).  Arb (python-flint), 300 bits.
Added by the final referee (2026-11-19); the same check was done by the gate referee of Theorem 1 in scratch.
Usage: python check_y1.py"""
from flint import fmpz_poly, arb, ctx
ctx.prec = 300
worst = None
for q in range(2, 61):
    W = q * q
    quart = [65536, 294912 * W - 16384, 135168 * W ** 2 + 159744 * W + 1536, 16384 * W ** 3 - 20480 * W ** 2 + 5248 * W - 64, 1]
    ch = fmpz_poly([c * (-1) ** k for k, c in enumerate(quart)])          # chi_w(-lambda)
    rs = [r for r, m in ch.complex_roots()]
    assert len(rs) == 4 and all(m == 1 for r, m in ch.complex_roots())
    assert all(r.imag.contains(0) for r in rs)                            # real (Arb root isolation of an integer poly)
    lam = sorted([r.real for r in rs], key=lambda x: -float(x.mid()))
    assert all(l > 0 for l in lam) and all(lam[i] > lam[i + 1] for i in range(3))
    y1 = 1 / lam[0]
    bound = arb(q * q) / 2 ** 20
    assert y1 < bound
    ratio = y1 / bound
    if worst is None or float(ratio.upper()) > worst[0]:
        worst = (float(ratio.upper()), q)
print('q = 2..60: roots real, positive, distinct; y_1 < 2^-20 q^2 at every q; max y_1/(2^-20 q^2) = %.3e at q = %d'
      % worst)
