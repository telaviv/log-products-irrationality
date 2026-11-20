"""budget_4t.py -- m_needed(q) for q = 26..53 (Arb).  lambda2 = second root of the lambda-quartic with W = q^2 (= the exact
characteristic polynomial of the order-4 recurrence for every q >= 2, PROOF_I1 sec. 7; roots real, distinct, one sign).
L = log(4 y2), y2 = 1/lambda2; tau = 6 (m^2 - 1)/m^2 (rows (0, 6, ..., 6), e = 0); M(m) = (1 - 1/m) L - tau.
m_needed(q) = least m with M(m) > 0 (certified: the Arb ball of M lies in (0, oo)); M(m) is increasing in m once positive
iff L > 6 (M(m) > 0  <=>  L > 6 (1 + 1/m)), so m_needed = floor(6/(L-6)) + 1, and the certificate is M(m_needed) > 0.
Writes budget_4t.json.  Usage: python budget_4t.py"""
import json
from flint import fmpz_poly, arb, ctx
ctx.prec = 300
out = {}
for q in range(26, 54):
    W = q * q
    # quartic(lam) coefficients (geom.py); the recurrence's char. polynomial is const * quartic(-lam)
    quart = [65536, 294912 * W - 16384, 135168 * W ** 2 + 159744 * W + 1536, 16384 * W ** 3 - 20480 * W ** 2 + 5248 * W - 64, 1]
    ch = fmpz_poly([c * (-1) ** k for k, c in enumerate(quart)])          # quartic(-lam)
    rs = [r for r, m in ch.complex_roots()]
    assert len(rs) == 4 and all(r.imag.contains(0) for r in rs)
    lam = sorted([r.real for r in rs], key=lambda x: -float(x.mid()))
    assert all(l > 0 for l in lam) and all(lam[i] > lam[i + 1] for i in range(3))
    L = arb(4).log() - lam[1].log()
    assert L > 6
    mn = int(float((6 / (L - 6)).mid())) + 1
    M = lambda m: L * (1 - arb(1) / m) - arb(6) * (m * m - 1) / (m * m)
    assert M(mn) > 0 and not (M(mn - 1) > 0)
    out[q] = dict(L=L.str(12), m_needed=mn, M_at_m_needed=M(mn).str(8), M_at_m_minus_1=M(mn - 1).str(8),
                  M14=M(14).str(8), M13=M(13).str(8), y2=(1 / lam[1]).str(12), y3_over_y2=(lam[1] / lam[2]).str(6))
    print(q, 'L=%s  m_needed=%d  M(m_needed)=%s  M(m-1)=%s  M(14)=%s' % (L.str(10), mn, M(mn).str(6), M(mn - 1).str(6), M(14).str(6)))
json.dump(out, open('budget_4t.json', 'w'), indent=1)
