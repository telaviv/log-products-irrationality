"""contact.py -- exact check: Q_n*theta - P3_n = O(z^-(2n+1)) (orthogonality k<2n) and Q_n*f_j - P_j orders, small n."""
import sys
from fractions import Fraction
sys.argv = ['rec.py', '32']
from rec import qcoefs, logcoef
for n in range(1, 13):
    cs = qcoefs(n); D = 6 * n; K = D + 2 * n + 6
    a, b, t = logcoef(K)
    out = []
    for c in (a, b, t):
        # Q(z) f(z), f = sum_k c_k z^-k ; coefficient of z^-s (s>=1) = sum_i q_i c_{i+s}
        first = None
        for s in range(1, 4 * n + 4):
            v = sum(Fraction(cs[i]) * c[i + s] for i in range(D + 1) if cs[i] and i + s <= K)
            if v != 0: first = s; break
        out.append(first)
    print('n=%d first nonzero z^-s in R1,R2,R3: %s  (2n+1=%d)' % (n, out, 2 * n + 1))
