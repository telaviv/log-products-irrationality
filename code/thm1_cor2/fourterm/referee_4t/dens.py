# empirical: log den(P_{j,n}(q))/n for j=1,2,3 at q=33, vs log y2; are P1,P2 denominators as small as P3's?
import sys, math
sys.argv = ['rec.py', '33']
from rec import exact_values
vals = exact_values(90)
from math import lcm
for n in (20, 40, 60, 80, 90):
    Q, P1, P2, P3 = vals[n]
    d = [P.denominator for P in (P1, P2, P3)]
    L = lcm(*range(1, 4*n+1)) * lcm(*range(1, 2*n+1))
    print(n, 'log den/n  P1 %.4f P2 %.4f P3 %.4f  common %.4f   log(lcm4n*lcm2n)/n %.4f' % (
        *[math.log(x)/n for x in d], math.log(lcm(*d))/n, math.log(L)/n))
print('log y2(33) = %.4f' % math.log(256.11942))
