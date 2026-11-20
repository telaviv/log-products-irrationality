import sys, math
sys.argv = ['rec.py', '33']
from rec import exact_values
from sympy import factorint
vals = exact_values(120)
for n in (40, 60, 90, 100, 110, 120):
    Q, P1, P2, P3 = vals[n]
    d3 = P3.denominator; c = math.lcm(P1.denominator, P2.denominator, d3)
    print(n, 'extra factor lcm/den(P3) =', factorint(c // d3), ' rate %.4f' % (math.log(c // d3) / n))
