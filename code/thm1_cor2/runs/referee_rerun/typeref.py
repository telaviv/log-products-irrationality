# independent check of den(P_{1,n}), den(P_{3,n}) | lcm(1..4n) lcm(1..2n), from Lysov (21) and the P1 formula, as polynomials in x
from math import comb, lcm
from fractions import Fraction
import sys
def D(k):
    r = 1
    for i in range(1, k + 1): r = lcm(r, i)
    return r
def c(n, l):
    s = sum(comb(n, l - k) * comb(2*n + k, 2*k) * comb(n + k, k) * comb(4*n + 2*k, 2*n + k) for k in range(max(0, l - n), min(2*n, l) + 1))
    return (-1) ** (n + l) * comb(n + l, l) * s
for n in map(int, sys.argv[1:]):
    cs = [c(n, l) for l in range(3*n + 1)]
    tau = [None] + [sum(Fraction((-1) ** k, k * (2*m - k)) for k in range(1, 2*m)) for m in range(1, 3*n + 1)]
    P3 = {}
    for l in range(1, 3*n + 1):
        for m in range(1, l + 1):
            P3[l - m] = P3.get(l - m, 0) + cs[l] * tau[m]
    P1 = {}
    for l in range(1, 3*n + 1):
        for m in range(1, 2*l + 1):
            P1[2*l - m] = P1.get(2*l - m, 0) + Fraction(cs[l], m)
    den3 = 1; den1 = 1
    for v in P3.values(): den3 = lcm(den3, Fraction(v).denominator)
    for v in P1.values(): den1 = lcm(den1, Fraction(v).denominator)
    M = D(4*n) * D(2*n)
    print('n=%d den(P3) | D4n D2n: %s   den(P1) | D4n D2n: %s   log den3/n = %.3f' % (n, M % den3 == 0, M % den1 == 0, __import__('math').log(den3) / n))
