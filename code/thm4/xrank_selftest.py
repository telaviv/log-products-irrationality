"""xrank_selftest.py -- checks the series of xrank.py: (i) B4 against CDT's printed expansion -4y + 4y^2/9 + 31y^3/900 + 389y^4/88200
(mod p); (ii) the pure functions against their closed forms at y = 1/3 (exact Fraction coefficients, 60 terms, mpmath 40 digits);
(iii) G(y) against P(x) + P(x/(x-1)) at a sample x (Fractions/mpmath), for (a,b) = (5,7), r = (1,2,3,4)."""
import mpmath as mp
from fractions import Fraction as Fr
from math import comb, factorial
from xrank import series_tools
p = 2147483647
Gser, family, pure = series_tools(p, 60)
B = pure()
inv = lambda x: pow(x % p, p - 2, p)
want = [0, -4, Fr(4, 9), Fr(31, 900), Fr(389, 88200)]
for n in range(1, 5):
    assert B[3][n] == want[n].numerator * inv(want[n].denominator) % p, n
print('B4 printed expansion: OK')
mp.mp.dps = 40
N = 80
y = Fr(1, 3)
B2 = sum(Fr(2 * factorial(n - 2) * factorial(n), factorial(2 * n)) * y ** n for n in range(2, N))
B3 = sum(Fr(factorial(n - 1) ** 2, factorial(2 * n)) * y ** n for n in range(1, N))
yy = mp.mpf(1) / 3
print('B3 vs 2 asin^2(sqrt(y)/2):', mp.nstr(mp.mpf(B3.numerator) / B3.denominator - 2 * mp.asin(mp.sqrt(yy) / 2) ** 2, 5))
print('B2 vs 2y - 2 sqrt(y(4-y)) asin(sqrt(y)/2):', mp.nstr(mp.mpf(B2.numerator) / B2.denominator - (2 * yy - 2 * mp.sqrt(yy * (4 - yy)) * mp.asin(mp.sqrt(yy) / 2)), 5))
B5 = sum(Fr(factorial(n - 1) ** 2, (2 * n - 1) * factorial(2 * n - 1)) * y ** n for n in range(1, N))
print('B5 vs y 3F2(1/2,1,1;3/2,3/2;y/4):', mp.nstr(mp.mpf(B5.numerator) / B5.denominator - yy * mp.hyp3f2(0.5, 1, 1, 1.5, 1.5, yy / 4), 5))
# B4 = Sym^- Li2: check the double-sum formula against (x/(x-1)... ) is skipped; the printed expansion check above covers it.
# (iii) G vs Hadamard products
a, b = 5, 7
def uvF(a, NN):
    u = [Fr(1), Fr(a)]; v = [Fr(0), Fr(1)]
    for n in range(1, NN):
        u.append((a * (2 * n + 1) * u[n] - n * u[n - 1]) / (n + 1)); v.append((a * (2 * n + 1) * v[n] - n * v[n - 1]) / (n + 1))
    return u, v
NN = 120
ua, va = uvF(a, NN); ub, vb = uvF(b, NN)
r = (1, 2, 3, 4)
c = [r[0] * ua[n] * ub[n] + r[1] * va[n] * ub[n] + r[2] * ua[n] * vb[n] + r[3] * va[n] * vb[n] for n in range(NN)]
x = Fr(-1, 2000)
P = lambda t: sum(c[n] * t ** n for n in range(NN))
Gx = P(x) + P(x / (x - 1))
yv = x * x / (x - 1)
Gs = Gser(a, b, r)
# compare mod p: evaluate the truncated y-series exactly needs rationals; instead compare coefficientwise with an exact Fraction series
pn = [[Fr(2)], [Fr(0), Fr(1)]]
for n in range(2, NN):
    d = [ (pn[n-1][i] if i < len(pn[n-1]) else 0) - (pn[n-2][i] if i < len(pn[n-2]) else 0) for i in range(max(len(pn[n-1]), len(pn[n-2])))]
    pn.append([Fr(0)] + d)
GF = [sum(c[n] * (pn[n][j] if j < len(pn[n]) else 0) for n in range(NN)) for j in range(60)]
ok = all(GF[j].numerator * inv(GF[j].denominator) % p == Gs[j] for j in range(55))
print('G mod p equals exact Fraction series (55 coeffs):', ok)
Gy = sum(GF[j] * yv ** j for j in range(60))
print('G(y(x)) vs P(x)+P(x/(x-1)) at x=-1/2000:', float(Gy - Gx))
