"""lemmaS2.py (pilot_F2; copy of range_pilot/lemmaS_check.py with eps = 1e-5, 3e-5, 1e-4) -- -- numerical spot check (mpmath, 60 digits) of Lemma S: for m != n in Z minus {-1,0}, |m/n - 1| < eps,
a = 2m+1, b = 2n+1, the three genuine singular values y = x^2/(x-1) (x in {uv, u/v, v/u}, u = a + sgn(a) sqrt(a^2-1)) are real
with |y| >= 1/eps - 6, and alpha (x = 1/(uv)) has |alpha| <= eps^4/16.  Worst cases are |n| just above 1/eps.  [Numerical]"""
import mpmath as mp
mp.mp.dps = 60
worst = None
for eps in (mp.mpf('1e-5'), mp.mpf('3e-5'), mp.mpf('1e-4')):
    N0 = int(1 / eps)
    for n in list(range(N0 + 1, N0 + 4)) + [3 * N0 + 1, 10 * N0 + 7]:
        for sgn in (1, -1):
            for d in (-1, 1, 2, -2, int(eps * n), -int(eps * n)):
                nn = sgn * n; mm = nn + sgn * d
                if mm in (0, -1) or mm == nn or abs(mp.mpf(mm) / nn - 1) >= eps:
                    continue
                a, b = 2 * mm + 1, 2 * nn + 1
                u = a + mp.sign(a) * mp.sqrt(a * a - 1); v = b + mp.sign(b) * mp.sqrt(b * b - 1)
                ys = [x * x / (x - 1) for x in (u * v, u / v, v / u)]
                al = (1 / (u * v)) ** 2 / (1 / (u * v) - 1)
                r = min(abs(y) for y in ys) - (1 / eps - 6)
                ra = abs(al) / (eps ** 4 / 16)
                if worst is None or r < worst[0]:
                    worst = (r, float(eps), mm, nn)
                assert r > 0 and ra < 1, (eps, mm, nn, r, ra)
print('Lemma S spot check passed; min over samples of (min|y| - (1/eps - 6)) =', mp.nstr(worst[0], 8), 'at', worst[1:])
