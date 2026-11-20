# magnitude survey of the quantities that mono13.circle_max compares as Python floats (overflow/underflow audit)
import sys, os, math, cmath
q = int(sys.argv[1]); HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.argv = ['mono13.py', str(q), '400']
import mono13 as M
from flint import acb, arb, ctx
ctx.prec = 256
R = M.R
lam = sorted([1/abs(complex(r)) for r in []]) if False else None
yt = sorted([z.real for z in M.rootsf if abs(z.imag) < 1e-30 and z.real > 0])
print('positive real singular points', ['%.4g' % x for x in yt])
cs = []
for x in [2**-21 * 1.2 ** k for k in range(0, 90)]: cs.append(complex(x, 0))
base = 0.42 * max(t for t in yt if t < 1000) / 1.17  # roughly
for k in range(60): cs.append(cmath.rect(66 if q == 26 else 101, 2 * math.pi * k / 60 + 0.01))
for x in yt[-3:]:
    for k in range(24): cs.append(x + 0.3 * x / 10 * cmath.exp(1j * (2 * math.pi * k / 24)))
lo, hi, vlo, vhi = 1e999, -1e999, 1e999, -1e999
for c in cs:
    ca = acb(c.real, c.imag)
    d = M.delta_of(ca); sig = arb(d.lower()) * 0.35
    if float(sig) <= 0: continue
    lp = ctx.prec
    for k in range(0, M.NARC, 16):
        ang = 2 * math.pi * (2 * k + 1) / (2 * M.NARC)
        y = ca + acb(sig * arb(ang).cos(), sig * arb(ang).sin())
        bR = abs(arb(M.Bc[R][M.RHO]))
        for r in M.roots: bR = bR * abs(y - r)
        ar = bR * abs(y) ** R
        l = float(ar.mid().log()) / math.log(10); lo = min(lo, l); hi = max(hi, l)
        for i in range(R):
            num = M.Blc[i]
            for r in M.Broots[i]: num = num * abs(y - r)
            v = num / bR * abs(y) ** (i - R)
            if v.mid() != 0:
                lv = float(v.mid().log()) / math.log(10); vlo = min(vlo, lv); vhi = max(vhi, lv)
print('q=%d log10|A_R| on circles in [%.1f, %.1f]; log10 M_i values in [%.1f, %.1f]' % (q, lo, hi, vlo, vhi))
