"""Exact factorisation of Disc_y(R10), Res_y(R4,R10), Res_y(R10, y(y-4)) (python-flint fmpz_mpoly).
Writes runs/factors.txt: one irreducible factor per line as 'mult | deg_a | total_deg | poly' (poly in a, b)."""
import os, time
from flint import fmpz_mpoly_ctx
from lib_f1 import R10, R4, CTX, Y_

os.makedirs(os.path.join(os.path.dirname(__file__), 'runs'), exist_ok=True)
out = open(os.path.join(os.path.dirname(__file__), 'runs', 'factors.txt'), 'w')
allf = {}
for name, P in [('Disc(R10)', lambda: R10.discriminant('y')),
                ('Res(R4,R10)', lambda: R4.resultant(R10, 'y')),
                ('Res(R10,y(y-4))', lambda: R10.resultant(Y_ * (Y_ - 4), 'y'))]:
    t = time.time()
    c, fl = P().factor()
    print(name, 'content', c, 'time %.1f' % (time.time() - t), flush=True)
    out.write('# %s content %s\n' % (name, c))
    for f, e in fl:
        da = max(m[0] for m in f.monoms())
        print('   deg_a %d total %d mult %d %s' % (da, f.total_degree(), e, str(f) if f.total_degree() <= 4 else ''))
        out.write('%d | %d | %d | %s\n' % (e, da, f.total_degree(), str(f)))
out.close()
