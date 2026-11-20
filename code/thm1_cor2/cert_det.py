"""cert_det.py -- from a mono13 pickle (ball strings, which enclose the computed balls): greedy choice of 14 words with
numerically independent continuation vectors (midpoints, mpmath), then the rigorous determinant of the 14 ball vectors
(python-flint acb at 2200 bits).  Usage: python cert_det.py pickle [gens]"""
import pickle, sys, json, mpmath as mp
from flint import acb, arb, acb_mat, ctx
ctx.prec = 2200; mp.mp.dps = 80
d = pickle.load(open(sys.argv[1], 'rb')); gens = sys.argv[2] if len(sys.argv) > 2 else None
def mpP(s):
    t = s.strip('[]').strip()
    return mp.mpf(0) if t.startswith('+/-') else mp.mpf(t.split(' +/-')[0])
Gb = {k: acb_mat([[acb(arb(a), arb(b)) for a, b in row] for row in M]) for k, M in d['G'].items() if gens is None or k in gens}
Gm = {k: mp.matrix([[mp.mpc(mpP(a), mpP(b)) for a, b in row] for row in M]) for k, M in d['G'].items() if k in Gb}
vb = acb_mat([[acb(arb(a), arb(b))] for a, b in d['v']]); vm = mp.matrix([[mp.mpc(mpP(a), mpP(b))] for a, b in d['v']])
n = 14
def nrm(x): return x / mp.sqrt(sum(abs(x[i]) ** 2 for i in range(n)))
def rank_ok(cols):
    M = mp.matrix(n, len(cols))
    for j, c in enumerate(cols):
        for i in range(n): M[i, j] = c[i]
    s = sorted([abs(t) for t in mp.svd_c(M, compute_uv=False)], reverse=True)
    return s[-1] / s[0] > mp.mpf(10) ** -45, s
basis = [('', vm)]; queue = [('', vm)]; tried = 0
while queue and len(basis) < n:
    w, x = queue.pop(0)
    for g in Gm:
        y = Gm[g] * x; tried += 1
        ok, s = rank_ok([nrm(b) for _, b in basis] + [nrm(y)])
        if ok:
            basis.append((g + w, y)); queue.append((g + w, y))
            print('dim %2d word %-10s smin/smax %s' % (len(basis), g + w, mp.nstr(s[-1] / s[0], 3)), flush=True)
            if len(basis) == n: break
words = [w for w, _ in basis]
print('words', words, 'tried', tried)
cache = {}
def ap(w):
    if w in cache: return cache[w]
    x = vb if w == '' else Gb[w[0]] * ap(w[1:]); cache[w] = x; return x
Mx = acb_mat(n, n)
for j, w in enumerate(words):
    x = ap(w)
    for i in range(n): Mx[i, j] = x[i, 0]
det = Mx.det()
print('CERTIFICATE det =', det.str(15), ' excludes 0:', not det.contains(0))
json.dump({'pickle': sys.argv[1], 'words': words, 'det': det.str(30), 'excludes_0': not det.contains(0)},
          open(sys.argv[1].replace('.pkl', '_det.json'), 'w'))
