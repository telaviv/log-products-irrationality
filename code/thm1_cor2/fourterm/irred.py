"""irred.py -- certified irreducibility of the monodromy action on S = Sol(L13) (13-dim, homogeneous solutions of the
minimal operator of g = y^5 H, coordinates = scaled Taylor coefficients at the base point), from a mono13 pickle
(generator matrices gA, g2, g3 as Arb ball strings, 14 x 14 on the affine module; v = data of g = y^5 H_3).

Criterion (reflection lemma; RESULT_4T sec. 3).  g2|_S = 1 + u (x) phi with u = (g2 - 1) v (in S, nonzero) [g2 is a
reflection on the affine module: proved from the local exponents 0..32, 1/2 of L at y2].  S is irreducible under
<gA, g2, g3> if (i) the vectors g_w u (13 words w) are independent and (ii) the covectors psi g_w (13 words) are
independent, psi = a nonzero row of (g2 - 1)|_S (a nonzero multiple of phi).  Both are 13 x 13 Arb determinants.

Also printed (sanity, not used): the second singular value of (g2 - 1)|_S, of (g3 - 1)|_S, and |g2^2 - 1|.
Usage: python irred.py pickle [out.json]"""
import pickle, sys, json, mpmath as mp
from flint import acb, arb, acb_mat, ctx
ctx.prec = 2200; mp.mp.dps = 120
d = pickle.load(open(sys.argv[1], 'rb'))
OUT = sys.argv[2] if len(sys.argv) > 2 else sys.argv[1].replace('.pkl', '_irred.json')
def mpP(s):
    t = s.strip('[]').strip()
    return mp.mpf(0) if t.startswith('+/-') else mp.mpf(t.split(' +/-')[0])
R = 13
Gfull = {k: acb_mat([[acb(arb(a), arb(b)) for a, b in row] for row in M]) for k, M in d['G'].items() if k in 'A23'}
for k, M in Gfull.items():      # affine structure: last row exactly (0, ..., 0, 1)
    assert all(M[R, j].contains(0) for j in range(R)) and M[R, R].contains(1), k   # exact by construction
G = {k: acb_mat([[M[i, j] for j in range(R)] for i in range(R)]) for k, M in Gfull.items()}
Gm = {k: mp.matrix([[mp.mpc(mpP(str(M[i, j].real)), mpP(str(M[i, j].imag))) for j in range(R)] for i in range(R)]) for k, M in G.items()}
v = acb_mat([[acb(arb(a), arb(b))] for a, b in d['v']])
gv = Gfull['2'] * v
u = acb_mat([[gv[i, 0] - v[i, 0]] for i in range(R)])
assert (gv[R, 0] - v[R, 0]).contains(0)
I = acb_mat(R, R)
for i in range(R): I[i, i] = acb(1)
D2 = G['2'] - I
# psi: the row of D2 at the index of the largest |u_i| (D2 = u phi^T, so row i0 = u_i0 phi^T)
i0 = max(range(R), key=lambda i: float(abs(u[i, 0]).mid()))
psi = acb_mat([[D2[i0, j] for j in range(R)]])
assert any(not psi[0, j].contains(0) for j in range(R)), 'psi row not certified nonzero'
um = mp.matrix([[mp.mpc(mpP(str(u[i, 0].real)), mpP(str(u[i, 0].imag)))] for i in range(R)])
psim = mp.matrix([[mp.mpc(mpP(str(psi[0, j].real)), mpP(str(psi[0, j].imag))) for j in range(R)]])
def nrm(x):
    s = mp.sqrt(sum(abs(t) ** 2 for t in x)); return [t / s for t in x]
def smin_ratio(cols):
    M = mp.matrix(R, len(cols))
    for j, c in enumerate(cols):
        for i in range(R): M[i, j] = c[i]
    s = sorted([abs(t) for t in mp.svd_c(M, compute_uv=False)], reverse=True)
    return s[-1] / s[0]
def greedy(start, act):
    """breadth-first words; keep a word if it raises the numerical rank (smin/smax > 1e-60)"""
    basis = [('', start)]; queue = [('', start)]; tried = 0
    while queue and len(basis) < R:
        w, x = queue.pop(0)
        for g in 'A23':
            y = act(g, x); tried += 1
            if smin_ratio([nrm(b) for _, b in basis] + [nrm(y)]) > mp.mpf(10) ** -60:
                basis.append((g + w, y)); queue.append((g + w, y))
                if len(basis) == R: break
    return [w for w, _ in basis], tried
colact = lambda g, x: list(Gm[g] * mp.matrix(x))
rowact = lambda g, x: list(mp.matrix([x]) * Gm[g])          # covector x -> x g
wu, tu = greedy(list(um), colact)
wp, tp = greedy(list(psim), rowact)
print('words (i)  :', wu, 'tried', tu)
print('words (ii) :', wp, 'tried', tp)
assert len(wu) == R and len(wp) == R
# rigorous determinants.  word 'abc' acting on a column: g_a g_b g_c x ; on a row: x g_a g_b g_c with the same convention
# as the greedy step (new letter applied last: row word g+w means (x g_w) g_g).
cache = {}
def apc(w):
    if w in cache: return cache[w]
    x = u if w == '' else G[w[0]] * apc(w[1:]); cache[w] = x; return x
cacher = {}
def apr(w):
    if w in cacher: return cacher[w]
    x = psi if w == '' else apr(w[1:]) * G[w[0]]; cacher[w] = x; return x
Mu = acb_mat(R, R); Mp = acb_mat(R, R)
for j, w in enumerate(wu):
    x = apc(w)
    for i in range(R): Mu[i, j] = x[i, 0]
for j, w in enumerate(wp):
    x = apr(w)
    for i in range(R): Mp[j, i] = x[0, i]
du, dp = Mu.det(), Mp.det()
print('det (i)  =', du.str(15), ' excludes 0:', not du.contains(0))
print('det (ii) =', dp.str(15), ' excludes 0:', not dp.contains(0))
# sanity (not used): rank-1 structure of g2 - 1, g3 - 1 on S, and g2^2 = 1
def svals(Mb):
    Mm = mp.matrix([[mp.mpc(mpP(str(Mb[i, j].real)), mpP(str(Mb[i, j].imag))) for j in range(R)] for i in range(R)])
    return sorted([abs(t) for t in mp.svd_c(Mm, compute_uv=False)], reverse=True)
s2 = svals(D2); s3 = svals(G['3'] - I); sq = svals(G['2'] * G['2'] - I)
print('sanity: svals(g2-1) top2 %s %s | svals(g3-1) top2 %s %s | max sval(g2^2-1) %s' % (
    mp.nstr(s2[0], 5), mp.nstr(s2[1], 5), mp.nstr(s3[0], 5), mp.nstr(s3[1], 5), mp.nstr(sq[0], 5)))
res = dict(pickle=sys.argv[1], q=d['q'], words_i=wu, words_ii=wp, det_i=du.str(30), det_ii=dp.str(30),
           excludes_0_i=not du.contains(0), excludes_0_ii=not dp.contains(0), psi_row=i0,
           sanity=dict(g2m1_sv12=[mp.nstr(s2[0], 6), mp.nstr(s2[1], 6)], g3m1_sv12=[mp.nstr(s3[0], 6), mp.nstr(s3[1], 6)],
                       g2sq_m1_svmax=mp.nstr(sq[0], 6)),
           irreducible_certified=(not du.contains(0)) and (not dp.contains(0)))
json.dump(res, open(OUT, 'w'), indent=1)
print('IRREDUCIBLE (certified):', res['irreducible_certified'])
