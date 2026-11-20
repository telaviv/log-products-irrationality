"""i1_opalg.py -- exact operator algebra: from (R1) D_n Q_n = 0 and (R2) Q_{n+1} = alpha(n)^{-1} beta(n,w,theta) Q_n
derive the order-4 recurrence  L4 = sum_{i=0}^{4} p_i(n,w) S_n^i,  p_i in Z[n,w],  with  L4 Q_n(w) = 0  for every n >= 0,
symbolically in w = z^2.  Then compare L4(n, q^2) with the pilot's runs/rec_q<q>.json for q = 26..32.

Algebra: R = Q(n,w)<theta>, theta w = w (theta + 1).  N_0 = 1,  N_{i+1} = ( E_{n+i} o N_i ) rem_right D_n, where
E_m = alpha(m)^{-1} beta(m,w,theta).  Then Q_{n+i} = N_i Q_n (i = 0..4), N_i of theta-order <= 3.  p = cofactor vector
(signed 4x4 minors) of the 4x5 matrix [N_0 .. N_4]; so sum_i p_i N_i = 0 identically, hence L4 Q_n = 0.
Denominators that occur: alpha(n+t) (t = 0..3; no root n >= 0) and powers of (w-1) (the leading coefficient of D_n is
-2(w-1)^2, free of n); the (w-1)-powers cancel exactly (every N_i is printed with E = 0, i.e. polynomial in w).  So the
identity specialises to every integer n >= 0 as an identity of polynomials in w.

Output: runs/i1_L4_symbolic.json (p_i as {(i_n, i_w): coeff}), and a comparison line per q.
Usage: python i1_opalg.py
"""
import json, os, time
from math import comb
import sympy as sp
from flint import fmpq_mpoly_ctx, fmpq, fmpz_poly, fmpq_poly

HERE = os.path.dirname(os.path.abspath(__file__))
ctx = fmpq_mpoly_ctx.get(('n', 'w'), 'lex')
N, W = ctx.gens()
one = ctx.from_dict({(0, 0): 1})
zero = ctx.from_dict({})
rel = json.load(open(os.path.join(HERE, 'runs', 'i1_relations.json')))
ns, ls, ws = sp.symbols('n l w')


def to_mp(expr):
    P = sp.Poly(sp.expand(expr), ns, ws)
    return ctx.from_dict({(int(m[0]), int(m[1])): fmpq(int(sp.Rational(c).p), int(sp.Rational(c).q)) for m, c in P.terms()}) \
        if not P.is_zero else zero


def shift_n(p, t):
    """p(n,w) -> p(n+t,w)"""
    return p.compose(N + t, W)


def theta(p):
    return W * p.derivative('w')


# ---- D_n = sum_m w^{2-m} a_m(n, theta-m) as sum_a delta_a theta^a
th = sp.Symbol('th')
a_list = [sp.sympify(x) for x in rel['a']]
Dop = {}
for m in range(3):
    e = sp.expand(a_list[m].subs(ls, th - m))
    P = sp.Poly(e, th)
    for (deg,), cf in P.terms():
        Dop[deg] = Dop.get(deg, zero) + to_mp(cf * ws ** (2 - m))
Dop = [Dop.get(i, zero) for i in range(5)]
assert Dop[4] == to_mp(-2 * (ws - 1) ** 2), Dop[4]
alpha = to_mp(sp.sympify(rel['alpha']))
beta = [to_mp(sp.sympify(x)) for x in rel['beta']]
wm1 = W - 1


def theta_pow_on(h, E, j):
    """theta^j applied to h/(w-1)^E  ->  numerator over (w-1)^(E+j)"""
    for _ in range(j):
        h = theta(h) * wm1 - E * W * h
        E += 1
    return h


def compose(Bop, Nop, E):
    """Bop = sum_b B_b theta^b (polynomial coeffs), Nop = (w-1)^{-E} sum_a eta_a theta^a.
    returns (list of numerators, new exponent) with result = (w-1)^{-E'} sum_c res_c theta^c"""
    bmax = len(Bop) - 1
    E2 = E + bmax
    res = {}
    for b, Bb in enumerate(Bop):
        if Bb == 0:
            continue
        for a, eta in enumerate(Nop):
            if eta == 0:
                continue
            for j in range(b + 1):
                h = theta_pow_on(eta, E, j) * wm1 ** (bmax - j)   # over (w-1)^(E+bmax)
                c = a + b - j
                res[c] = res.get(c, zero) + comb(b, j) * Bb * h
    top = max(res)
    return [res.get(c, zero) for c in range(top + 1)], E2


def reduce_D(op, E, Dn):
    """op = (w-1)^{-E} sum op_c theta^c  ->  right remainder mod Dn (order 4, lead -2(w-1)^2)."""
    op = list(op)
    while len(op) > 4:
        c = len(op) - 1
        h = op[c]
        if h == 0:
            op.pop(); continue
        s = c - 4
        # theta^s o Dn  as list
        TD = {}
        for a, da in enumerate(Dn):
            for j in range(s + 1):
                dj = da
                for _ in range(j):
                    dj = theta(dj)
                TD[a + s - j] = TD.get(a + s - j, zero) + comb(s, j) * dj
        new = [x * wm1 ** 2 for x in op]
        for cc, val in TD.items():
            new[cc] = new[cc] + fmpq(1, 2) * h * val   # subtract h/(-2(w-1)^2) * (theta^s Dn), scaled by (w-1)^2
        assert new[c] == 0
        new.pop()
        op = new
        E += 2
    while len(op) < 4:
        op.append(zero)
    return op, E


if __name__ == '__main__':
    t0 = time.time()
    Ns = [([one, zero, zero, zero], 0, one)]   # (numerators, E, n-denominator)
    for i in range(4):
        eta, E, dn = Ns[-1]
        Bi = [shift_n(b, i) for b in beta]
        Dn = Dop
        comp, E2 = compose(Bi, eta, E)
        red, E3 = reduce_D(comp, E2, Dn)
        # remove common (w-1) factors
        while E3 > 0 and all(x == 0 or x.compose(N, ctx.from_dict({(0, 0): 1})) == 0 for x in red):
            red = [x / wm1 if x != 0 else zero for x in red]
            E3 -= 1
        Ns.append((red, E3, dn * shift_n(alpha, i)))
        print('N_%d: E=%d, n-deg %s, w-deg %s  (%.1fs)' % (i + 1, E3, max(x.degrees()[0] for x in red if x != 0),
                                                       max(x.degrees()[1] for x in red if x != 0), time.time() - t0), flush=True)
    # common denominator: Dall = dn_4 * (w-1)^Emax ; column i scaled
    Emax = max(E for _, E, _ in Ns)
    d4 = Ns[4][2]
    cols = []
    for eta, E, dn in Ns:
        # N_i = eta / (dn (w-1)^E) ;  multiply by d4 (w-1)^Emax
        fac_n = d4 / dn   # exact polynomial division in n
        cols.append([x * fac_n * wm1 ** (Emax - E) for x in eta])
    M = [[cols[i][a] for i in range(5)] for a in range(4)]

    def det(mat):
        k = len(mat)
        if k == 1:
            return mat[0][0]
        s = zero
        for j in range(k):
            if mat[0][j] == 0:
                continue
            minor = [row[:j] + row[j + 1:] for row in mat[1:]]
            term = mat[0][j] * det(minor)
            s = s + term if j % 2 == 0 else s - term
        return s
    p = []
    for i in range(5):
        sub = [row[:i] + row[i + 1:] for row in M]
        d = det(sub)
        p.append(d if i % 2 == 0 else -d)
    print('minors done (%.1fs)' % (time.time() - t0), flush=True)
    g = p[0]
    for x in p[1:]:
        g = g.gcd(x)
    p = [x / g for x in p]
    # make integral and primitive
    den = 1
    from math import lcm, gcd
    for x in p:
        for cf in x.coeffs():
            den = lcm(den, int(cf.q))
    p = [x * den for x in p]
    cg = 0
    for x in p:
        for cf in x.coeffs():
            cg = gcd(cg, int(cf.p))
    p = [x / cg for x in p]
    # verify  sum_i p_i N_i = 0  exactly
    for a in range(4):
        s = zero
        for i in range(5):
            s = s + p[i] * M[a][i]
        assert s == 0
    print('identity sum p_i N_i = 0 verified exactly')
    print('content removed:', g.degrees(), ' p_i degrees (n,w):', [x.degrees() for x in p])
    dw = int(max(x.degrees()[1] for x in p))
    out = dict(p=[{'%d,%d' % m: str(cf) for m, cf in zip(x.monoms(), x.coeffs())} for x in p], deg_w=dw,
               deg_n=int(max(x.degrees()[0] for x in p)), removed_gcd_degrees=[int(t) for t in g.degrees()])
    json.dump(out, open(os.path.join(HERE, 'runs', 'i1_L4_symbolic.json'), 'w'))
    # ---- compare with the pilot's per-q recurrences
    comps = []
    for q in range(26, 33):
        d = json.load(open(os.path.join(HERE, 'runs', 'rec_q%d.json' % q)))
        P = [[int(x) for x in pp] for pp in d['polys']]
        spec = []
        for x in p:
            y = x.compose(N, ctx.from_dict({(0, 0): q * q}))
            # univariate in n
            coeffs = [0] * (y.degrees()[0] + 1)
            for m, cf in zip(y.monoms(), y.coeffs()):
                coeffs[m[0]] = cf
            spec.append(fmpq_poly(coeffs))
        rec = [fmpq_poly(pp) for pp in P]
        # proportional over Q(n)?  spec_i * rec_4 == rec_i * spec_4
        prop = all(spec[i] * rec[4] == rec[i] * spec[4] for i in range(5))
        ratio = spec[4] / rec[4] if prop else None
        const = None
        if prop:
            rr = spec[4] / rec[4]          # fmpq_poly division -> must be exact constant for 'same recurrence'
            const = rr if (spec[4] == rec[4] * rr) else None
        print('q=%d: L4(n,q^2) proportional to rec_q%d.json over Q(n): %s; constant factor: %s' % (
            q, q, prop, const), flush=True)
        comps.append(dict(q=q, proportional=prop, constant=str(const)))
    out['compare_rec_q'] = comps
    json.dump(out, open(os.path.join(HERE, 'runs', 'i1_L4_symbolic.json'), 'w'))
    print('done (%.1fs)' % (time.time() - t0))
