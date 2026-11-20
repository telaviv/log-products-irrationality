"""K-3 closure, step C: proof of the determinant identity of RESULT.md section 5, step 3, by exact univariate computer
algebra on enough lines b = b0 (no interpolation).

From f1_detsym.py: det M(a, b) = det M'(t, b) / t^200, t = a - b, with M' polynomial and the a-priori bound
deg_b det M' <= 690 (in the coordinates (t, b); sum over columns of the largest b-degree of an entry).
Claim: det M (a - b)^86 = kappa * Pi(a, b)   <=>   P(t, b) := det M'(t, b) - kappa * t^114 * Pi(t + b, b) = 0,
Pi = a^6 b^6 (a+b)^24 (...)^4 (...)^2 Phi4 Phi6 Phi14^3 Phi26 Phi79^3 (f1_phi6_grid.rest times Phi6), kappa as in RESULT.md.
deg_b P <= max(690, totdeg Pi) =: DB.  For each of DB + 1 distinct integers b0 we compute det M'(t, b0) in Q[t] exactly
(fraction-free Bareiss over Z[t]) and compare it with kappa t^114 Pi(t + b0, b0).  Equality on DB + 1 lines gives P = 0.
Cross-check: det M'/t^200 agrees with lib_det.system (the code that produced the filed interpolation data) at sample points.
Usage: python f1_identity.py [b_start b_stop]   (results appended to runs/identity.json; resumable)
"""
import sys, json, os, time, pickle
from flint import fmpq, fmpz, fmpq_poly, fmpz_poly, fmpq_mpoly_ctx
from lib_fac import load
from lib_f1 import CTX, A_, B_

TC = fmpq_mpoly_ctx.get(('t', 'b', 'y'), 'lex')
t_, b_, y_ = TC.gens()
KAPPA = fmpq(645700815, 37778931862957161709568)
OUT = 'runs/identity.json'


def entry_tables(d):
    """each entry -> list over j of fmpq_poly in t (coefficient of b^j)."""
    cols = []
    for col in d['cols']:
        cc = []
        for dct in col:
            tab = {}
            for m, (cp, cq) in dct.items():
                tab.setdefault(m[1], {})[m[0]] = fmpq(cp, cq)
            polys = []
            for j in range(max(tab) + 1 if tab else 0):
                row = tab.get(j, {})
                polys.append(fmpq_poly([row.get(i, fmpq(0)) for i in range(max(row) + 1)]) if row else fmpq_poly([]))
            cc.append(polys)
        cols.append(cc)
    return cols


def spec_entry(polys, b0):
    r = fmpq_poly([])
    for p in reversed(polys):
        r = r * b0 + p
    return r


def bareiss(M):
    """det of an n x n matrix over Z[t] (lists of fmpz_poly), fraction-free, exact divisions asserted."""
    n = len(M)
    M = [row[:] for row in M]
    sign = 1
    prev = fmpz_poly([1])
    for k in range(n - 1):
        if M[k][k] == 0:
            sw = next((r for r in range(k + 1, n) if M[r][k] != 0), None)
            if sw is None:
                return fmpz_poly([])
            M[k], M[sw] = M[sw], M[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                q, r = divmod(M[k][k] * M[i][j] - M[i][k] * M[k][j], prev)
                assert r == 0
                M[i][j] = q
        prev = M[k][k]
    return sign * M[n - 1][n - 1]


def pi_poly(F, PHI6, b0):
    """kappa * t^114 * Pi(t + b0, b0) as fmpq_poly in t."""
    T = fmpq_poly([0, 1])
    a = T + b0
    b = fmpq(b0)

    def ev(P):
        # fmpz_mpoly in (a, b, y), y-free -> fmpq_poly in t
        r = fmpq_poly([])
        for (ea, eb, ey), c in zip(P.monoms(), P.coeffs()):
            r += fmpq(int(c)) * b ** eb * a ** ea
        return r
    v = a ** 6 * b ** 6 * (a + b) ** 24 * (4*a*a-4*a*b+1) ** 4 * (4*b*b-4*a*b+1) ** 4 * (4*a*a-4*a*b+4*b*b-3) ** 4
    v *= (16*a*a-40*a*b+16*b*b+9) ** 2 * (80*a*a-128*a*b+80*b*b-45) ** 2
    v *= ev(F['Phi4']) * ev(PHI6) * ev(F['Phi14']) ** 3 * ev(F['Phi26']) * ev(F['Phi79']) ** 3
    return KAPPA * T ** 114 * v


if __name__ == '__main__':
    d = pickle.load(open('runs/detsym.pkl', 'rb'))
    assert d['E'] == 200
    cols = entry_tables(d)
    degb = sum(max(len(p) - 1 for p in col) for col in cols)
    F = load()
    PHI6 = eval(open('runs/phi6.txt').read().replace('^', '**'), {'a': A_, 'b': B_})
    PiTot = sum(F[k].total_degree() * e for k, e in (('Phi4', 1), ('Phi14', 3), ('Phi26', 1), ('Phi79', 3)))
    PiTot += PHI6.total_degree() + 12 + 24 + 8 + 8 + 8 + 4 + 4  # a^6 b^6, (a+b)^24, the five quadratic powers
    degb, PiTot = int(degb), int(PiTot)
    DB = max(degb, PiTot)
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    res.setdefault('meta', {})
    res['meta'].update({'deg_b_bound_detMp': degb, 'totdeg_Pi': PiTot, 'DB': DB, 'lines_needed': DB + 1})
    print('deg_b(det M\') <=', degb, ' totdeg Pi =', PiTot, ' -> need', DB + 1, 'lines', flush=True)
    lo, hi = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (-(DB // 2) - 1, DB - DB // 2 + 1)
    # cross-check with lib_det at a few points (only on the first run)
    if 'xcheck' not in res['meta']:
        from lib_det import system
        ok = True
        for (x, z) in [(7, 5), (-9, 3), (101, 99), (2001, 2003)]:
            tt = fmpq(x - z)
            dm = system(x, z)[0].det()
            col = [[spec_entry(e, fmpq(z))(tt) for e in c] for c in cols]
            from flint import fmpq_mat
            dsym = fmpq_mat([[col[k][i] for k in range(10)] for i in range(10)]).det() / tt ** 200
            ok &= (dsym == dm)
        res['meta']['xcheck'] = bool(ok)
        print('symbolic det M = lib_det det at 4 points:', ok, flush=True)
    lines = res.setdefault('lines', {})
    t0 = time.time()
    for b0 in range(lo, hi):
        if str(b0) in lines:
            continue
        colsp = []
        scal = fmpz(1)
        for c in cols:
            ps = [spec_entry(e, fmpq(b0)) for e in c]
            den = fmpz(1)
            for p in ps:
                dd = p.denom()
                den = den * dd // den.gcd(dd)
            scal *= den
            colsp.append([(p * den).numer() if p != 0 else fmpz_poly([]) for p in ps])
        M = [[colsp[k][i] for k in range(10)] for i in range(10)]
        D = fmpq_poly(bareiss(M)) / scal
        ok = (D == pi_poly(F, PHI6, b0))
        lines[str(b0)] = {'ok': bool(ok), 'deg': int(D.degree())}
        if not ok:
            print('MISMATCH at b0 =', b0, flush=True)
        if len(lines) % 25 == 0:
            json.dump(res, open(OUT, 'w'))
            print(len(lines), 'lines,', b0, '%.0fs' % (time.time() - t0), flush=True)
    json.dump(res, open(OUT, 'w'))
    n_ok = sum(1 for v in lines.values() if v['ok'])
    print('lines checked:', len(lines), ' all equal:', n_ok == len(lines), ' needed:', DB + 1,
          ' PROVED:' if n_ok == len(lines) and len(lines) >= DB + 1 else ' not yet', flush=True)
