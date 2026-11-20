"""rec.py -- exact common n-recurrence of Lysov's (Q_n, P1_n, P2_n, P3_n)(q) at a fixed integer q.

Q_n(z) = sum_l (-1)^(n+l) C(n+l,l) s_l z^(2l), s_l = sum_k C(n,l-k) C(2n+k,2k) C(n+k,k) C(4n+2k,2n+k)   (Lysov; n6_control.lysov)
P_j = polynomial part (in z) of Q_n * f_j,  f1 = log(1+1/z), f2 = log(1-1/z), f3 = theta = f1*f2.
At z = q:  P_j(q) = sum_i q_i * T_i,  T_i = sum_{k=1..i} c_k q^(i-k), f_j = sum_k c_k z^-k.

Method: residues mod many 62-bit primes; nullspace of the stacked (4 sequences) order/degree matrix; normalise the last
nonzero entry to 1; CRT + rational reconstruction; stop when two successive prime counts agree.  Then exact verification
over Q (fmpq) on exact values for n <= NV.  Output runs/rec_q<q>.json: p_j(n) coefficient lists (integers), j = 0..order.
Usage: python rec.py q [order] [deg] [NV]
"""
import sys, json, time, os
from math import comb
from fractions import Fraction
from flint import nmod_mat, fmpz, fmpq, fmpz_poly
HERE = os.path.dirname(os.path.abspath(__file__))

q = int(sys.argv[1]); ORD = int(sys.argv[2]) if len(sys.argv) > 2 else 4
DEG = int(sys.argv[3]) if len(sys.argv) > 3 else 34
NV = int(sys.argv[4]) if len(sys.argv) > 4 else 120
W = q * q

def logcoef(K):
    """c_k (k<=K) of f1, f2, theta in u=1/z, as Fractions."""
    a = [Fraction(0)] + [Fraction((-1) ** (k + 1), k) for k in range(1, K + 1)]
    b = [Fraction(0)] + [Fraction(-1, k) for k in range(1, K + 1)]
    t = [Fraction(0)] * (K + 1)
    for i in range(1, K + 1):
        for j in range(1, K + 1 - i):
            t[i + j] += a[i] * b[j]
    return a, b, t

def qcoefs(n):
    cs = [0] * (6 * n + 1)
    for l in range(3 * n + 1):
        s = 0
        for k in range(max(0, l - n), min(2 * n, l) + 1):
            s += comb(n, l - k) * comb(2 * n + k, 2 * k) * comb(n + k, k) * comb(4 * n + 2 * k, 2 * n + k)
        cs[2 * l] = (-1) ** (n + l) * comb(n + l, l) * s
    return cs

def exact_values(NMAX):
    K = 6 * NMAX
    a, b, t = logcoef(K)
    T = []
    for c in (a, b, t):
        Ti = [Fraction(0)] * (K + 1)
        for i in range(1, K + 1):
            Ti[i] = Ti[i - 1] * q + c[i]
        T.append(Ti)
    out = []
    for n in range(NMAX + 1):
        cs = qcoefs(n)
        Qv = sum(cs[i] * q ** i for i in range(len(cs)))
        Ps = [sum(cs[i] * Tj[i] for i in range(len(cs)) if cs[i]) for Tj in T]
        out.append((Qv, *Ps))
    return out

def matrix_modp(vals, p, nrows_per):
    unk = (ORD + 1) * (DEG + 1)
    rows = []
    for s in range(4):
        seq = [int(v[s].numerator % p) * pow(int(v[s].denominator % p), -1, p) % p if isinstance(v[s], Fraction) else v[s] % p
               for v in vals]
        for n in range(nrows_per):
            row = [0] * unk
            pw = 1
            for e in range(DEG + 1):
                for j in range(ORD + 1):
                    row[j * (DEG + 1) + e] = pw * seq[n + j] % p
                pw = pw * n % p
            rows.append(row)
    return nmod_mat(rows, p)

def ratrec(a, m):
    # rational reconstruction of a mod m
    a %= m
    r0, r1, s0, s1 = m, a, 0, 1
    import math
    B = math.isqrt(m // 2)
    while r1 > B:
        qq = r0 // r1
        r0, r1 = r1, r0 - qq * r1
        s0, s1 = s1, s0 - qq * s1
    if s1 == 0 or abs(s1) > B: return None
    return Fraction(r1, s1) if s1 > 0 else Fraction(-r1, -s1)

if __name__ == '__main__':
    t0 = time.time()
    unk = (ORD + 1) * (DEG + 1)
    nrows_per = unk // 4 + 12
    NMAX = max(nrows_per + ORD, NV + ORD)
    vals = exact_values(NMAX)
    print('exact values to n=%d (%.1fs)' % (NMAX, time.time() - t0), flush=True)
    # primes
    from sympy import prevprime
    p = 2 ** 62
    residues, mods = [], []
    prev = None
    sol = None
    M = 1
    for it in range(400):
        p = prevprime(p)
        A = matrix_modp(vals, p, nrows_per)
        X, nul = A.nullspace()
        if nul != 1:
            print('prime', p, 'nullity', nul); continue
        v = [int(X[i, 0]) for i in range(unk)]
        # normalise: last nonzero -> 1 (fix the index from the first prime)
        if not mods:
            piv = max(i for i in range(unk) if v[i])
        if v[piv] == 0:
            print('bad prime'); continue
        inv = pow(v[piv], -1, p)
        v = [x * inv % p for x in v]
        # CRT
        if not mods:
            acc = v; M = p
        else:
            acc = [(a + M * (((b - a) * pow(M, -1, p)) % p)) for a, b in zip(acc, v)]
            M *= p
        mods.append(p)
        rec = [ratrec(a, M) for a in acc]
        if all(r is not None for r in rec):
            if prev is not None and rec == prev:
                sol = rec; break
            prev = rec
        else:
            prev = None
    print('primes used', len(mods), '(%.1fs)' % (time.time() - t0), flush=True)
    import math
    den = 1
    for r in sol: den = den * r.denominator // math.gcd(den, r.denominator)
    ints = [int(r * den) for r in sol]
    g = 0
    for x in ints: g = math.gcd(g, x)
    ints = [x // g for x in ints]
    polys = [[ints[j * (DEG + 1) + e] for e in range(DEG + 1)] for j in range(ORD + 1)]
    # exact verification on all four sequences for n <= NV
    bad = 0
    for s in range(4):
        for n in range(0, NMAX - ORD + 1):
            tot = Fraction(0)
            for j in range(ORD + 1):
                pj = sum(polys[j][e] * n ** e for e in range(DEG + 1))
                tot += pj * Fraction(vals[n + j][s])
            if tot != 0: bad += 1
    print('exact verification n<=%d, 4 sequences: failures %d' % (NMAX - ORD, bad))
    mx = max(abs(x) for x in ints)
    print('max |coef| digits', len(str(mx)))
    for j in range(ORD + 1):
        print(' p_%d factor:' % j, str(fmpz_poly(polys[j]).factor())[:300])
    json.dump(dict(q=q, order=ORD, deg=DEG, polys=[[str(x) for x in pp] for pp in polys], verified_to=NMAX - ORD,
                   failures=bad, primes=len(mods)), open(os.path.join(HERE, 'runs', 'rec_q%d.json' % q), 'w'))
    print('done (%.1fs)' % (time.time() - t0))
