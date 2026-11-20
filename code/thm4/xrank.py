"""xrank.py -- falsifier for the cross-pair independence lemma (pilot_F2 step 1).
Exact mod-p rank of the union {B1..B7} + {10 G-family functions of each pair}, pairs given on the command line.
For pair (a, b) (odd a, b; a = 2m+1, b = 2n+1) and a relation-type vector r = (r0, ra, rb, rab):
  c_n = r0 u_n(a)u_n(b) + ra v_n(a)u_n(b) + rb u_n(a)v_n(b) + rab v_n(a)v_n(b)   (CDT (14.1.8); u = Legendre, v = H-coefficients)
  G(y) = sum_n c_n p_n(y), p_0 = 2, p_1 = y, p_n = y (p_{n-1} - p_{n-2})       (symmetrisation x + x/(x-1) = y; RESULT.md Sec. 5)
G-family (CDT Lemma 14.3.1): G, G', G'', G''', int yG, int G, int (G - T_0)/y, ..., int (G - T_3)/y^4 (primitives vanishing at 0).
Pure functions (CDT (10.1.2), (10.1.3), (10.2.1)):
  B1 = 1; B2 = sum_{n>=2} 2 (n-2)! n!/(2n)! y^n; B3 = sum_{n>=1} (n-1)!^2/(2n)! y^n;
  B4 = 4 sum_{n>=0} y^{n+1}/16^n sum_k C(2k,k) C(2n-2k,n-k)/((2k-1)(2n-2k+1)^2);  B5 = sum_{n>=1} (n-1)!^2/((2n-1)(2n-1)!) y^n;
  B6 = int B3/y;  B7 = int B4/y.
Relation search: sum_i P_i(y) f_i(y) = 0 mod y^R with deg P_i <= D.  All series are p-integral up to the truncation (p > 2*Nser + 8),
so nullity 0 mod p implies: no relation over Q with polynomial coefficients of degree <= D [exact, for the (pairs, r) tested].
Usage: python xrank.py p Nser D seed spec   where spec = "a,b:r0,ra,rb,rab;a,b:..." ; r entries may be 'R' (random mod p)."""
import sys, random
import numpy as np
from flint import nmod_mat

def series_tools(p, Nser):
    inv = lambda x: pow(x % p, p - 2, p)
    NN = 2 * Nser + 4
    def uv(a):
        u = [1, a % p]; v = [0, 1]
        for n in range(1, NN):
            i = inv(n + 1)
            u.append(((a * (2 * n + 1) * u[n] - n * u[n - 1]) * i) % p)
            v.append(((a * (2 * n + 1) * v[n] - n * v[n - 1]) * i) % p)
        return u, v
    # p_n(y) truncated at y^Nser, as numpy rows (built once)
    P = np.zeros((NN, Nser), dtype=np.int64)
    P[0, 0] = 2; P[1, 1] = 1
    for n in range(2, NN):
        d = (P[n - 1] - P[n - 2]) % p
        P[n, 1:] = d[:-1]
    def Gser(a, b, r):
        ua, va = uv(a); ub, vb = uv(b)
        r0, ra, rb, rab = r
        c = np.array([(r0 * ua[n] * ub[n] + ra * va[n] * ub[n] + rb * ua[n] * vb[n] + rab * va[n] * vb[n]) % p for n in range(NN)],
                     dtype=object)
        G = [0] * Nser
        # G = c . P  (mod p), done in python ints column by column to avoid overflow
        Pt = P.T
        for j in range(Nser):
            col = Pt[j]
            nz = np.nonzero(col)[0]
            G[j] = int(sum(int(c[n]) * int(col[n]) for n in nz) % p)
        return G
    def deriv(s): return [((n + 1) * s[n + 1]) % p for n in range(len(s) - 1)] + [0]
    def int_ymul(s, j):
        out = [0] * len(s)
        for n in range(j + 1, len(s)):
            out[n] = (s[n - j - 1] * inv(n)) % p
        return out
    def int_neg(s, k):
        out = [0] * len(s)
        for n in range(k, len(s)):
            out[n - k + 1] = (s[n] * inv(n - k + 1)) % p
        return out
    def family(G):
        G1 = deriv(G); G2 = deriv(G1); G3 = deriv(G2)
        return [G, G1, G2, G3, int_ymul(G, 1), int_ymul(G, 0), int_neg(G, 1), int_neg(G, 2), int_neg(G, 3), int_neg(G, 4)]
    def pure():
        fact = [1] * (2 * Nser + 2)
        for i in range(1, len(fact)):
            fact[i] = fact[i - 1] * i % p
        ifact = [inv(f) for f in fact]
        C = lambda n, k: fact[n] * ifact[k] % p * ifact[n - k] % p
        B1 = [1] + [0] * (Nser - 1)
        B2 = [0] * Nser; B3 = [0] * Nser; B4 = [0] * Nser; B5 = [0] * Nser
        for n in range(2, Nser):
            B2[n] = 2 * fact[n - 2] * fact[n] % p * ifact[2 * n] % p
        for n in range(1, Nser):
            B3[n] = fact[n - 1] ** 2 % p * ifact[2 * n] % p
            B5[n] = fact[n - 1] ** 2 % p * inv(2 * n - 1) % p * ifact[2 * n - 1] % p
        i16 = inv(16)
        for n in range(0, Nser - 1):
            s = 0
            for k in range(n + 1):
                s += C(2 * k, k) * C(2 * n - 2 * k, n - k) % p * inv((2 * k - 1) * (2 * n - 2 * k + 1) ** 2)
            B4[n + 1] = 4 * pow(i16, n, p) * (s % p) % p
        B6 = [0] * Nser; B7 = [0] * Nser
        for n in range(1, Nser):
            B6[n] = B3[n] * inv(n) % p
            B7[n] = B4[n] * inv(n) % p
        return [B1, B2, B3, B4, B5, B6, B7]
    return Gser, family, pure

def nullity(funcs, D, rows, p):
    m = len(funcs); cols = m * (D + 1)
    ent = []
    for row in range(rows):
        for s in funcs:
            for j in range(D + 1):
                ent.append(s[row - j] if row - j >= 0 else 0)
    M = nmod_mat(rows, cols, ent, p)
    return cols - M.rank(), cols

def parse(spec, p, rng):
    out = []
    for item in spec.split(';'):
        ab, rs = item.split(':')
        a, b = [int(t) for t in ab.split(',')]
        r = [rng.randrange(1, p) if t == 'R' else int(t) % p for t in rs.split(',')]
        out.append((a, b, r))
    return out

if __name__ == '__main__':
    p = int(sys.argv[1]); Nser = int(sys.argv[2]); D = int(sys.argv[3]); seed = int(sys.argv[4]); spec = sys.argv[5]
    assert p > 2 * Nser + 8
    rng = random.Random(seed)
    Gser, family, pure = series_tools(p, Nser)
    pairs = parse(spec, p, rng)
    B = pure()
    funcs = list(B)
    for a, b, r in pairs:
        funcs += family(Gser(a, b, r))
    rows = Nser - 8
    nl, cols = nullity(funcs, D, rows, p)
    print(f'p={p} Nser={Nser} D={D} pairs={[(a, b) for a, b, _ in pairs]} r={[r if any(x > 10**6 for x in r) is False else "rand" for a, b, r in pairs]} '
          f'm={len(funcs)} cols={cols} rows={rows} nullity={nl}', flush=True)
