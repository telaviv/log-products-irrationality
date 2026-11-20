"""Common code for the F1 check (CDT Lemma 14.3.1 determinant).

R10, R4: transcribed from arXiv:2408.15403v2 p.180-181 (files r10.txt, r4.txt; transcription by the gate referee,
checked against CDT's printed R10(0) = 4a^2b^2 and R10(4)).
A(a,x) = 1/sqrt(1-2ax+x^2) = sum P_n(a) x^n (Legendre), P_A = A(a,.) * A(b,.) (Hadamard), and
G_A(y) = P_A(x) + P_A(x/(x-1)),  y = x^2/(x-1).  Since x + w = x w = y for w = x/(x-1),
x^n + w^n = p_n(y) with p_0 = 2, p_1 = y, p_n = y (p_{n-1} - p_{n-2}).
"""
import re, os
from flint import fmpq, fmpz, fmpq_poly, fmpz_mpoly_ctx

HERE = os.path.dirname(os.path.abspath(__file__))
CTX = fmpz_mpoly_ctx.get(('a', 'b', 'y'), 'lex')
A_, B_, Y_ = CTX.gens()


def parse(fn):
    s = open(os.path.join(HERE, fn), encoding='utf-8', errors='replace').read().replace('\n', ' ')
    s = s.replace('−', '-')
    s = re.sub(r'\s+', ' ', s).strip().rstrip('.,')
    toks = re.findall(r'([+-]?)\s*(\d*)((?:[aby]\d*)*)', s)
    P = CTX.from_dict({})
    for sg, c, mon in toks:
        if c == '' and mon == '':
            continue
        coef = int(c) if c else 1
        if sg == '-':
            coef = -coef
        t = CTX.from_dict({(0, 0, 0): coef})
        for v, e in re.findall(r'([aby])(\d*)', mon):
            t *= {'a': A_, 'b': B_, 'y': Y_}[v] ** (int(e) if e else 1)
        P += t
    return P


R10 = parse('r10.txt')
R4 = parse('r4.txt')


def spec_y(P, a, b):
    """fmpz_mpoly in (a,b,y) -> fmpq_poly in y at numeric a, b (fmpq)."""
    coeffs = {}
    for (ea, eb, ey), c in zip(P.monoms(), P.coeffs()):
        coeffs[ey] = coeffs.get(ey, fmpq(0)) + fmpq(c) * fmpq(a) ** ea * fmpq(b) ** eb
    d = max(coeffs) if coeffs else 0
    return fmpq_poly([coeffs.get(i, fmpq(0)) for i in range(d + 1)])


def legendre(a, N):
    P = [fmpq(1), fmpq(a)]
    for n in range(1, N):
        P.append((fmpq(2 * n + 1) * a * P[n] - n * P[n - 1]) / (n + 1))
    return P[:N + 1]


def GA_series(a, b, K):
    """Coefficients g_0..g_{K-1} of G_A(y) (exact)."""
    N = 2 * K + 2
    Pa, Pb = legendre(fmpq(a), N), legendre(fmpq(b), N)
    # p_n(y) truncated mod y^K as lists
    pm2 = [fmpq(2)] + [fmpq(0)] * (K - 1)          # p_0
    pm1 = [fmpq(0), fmpq(1)] + [fmpq(0)] * (K - 2)  # p_1
    g = [fmpq(0)] * K
    for n in range(N + 1):
        if n == 0:
            cur = pm2
        elif n == 1:
            cur = pm1
        else:
            d = [pm1[i] - pm2[i] for i in range(K)]
            cur = [fmpq(0)] + d[:K - 1]
            pm2, pm1 = pm1, cur
        c = Pa[n] * Pb[n]
        for i in range(K):
            if cur[i] != 0:
                g[i] += c * cur[i]
    return g
