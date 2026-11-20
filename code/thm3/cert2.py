"""cert2.py map.json eps [delta_max M0 prec] -- CERTIFIED (python-flint Arb + exact fmpq) verification that
psi(z) = z exp(sum_k c_k z^k), c_k := round(2^40 c_k)/2^40 (dyadic; this rounded map is the certified object),
phi = h o psi, h(q) = -256 q prod_{n>=1}(1+q^n)^24 = -256 (eta(2 tau)/eta(tau))^24 (q = e^{2 pi i tau}),
satisfies, for all pairs m, n in Z minus {-1,0} with 0 < |m/n - 1| < eps (CDT Sec. 14.4: genuine singular values real with
|y| >= Y1 = 1/eps - 6 (Lemma S, RESULT.md; CDT's |y| > 1/eps is off by O(1)); alpha real, |alpha| <= Y0 = eps^4/16):
  (a) |psi| < 1 on T;  (b) |phi| > 2 Y0 on T;  (c) phi(T) misses the real rays |y| >= Y1;
and computes a rigorous upper bound R+ for the rearrangement integral R = iint max(g(s), g(t)) ds dt, g(s) = log|phi(e^{2 pi i s})|
(CDT Thm 6.0.2, trivial partition), and L = log|phi'(0)| = log 256 + c_0 exactly.  Certified iff R+ < 17 (L - tau).

Method: dyadic partition of s in [0,1); on each interval I with midpoint s_m and half-width rho (in s):
  log phi(s) = log phi(s_m) + w,  |w| <= delta := 2 pi rho * D,   D >= sup_I |d/dtheta log phi|,
  d/dtheta log phi = i * (q h'/h)(psi) * (1 + z g'(z)),  |q h'/h| = |1 + 24 sum n q^n/(1+q^n)| <= 1 + 24 S(r),
  S(r) = sum_n n r^n/(1 - r^n) (r >= sup_I |psi|; exact upper bound with a geometric tail),
  |1 + z g'(z)| <= |1 + z_m g'(z_m)| + rz * sum k^2 |c_k| (1+rz)^(k-1),  rz = 2 pi rho  (|z - z_m| <= rz),
  sup_I |psi| <= |psi(z_m)| + rz * sup|psi'|, |psi'| = |e^g (1 + z g')| <= exp(Re g(z_m) + rz G1) * sup|1 + z g'|.
(The branch of log phi along I is the continuous one; |w| <= delta follows from integrating the derivative bound.)
Intervals are bisected until (a), (b), (c) are decided and delta <= delta_max.  Output runs/cert2_<tag>.json."""
import sys, json, math, time
from fractions import Fraction
from flint import arb, acb, fmpq, ctx


def main():
    src = sys.argv[1]
    eps = Fraction(sys.argv[2])
    dmax = float(sys.argv[3]) if len(sys.argv) > 3 else 0.03
    M0 = int(sys.argv[4]) if len(sys.argv) > 4 else 2 ** 13
    ctx.prec = int(sys.argv[5]) if len(sys.argv) > 5 else 96
    cf = json.load(open(src))['c']
    cq = [fmpq(int(round(x * 2 ** 40)), 2 ** 40) for x in cf]
    K = len(cq)
    ca = [acb(arb(c)) for c in cq]
    cd = [acb(arb(k * cq[k])) for k in range(1, K)]          # g' coefficients
    pi = arb.pi()
    twopi = 2 * pi
    twopii = acb(0, twopi)
    Y1 = arb(fmpq(eps.denominator, eps.numerator) - 6)            # Lemma S (RESULT.md): genuine |y| >= 1/eps - 6
    Y0 = arb(fmpq(eps.numerator ** 4, 16 * eps.denominator ** 4))
    logY1 = Y1.log()
    log2Y0 = (2 * Y0).log()
    TAU = arb(fmpq(1032659, 242760))
    log256 = arb(256).log()
    L = log256 + arb(cq[0])

    def S_upper(r):
        """rigorous upper bound of sum_{n>=1} n r^n / (1 - r^n) for 0 < r < 1 (r an arb upper bound)."""
        s = arb(0)
        n = 1
        rn = r
        while True:
            term = n * rn / (1 - rn)
            s += term
            if n > 5 and term < arb(1e-12):
                rn1 = rn * r
                # sum_{k>n} k r^k/(1-r^k) <= (1/(1-r^{n+1})) sum_{k>n} k r^k = r^{n+1}((n+1) - n r)/((1-r)^2 (1-r^{n+1}))
                s += rn1 * ((n + 1) - n * r) / ((1 - r) ** 2 * (1 - rn1))
                return s.upper()
            n += 1
            rn = rn * r

    def poly(cs, z):
        v = acb(0)
        for c in reversed(cs):
            v = v * z + c
        return v

    absq = [abs(arb(c)) for c in cq]

    def Gk(p, rz):
        s = arb(0)
        for k in range(1, K):
            s += (k ** p) * absq[k] * (1 + rz) ** (k - 1)
        return s.upper()

    Gcache = {}

    def evaluate(a, b):
        mid = (a + b) / 2
        rho = (b - a) / 2
        key = str(rho)
        rz = (twopi * arb(rho)).upper()
        if key not in Gcache:
            Gcache[key] = (Gk(1, rz), Gk(2, rz))
        G1, G2 = Gcache[key]
        zm = acb(0, twopi * arb(mid)).exp()
        gm = poly(ca, zm)
        gpm = poly(cd, zm)
        psim = zm * gm.exp()
        w1 = (abs(1 + zm * gpm) + rz * G2).upper()               # sup_I |1 + z g'|
        dpsi = ((gm.real + rz * G1).exp() * w1).upper()          # sup_I |psi'|
        r = (abs(psim) + rz * dpsi).upper()
        if not (r < 1):
            return dict(ok=False, why='psi')
        D = (1 + 24 * S_upper(r)) * w1
        delta = (rz * D).upper()
        if float(psim.real.mid()) < 0:
            tau = ((-psim).log() + acb(0, pi)) / twopii
        else:
            tau = psim.log() / twopii
        e1 = tau.modular_eta()
        e2 = (2 * tau).modular_eta()
        lam = log256 + 24 * (abs(e2).log() - abs(e1).log())      # log|phi(s_m)|
        phim = -256 * (e2 / e1) ** 24
        return dict(ok=True, delta=delta, lam=lam, phim=phim, r=r)

    t0 = time.time()
    stack = [(fmpq(k, M0), fmpq(k + 1, M0)) for k in range(M0)]
    vals = []
    nsplit = 0
    maxr = 0.0
    minlam = 1e9
    maxdelta = 0.0
    why = {}
    while stack:
        a, b = stack.pop()
        e = evaluate(a, b)
        ok = e['ok']
        if ok:
            dl = e['delta']
            lam = e['lam']
            ok = bool(dl < dmax)
            if not ok:
                why['delta'] = why.get('delta', 0) + 1
            if ok and not bool(lam - dl > log2Y0):                 # (b)
                ok = False
                why['small'] = why.get('small', 0) + 1
            if ok and not bool(lam + dl < logY1):                  # (c): |phi| may reach Y1: need arg away from 0, pi
                ph = e['phim']
                # phi(s) = phim * exp(w), |w| <= delta: arg moves by <= delta (as delta < 1).  Decide via |Im phim| > |phim| sin(delta)
                sd = arb(dl).sin()
                ok = bool(abs(ph.imag) > abs(ph) * sd)
                if not ok:
                    why['ray'] = why.get('ray', 0) + 1
        else:
            why[e['why']] = why.get(e['why'], 0) + 1
        if not ok:
            if b - a < fmpq(1, 2 ** 34):
                print('FAIL: undecided at s =', float(a), {k: str(v) for k, v in e.items()})
                json.dump(dict(src=src, eps=str(eps), certified=False, fail_at=float(a)), open('runs/cert2_fail.json', 'w'))
                return
            m = (a + b) / 2
            stack.append((a, m))
            stack.append((m, b))
            nsplit += 1
            continue
        up = (e['lam'] + e['delta']).upper()
        v = fmpq(int(math.ceil(float(up) * 2 ** 40)) + 2, 2 ** 40)
        assert arb(v) > up
        vals.append((v, b - a))
        maxr = max(maxr, float(e['r']))
        minlam = min(minlam, float(e['lam']))
        maxdelta = max(maxdelta, float(e['delta']))
        if len(vals) % 50000 == 0:
            print(' ..', len(vals), 'intervals', f'{time.time() - t0:.0f}s', flush=True)
    vals.sort(key=lambda t: t[0])
    R = fmpq(0)
    W = fmpq(0)
    for v, w in vals:
        Wn = W + w
        R += v * (Wn * Wn - W * W)
        W = Wn
    assert W == 1
    Rb = arb(R)
    bound = Rb / (L - TAU)
    cert = bool(bound < 17) and bool(L > TAU) and bool(L > 4)
    res = dict(src=src, eps=str(eps), delta_max=dmax, M0=M0, prec=ctx.prec, intervals=len(vals), splits=nsplit,
               L=L.str(20), R_upper=float(R), bound_m=bound.str(15),
               margin_lower=(L - Rb / 17 - TAU).str(15), certified=cert, sup_abs_psi=maxr, min_log_abs_phi_mid=minlam,
               max_delta=maxdelta, log_2Y0=float(log2Y0), log_Y1=float(logY1), split_reasons=why, secs=time.time() - t0,
               c_dyadic=[str(c) for c in cq])
    print(json.dumps({k: v for k, v in res.items() if k != 'c_dyadic'}, indent=1))
    tag = src.replace(chr(92), '/').split('/')[-1].replace('.json', '') + f'_eps{float(eps):.1e}'
    json.dump(res, open(f'runs/cert2_{tag}.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
