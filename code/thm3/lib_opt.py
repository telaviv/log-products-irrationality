"""lib_opt.py -- maximise M = L - T/m - tau over psi = z exp(g), subject to admissibility (winding 0 about every obstruction,
max|psi| < QLIM), with a soft log-barrier on the distance to obstructions.  [Numerical]"""
import numpy as np, math, json, time
from scipy.optimize import minimize
from lib_map import *

TAU = 1032659 / 242760
MFUN = 17

def barrier_and_grad(c, E, N, p=8):
    z = np.exp(2j * np.pi * (np.arange(N) + 0.37) / N)
    K = len(c)
    ps = Map(c).psi(z)
    V = np.vander(z, K, increasing=True)
    val = 0.0; gr = np.zeros(K)
    for e in E:
        dv = ps - e
        ad = np.abs(dv)
        u = ad ** (-p)
        S = u.mean()
        d = S ** (-1.0 / p)
        val += math.log(d)
        # d log d / dc_n = (1/d) dd/dc_n ; dd/dc_n = S^{-1/p-1} mean(|dv|^{-p-2} Re(conj(dv) psi z^n))
        wgt = ad ** (-p - 2) * np.conj(dv) * ps
        dd = S ** (-1.0 / p - 1) * np.real(wgt @ V) / N
        gr += dd / d
    return val, gr

class Problem:
    def __init__(self, E, N=1024, mu=1e-3, qlim=0.95, m=MFUN, tau=TAU, Nw=16384, dsafe=1e-4):
        self.E = np.asarray(E, dtype=complex); self.N = N; self.mu = mu; self.qlim = qlim
        self.m = m; self.tau = tau; self.Nw = Nw; self.dsafe = dsafe
        self.best = None; self.nev = 0
        self.Eb = self.E[np.abs(self.E) < qlim]
    def admissible(self, c):
        rad = float(np.max(np.abs(Map(c).psi(np.exp(2j * np.pi * np.arange(4096) / 4096)))))
        if not np.isfinite(rad) or rad >= self.qlim:
            return False, rad
        near = self.E[np.abs(self.E) < rad + 0.02]
        if len(near) == 0:
            return True, rad
        wn, rad2 = winding(c, near, self.Nw)
        return bool(np.all(wn == 0)), rad2
    def f(self, c):
        self.nev += 1
        ok, rad = self.admissible(c)
        if not ok:
            return 1e3, np.zeros_like(c)
        T, gT = T_and_grad(c, self.N)
        M = LOG256 + c[0] - T / self.m - self.tau
        gM = -gT / self.m; gM[0] += 1.0
        near = self.Eb
        bv, bg = barrier_and_grad(c, near, self.N) if len(near) else (0.0, np.zeros_like(c))
        obj = M + self.mu * bv
        if self.best is None or M > self.best[0]:
            self.best = (M, c.copy(), T)
        return -obj, -(gM + self.mu * bg)

def optimise(c0, E, N=1024, mus=(1e-3, 3e-4, 1e-4, 3e-5, 1e-5), iters=200, log=None, K=None, restarts=4):
    c = np.array(c0, dtype=float)
    if K is not None and K > len(c):
        c = np.concatenate([c, np.zeros(K - len(c))])
    P = Problem(E, N)
    ok, _ = P.admissible(c)
    assert ok, 'start not admissible'
    for mu in mus:
        P.mu = mu
        for rs in range(restarts):
            r = minimize(P.f, c, jac=True, method='L-BFGS-B', options=dict(maxiter=iters, maxcor=30))
            cnew = r.x if r.fun < 1e2 else P.best[1]
            moved = np.max(np.abs(cnew - c))
            c = cnew
            if moved < 1e-7:
                break
        T, _ = T_and_grad(c, N, False)
        M = LOG256 + c[0] - T / MFUN - TAU
        msg = f'mu={mu:.0e} L={LOG256 + c[0]:.4f} T={T:.4f} M={M:+.4f} nev={P.nev} best={P.best[0]:+.4f}'
        print(msg, flush=True)
        if log: log.write(msg + '\n'); log.flush()
    return P.best[1], P

def ascend(c0, E, N=1024, mu=1e-3, iters=400, eta=1e-3, prec=1.0, log=None, qlim=0.95, every=50, m=MFUN, tau=TAU):
    """adaptive preconditioned gradient ascent on M + mu*barrier with step rejection (admissibility + monotone objective)."""
    P = Problem(E, N, mu=mu, qlim=qlim, m=m, tau=tau)
    c = np.array(c0, dtype=float)
    f, g = P.f(c)
    assert f < 1e2, 'start not admissible'
    K = len(c)
    D = 1.0 / (1.0 + np.arange(K)) ** prec
    for it in range(iters):
        d = -g * D
        while True:
            cn = c + eta * d
            fn, gn = P.f(cn)
            if fn < f:
                c, f, g = cn, fn, gn
                eta *= 1.3
                break
            eta *= 0.4
            if eta < 1e-12:
                break
        if eta < 1e-12:
            break
        if it % every == 0 or it == iters - 1:
            T, _ = T_and_grad(c, N, False)
            msg = f'  it={it} eta={eta:.2e} L={LOG256 + c[0]:.4f} T={T:.4f} M={LOG256 + c[0] - T / m - tau:+.4f} obj={-f:+.5f}'
            print(msg, flush=True)
            if log: log.write(msg + '\n'); log.flush()
    return c, P

def lbfgs(P, c0, iters=500, mem=30, log=None, every=25, tag=''):
    """minimise P.f (returns (f, g); f >= 1e2 means inadmissible) by L-BFGS with backtracking that also rejects inadmissible."""
    c = np.array(c0, dtype=float)
    f, g = P.f(c)
    assert f < 1e2
    S, Y = [], []
    t0 = time.time()
    for it in range(iters):
        q = g.copy(); al = []
        for s, y in reversed(list(zip(S, Y))):
            a = s @ q / (y @ s); al.append(a); q -= a * y
        if S:
            q *= (S[-1] @ Y[-1]) / (Y[-1] @ Y[-1])
        else:
            q *= 1e-3 / max(1e-12, np.linalg.norm(q))
        for (s, y), a in zip(zip(S, Y), reversed(al)):
            b = y @ q / (y @ s); q += s * (a - b)
        d = -q
        if d @ g >= 0:
            d = -g * 1e-3 / max(1e-12, np.linalg.norm(g)); S, Y = [], []
        st = 1.0; ok = False
        for _ in range(40):
            cn = c + st * d
            fn, gn = P.f(cn)
            if fn < 1e2 and fn <= f + 1e-4 * st * (d @ g):
                ok = True; break
            st *= 0.5
        if not ok:
            if S:
                S, Y = [], []; continue
            break
        s = cn - c; y = gn - g
        if s @ y > 1e-14:
            S.append(s); Y.append(y)
            if len(S) > mem: S.pop(0); Y.pop(0)
        c, f, g = cn, fn, gn
        if it % every == 0:
            T = T_and_grad(c, P.N, False)[0]
            msg = f'{tag} it={it} L={LOG256 + c[0]:.5f} T={T:.5f} M={LOG256 + c[0] - T / P.m - P.tau:+.5f} obj={-f:+.5f} st={st:.1e} {time.time()-t0:.0f}s'
            print(msg, flush=True)
            if log: log.write(msg + '\n'); log.flush()
    return c
