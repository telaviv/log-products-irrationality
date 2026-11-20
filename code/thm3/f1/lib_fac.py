"""Load the irreducible factors written by f1_factors.py, keyed by name."""
import os
from lib_f1 import CTX, A_, B_, Y_
HERE = os.path.dirname(os.path.abspath(__file__))
def load():
    F = {}
    for line in open(os.path.join(HERE, 'runs', 'factors.txt')):
        if line.startswith('#'): continue
        e, da, td, p = [s.strip() for s in line.split('|')]
        P = eval(p.replace('^', '**'), {'a': A_, 'b': B_, 'y': Y_})
        da, td = int(da), int(td)
        name = {14: 'Phi14', 26: 'Phi26', 79: 'Phi79', 4: 'Phi4'}.get(da) if td > 2 else p
        F[name] = P
    return F
def ev(P, a, b):
    from flint import fmpq
    s = fmpq(0)
    for (ea, eb, ey), c in zip(P.monoms(), P.coeffs()):
        s += fmpq(c) * fmpq(a) ** ea * fmpq(b) ** eb
    return s
