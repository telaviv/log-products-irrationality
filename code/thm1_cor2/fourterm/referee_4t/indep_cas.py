import sympy as sp
from math import comb
z, u = sp.symbols('z u')
def Qn(nn):
    return sum((-1)**(nn + l) * comb(nn + l, l) * sum(comb(nn, l - k) * comb(2*nn + k, 2*k) * comb(nn + k, k) * comb(4*nn + 2*k, 2*nn + k)
               for k in range(max(0, l - nn), min(2*nn, l) + 1)) * z**(2*l) for l in range(3*nn + 1))
fs = [sp.log(1 + u), sp.log(1 - u), sp.log(1 + u) * sp.log(1 - u)]
ser = [sp.series(f, u, 0, 24).removeO() for f in fs]
def polypart(expr):
    pz = sp.Poly(sp.expand(expr * z**40), z)
    return sp.expand(sum(c * z**(e - 40) for (e,), c in pz.terms() if e >= 40))
rows = []
for nn in range(4):
    Q = sp.expand(Qn(nn))
    rows.append([Q] + [polypart(Q * s.subs(u, 1 / z)) for s in ser])
C0 = sp.factor(sp.Matrix(rows).det())
print('Casoratian C(0)(z) =', C0)
print('C(0)(q) != 0 for q=2..60:', all(C0.subs(z, q) != 0 for q in range(2, 61)))
print('rational roots:', sp.Poly(sp.numer(sp.together(C0)), z).ground_roots())
