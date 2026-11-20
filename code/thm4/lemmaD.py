"""lemmaD.py -- Lemma D (pilot_F2): the genuine singular sets of inequivalent pairs are disjoint.
R4(a,b;y) = CDT (14.2.1)'s quartic whose roots are the y = x^2/(x-1) images of x = (a +- sqrt(a^2-1))(b +- sqrt(b^2-1)).
Checks (exact, sympy): (1) the transcription: R4(y(x)) * (x-1)^4 vanishes on the 4 products (as a polynomial identity modulo
x^2 - ... : verified by resultant and by 60-digit evaluation); (2) R4 is symmetric, even under (a,b) -> (-a,-b), constant term 1,
leading coefficient 4(a-b)^2; (3) its coefficients in sigma = a^2+b^2, pi = ab; (4) pi and sigma are rational functions of the
y^3 and y^4 coefficients: pi = -(c3 + 2 c4)/(4 c4), sigma = c4/4 + 2 pi; (5) R4(4) is not 0 at integer points of interest
(argument: 4 = y(2), and 2 is not an algebraic unit, whereas every product is a unit) -- here only a numeric spot check."""
import sympy as sp
import mpmath as mp
a, b, y, x = sp.symbols('a b y x')
R4 = (1 + 4*y - 8*a**2*y - 4*a*b*y - 8*b**2*y + 16*a**2*b**2*y + 4*y**2 - 12*a**2*y**2 + 16*a**4*y**2
      + 20*a*b*y**2 - 16*a**3*b*y**2 - 12*b**2*y**2 - 16*a*b**3*y**2 + 16*b**4*y**2 - 8*a**2*y**3 + 16*a*b*y**3
      - 16*a**3*b*y**3 - 8*b**2*y**3 + 32*a**2*b**2*y**3 - 16*a*b**3*y**3 + 4*a**2*y**4 - 8*a*b*y**4 + 4*b**2*y**4)
P = sp.Poly(R4, y)
c = [sp.expand(P.coeff_monomial(y**k)) for k in range(5)]
print('coefficients c0..c4:'); [print(' c%d =' % k, ck) for k, ck in enumerate(c)]
assert c[0] == 1 and sp.factor(c[4]) == 4*(a - b)**2
assert sp.expand(R4 - R4.subs({a: b, b: a}, simultaneous=True)) == 0
assert sp.expand(R4 - R4.subs({a: -a, b: -b}, simultaneous=True)) == 0
s, p_ = sp.symbols('sigma pi')
claimed = [1, 16*p_**2 - 8*s - 4*p_ + 4, 16*(s**2 - 2*p_**2) - 16*p_*s - 12*s + 20*p_ + 4, -16*p_*s + 32*p_**2 - 8*s + 16*p_, 4*s - 8*p_]
sub = {s: a**2 + b**2, p_: a*b}
for k in range(5):
    assert sp.expand(sp.sympify(claimed[k]).subs(sub) - c[k]) == 0, k
cs = claimed
print('in (sigma, pi):'); [print(' c%d =' % k, ck) for k, ck in enumerate(cs)]
assert sp.expand(cs[3] - (-2*cs[4] - 4*p_*cs[4])) == 0
print('identity c3 = -2 c4 - 4 pi c4: OK  =>  pi = -(c3 + 2 c4)/(4 c4), sigma = c4/4 + 2 pi (c4 = 4(a-b)^2 != 0)')
# (1) transcription check: numerically, roots of R4 are y(x) for the four products
mp.mp.dps = 60
for (A, B) in [(5, 7), (201, 203), (200001, 200005), (-200001, -200003)]:
    u = A + mp.sign(A) * mp.sqrt(A * A - 1); v = B + mp.sign(B) * mp.sqrt(B * B - 1)
    vals = []
    for X in (u * v, u / v, v / u, 1 / (u * v)):
        Y = X * X / (X - 1)
        vals.append(abs(sp.lambdify(y, R4.subs({a: A, b: B}), 'mpmath')(Y)) / max(1, abs(Y)) ** 4)
    R44 = R4.subs({a: A, b: B, y: 4})
    print(f' (a,b)=({A},{B}): max |R4(y(x))|/max(1,|y|)^4 over the 4 products = {mp.nstr(max(vals), 3)};  R4(4) = {R44}')
# exact: R4(y(x)) (x-1)^4 is divisible by the quartic in x whose roots are the products: x^4 - 4ab x^3 + ... is checked via resultant
# Products satisfy  X^2 - 2(ab + sqrt..)X + 1 ... ; we use instead: for x = uv, u+1/u = 2a, v + 1/v = 2b.  Let w = x + 1/x; then
# w = uv + 1/(uv) and w' = u/v + v/u satisfy w + w' = (u+1/u)(v+1/v) = 4ab and w w' = (u^2+u^-2) + (v^2+v^-2) = (4a^2-2)+(4b^2-2).
w = sp.symbols('w')
quad_w = w**2 - 4*a*b*w + (4*a**2 + 4*b**2 - 4)            # roots w, w'
quart_x = sp.expand(x**2 * quad_w.subs(w, x + 1/x))          # roots: the four products
num = sp.expand(sp.together(R4.subs(y, x**2/(x - 1)) * (x - 1)**4))
rem = sp.rem(sp.Poly(num, x), sp.Poly(quart_x, x))
print('exact: (x-1)^4 R4(x^2/(x-1)) mod (the quartic of the 4 products) =', rem.as_expr())
r4at4 = sp.factor(R4.subs(y, 4)); print('R4(4) =', r4at4)
