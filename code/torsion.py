"""
NOTE (legacy): this script conflates the two points at infinity of the
quartic model (multiples 3A/5A print nan/zoo); the conclusion 5A = O via
the 4A = -A route remains valid.  For the rigorous verification see
verify_divisors_k0.py and endpoint_torsion3.py.

Torsion check (proof-strategy support).

Curve S_0: y^2+(x^2+1)y+x^3 = 0. With u = 2y + x^2 + 1 it becomes the quartic
    u^2 = f(x) = x^4 - 4x^3 + 2x^2 + 1.
Points of interest: div(x) = [(0,0)] + [(0,-1)] - 2[P_inf]  on S_0,
i.e. quartic points  A = (x,u) = (0, 1)  (y=0)  and  (0,-1) (y=-1).

Claim needed for the modular-unit / BMZ proof route: A is 5-torsion (E(Q)=Z/5Z).
We verify EXACTLY (rational arithmetic) that 4A = -A, i.e. 5A = O.

Group law on u^2 = quartic (monic), neutral O = infinity_- (point at infinity
with u/x^2 -> -1), via monic parabolas u = g(x) = x^2+bx+c:
  * f - g^2 is a cubic with zeros at the x-coords of P, Q, and a third point P3;
    P+Q = iota(P3) = (x3, -g(x3)).
  * doubling: tangent parabola at P.
  * negation: -R = (r', g(r')) where g = x^2 - 2x + c passes through R and the
    third intersection is at infinity (x^3 coeff of f-g^2 vanishes).
"""
from sympy import Rational, symbols, expand, simplify

x = symbols('x')
f = x**4 - 4*x**3 + 2*x**2 + 1

def add(P, Q):
    """P+Q on u^2=f, neutral infinity_-. P=(r,s) affine."""
    r1, s1 = P; r2, s2 = Q
    if P == Q:  # tangent parabola
        b = (4*r1**3 - 12*r1**2 + 4*r1) / (2*s1) - 2*r1
        c = s1 - r1**2 - b*r1
    else:
        b = (s1 - s2 - (r1**2 - r2**2)) / (r1 - r2)
        c = s1 - r1**2 - b*r1
    g = x**2 + b*x + c
    h = expand(f - g**2)          # cubic with roots r1, r2, x3
    a3 = h.coeff(x, 3); a2 = h.coeff(x, 2)
    x3 = -a2/a3 - r1 - r2         # Vieta
    y3 = x3**2 + b*x3 + c
    return (simplify(x3), simplify(-y3))

def neg(R):
    r, s = R
    c = s - r**2 + 2*r            # g = x^2 - 2x + c through R, third point at infinity
    # f - g^2 = -(2+2c) x^2 + 4c x - c^2 ; roots r, r'
    rp = 2*c/(1+c) - r
    gp = rp**2 - 2*rp + c
    return (simplify(rp), simplify(gp))

A = (Rational(0), Rational(1))
print("A       =", A)
negA = neg(A)
print("-A      =", negA, "  (on curve:", simplify(negA[1]**2 - f.subs(x, negA[0])) == 0, ")")
A2 = add(A, A)
print("2A      =", A2)
A3 = add(A2, A)
print("3A      =", A3)
A4 = add(A2, A2)
print("4A      =", A4)
print("4A == -A ?", A4 == negA, " => 5A = O" if A4 == negA else "")
A5 = add(A4, A)
print("5A via 4A+A =", A5)
