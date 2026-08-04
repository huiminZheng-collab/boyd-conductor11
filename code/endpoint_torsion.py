"""
Endpoint torsion check: the split-path boundary points on the quartic
u^2 = x^4-4x^3+2x^2+1 are
  theta=0:   (x,y)=(1,-1)  -> (x,u)=(1,0)   = -A (already 5-torsion)
  theta=pi:  (x,y)=(-1,-1+sqrt2) -> (x,u)=(-1, 2*sqrt2)
  theta=pi/2: x=i -> u^2 = f(i) = 4i -> u = sqrt(2)*(1+i)
If these are torsion (=> cusps on X_1(11) by Manin-Drinfeld), the BMZ route applies.
We compute multiples over QQ(sqrt2) and QQ(i, sqrt2) with exact arithmetic.
"""
from sympy import Rational, sqrt, I, symbols, expand, simplify, AlgebraicNumber, to_number_field

x = symbols('x')
f = x**4 - 4*x**3 + 2*x**2 + 1

def add(P, Q):
    r1, s1 = P; r2, s2 = Q
    if P == Q:
        b = (4*r1**3 - 12*r1**2 + 4*r1) / (2*s1) - 2*r1
        c = s1 - r1**2 - b*r1
    else:
        b = (s1 - s2 - (r1**2 - r2**2)) / (r1 - r2)
        c = s1 - r1**2 - b*r1
    g = x**2 + b*x + c
    h = expand(f - g**2)
    a3 = h.coeff(x, 3); a2 = h.coeff(x, 2)
    x3 = -a2/a3 - r1 - r2
    y3 = x3**2 + b*x3 + c
    return (simplify(x3), simplify(-y3))

def check_torsion(P, maxn=20, name="P"):
    print(f"checking {name} = {P}")
    mults = {1: P}
    cur = P
    for n in range(2, maxn+1):
        cur = add(cur, P)
        mults[n] = cur
        xs, ys = cur
        # point at infinity shows up as zoo/nan in rational arithmetic
        if xs in (None,) or xs.has(symbols('zoo')) or str(xs) in ('zoo','nan'):
            print(f"  {n}P = (point at infinity)  => ORDER DIVIDES {n}")
            return n
    for n, (xs, ys) in mults.items():
        print(f"  {n}P = ({xs}, {ys})")
    print(f"  no torsion up to {maxn}")
    return None

s2 = sqrt(2)
check_torsion((-1, 2*s2), 20, "P_pi")
check_torsion((I, s2*(1+I)), 20, "P_{pi/2}")
