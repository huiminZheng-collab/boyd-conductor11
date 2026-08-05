"""
CERTIFICATION of the CORRECTED k=0 (conductor 11) divisors (referee 2nd pass).

Curve S_0: y^2 + (x^2+1) y + x^3 = 0, projective closure
    H: Y^2 Z + (X^2 + Z^2) Y + X^3 = 0.
Points at infinity (Z=0): X^2(X+Y) = 0  =>  O1 = [0:1:0], Qinf = [1:-1:0].

Claimed corrections under test:
  on S_0:  div(x) = [(0,0)] + [(0,-1)] - [O1] - [Qinf]
           div(y) = 3[(0,0)] - 2[O1] - [Qinf]
  on 11.a3 (Y^2+Y = X^3-X^2, A = (0,0), 5A = O):
           div(x) = [A] + [2A] - [O] - [3A]
           div(y) = 3[2A] - 2[3A] - [O]
  diamond: (x)<>(y) = 6(O) - 5(A) + 5(2A)  in Z[E]^-  (odd convention).

All arithmetic EXACT (sympy rationals, truncated power series with rational
coefficients).  Birational map formulas copied from code/kappa_exact.py:
  u = 2y + x^2 + 1,  F = u + x^2 - 2x,  G = x(F-1),
  X = (F+1)/2,  Y = (G-F-1)/2,   inverse: x = (X+Y)/(X-1).
"""
import sys
from fractions import Fraction
import sympy as sp

OUT = open("notes/attack12-divisors.txt", "w")
def log(*a):
    s = " ".join(str(v) for v in a)
    print(s)
    OUT.write(s + "\n")

PASS = []
def verdict(name, ok, detail=""):
    PASS.append(ok)
    log(f"VERDICT [{name}]: {'PASS' if ok else 'FAIL'}  {detail}")

# =====================================================================
log("=" * 72)
log("PART 1: points at infinity of the projective closure of S_0")
log("=" * 72)
X, Y, Z = sp.symbols('X Y Z')
H = Y**2 * Z + (X**2 + Z**2) * Y + X**3
log("homogenization H =", H)
Hz0 = sp.factor(H.subs(Z, 0))
log("H | Z=0  =", Hz0)
O1   = (0, 1, 0)
Qinf = (1, -1, 0)
on_curve = all(sp.expand(H.subs({X: p[0], Y: p[1], Z: p[2]})) == 0
               for p in (O1, Qinf))
# smoothness of the two points (partials w.r.t. chart coordinates)
smooth = True
for (px, py, pz) in (O1, Qinf):
    if py != 0:   # chart Y=1, coords (a=X/Y, b=Z/Y)
        a, b = sp.symbols('a b')
        g = sp.expand(H.subs({X: a, Y: 1, Z: b}))
        pa, pb = px / py, pz / py
        if (sp.diff(g, a).subs({a: pa, b: pb}) == 0 and
                sp.diff(g, b).subs({a: pa, b: pb}) == 0):
            smooth = False
log("both points lie on H:", on_curve, "; both smooth:", smooth)
verdict("P1: exactly two points at infinity",
        Hz0 == X**2 * (X + Y) and on_curve and smooth,
        "H|Z=0 = X^2*(X+Y)  =>  O1=[0:1:0], Qinf=[1:-1:0] only")

# =====================================================================
log("")
log("=" * 72)
log("PART 2: div(x), div(y) on S_0 via exact local expansions")
log("=" * 72)
t = sp.symbols('t')

def trunc_poly(expr, n, var=t):
    """truncate polynomial expr in var to order < n."""
    p = sp.Poly(sp.expand(expr), var)
    return sum(c * var**k for (k,), c in p.terms() if k < n)

def ord_poly(p, var=t):
    """order of vanishing of a (possibly zero) polynomial in var."""
    p = sp.expand(p)
    if p == 0:
        return sp.oo
    return min(k for (k,), c in sp.Poly(p, var).terms() if c != 0)

def ord_rational(num, den):
    return ord_poly(num) - ord_poly(den)

NSER = 14   # truncation order (>> any order we need)

# --- affine zeros of x and y on S_0 ---
xs, ys_ = sp.symbols('x y')
S0 = ys_**2 + (xs**2 + 1) * ys_ + xs**3
zeros_x = sp.solve(S0.subs(xs, 0), ys_)
zeros_y = sp.solve(S0.subs(ys_, 0), xs)
log("  solve: x=0 -> y =", zeros_x, ";  y=0 -> x =", zeros_y)
assert sorted(zeros_x) == [-1, 0] and zeros_y == [0]

# --- (0,0): uniformizer x (dF/dy = 2y+x^2+1 = 1 != 0 there) ---
# y(1 + y + x^2) = -x^3  =>  y = -x^3/(1 + y + x^2),  iterate in Q[[x]]
yser = sp.Integer(0)
for _ in range(NSER):
    yser = sp.expand(sp.series(-xs**3 / (1 + yser + xs**2), xs, 0, NSER).removeO())
ord_y_00 = ord_poly(yser, xs)
log(f"at (0,0):  y = {yser}")
log(f"  ord_(0,0)(x) = 1,  ord_(0,0)(y) = {ord_y_00}   (expect 3)")

# --- (0,-1): t0 = y+1 satisfies t0*(t0 - 1 + x^2) = x^2 - x^3 ---
#     => t0 = (t0^2 - x^2 + x^3)/(1 - x^2),  iterate in Q[[x]]
t0 = sp.Integer(0)
for _ in range(NSER):
    t0 = sp.expand(sp.series((t0**2 - xs**2 + xs**3) / (1 - xs**2), xs, 0, NSER).removeO())
ord_t0 = ord_poly(t0, xs)
log(f"at (0,-1): y+1 = {t0}")
log(f"  ord_(0,-1)(x) = 1,  ord_(0,-1)(y) = 0 (y=-1 there),  ord(y+1) = {ord_t0}")

# --- chart Y=1 at infinity: coords (a,b) = (X/Y, Z/Y); x = a/b, y = 1/b ---
a, b = sp.symbols('a b')
# equation:  b(1+b) = -a^2(1+a)
# O1: (a,b) = (0,0).  b = -a^2(1+a) - b^2, iterate.
bser = sp.Integer(0)
for _ in range(NSER):
    bser = trunc_poly(-a**2 * (1 + a) - bser**2, NSER, a)
ord_b_O1 = ord_poly(bser, a)
log(f"at O1:  b = {bser}   (ord = {ord_b_O1})")
log(f"  ord_O1(x) = ord(a) - ord(b) = 1 - {ord_b_O1} = {1 - ord_b_O1}  (expect -1)")
log(f"  ord_O1(y) = -ord(b) = {-ord_b_O1}  (expect -2)")

# Qinf: (a,b) = (-1,0).  c = a+1;  b(1+b) = -c(c-1)^2;  b = -c(c-1)^2 - b^2.
c = sp.symbols('c')
bser2 = sp.Integer(0)
for _ in range(NSER):
    bser2 = trunc_poly(-c * (c - 1)**2 - bser2**2, NSER, c)
ord_b_Q = ord_poly(bser2, c)
log(f"at Qinf (c = a+1):  b = {bser2}   (ord = {ord_b_Q})")
log(f"  ord_Qinf(x) = ord(c-1) - ord(b) = 0 - {ord_b_Q} = {-ord_b_Q}  (expect -1)")
log(f"  ord_Qinf(y) = -ord(b) = {-ord_b_Q}  (expect -1)")

ok2x = (1 - ord_b_O1 == -1) and (-ord_b_Q == -1)
ok2y = (ord_y_00 == 3) and (-ord_b_O1 == -2) and (-ord_b_Q == -1)
degx = 1 + 1 + (1 - ord_b_O1) + (-ord_b_Q)
degy = ord_y_00 + (-ord_b_O1) + (-ord_b_Q)
log(f"degree check: deg div(x) = {degx}, deg div(y) = {degy}  (both must be 0)")
verdict("P2a: div(x) = [(0,0)]+[(0,-1)]-[O1]-[Qinf]", ok2x and degx == 0)
verdict("P2b: div(y) = 3[(0,0)]-2[O1]-[Qinf]", ok2y and degy == 0)

# =====================================================================
log("")
log("=" * 72)
log("PART 3: 11.a3 side (Y^2+Y = X^3-X^2), Abel + transported orders")
log("=" * 72)
# --- exact group law, a1=0, a2=-1, a3=1, a4=a6=0 ---
def eneg(P):
    if P is None:
        return None
    return (P[0], -P[1] - 1)

def eadd(P, Q):
    if P is None:
        return Q
    if Q is None:
        return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2 and y1 != y2:
        return None
    if P == Q:
        lam = (3*x1**2 - 2*x1) / (2*y1 + 1)
        nu = (-x1**3 - y1) / (2*y1 + 1)
    else:
        lam = (y2 - y1) / (x2 - x1)
        nu = (y1*x2 - y2*x1) / (x2 - x1)
    x3 = lam**2 + 1 - x1 - x2
    y3 = -lam * x3 - nu - 1
    return (Fraction(x3), Fraction(y3))

def emul(n, P):
    R, Q = None, P
    while n:
        if n & 1:
            R = eadd(R, Q)
        Q = eadd(Q, Q)
        n >>= 1
    return R

A = (Fraction(0), Fraction(0))
A2 = emul(2, A); A3 = emul(3, A); A4 = emul(4, A); A5 = emul(5, A)
log(f"A = (0,0); 2A = {A2}; 3A = {A3}; 4A = {A4}; 5A = {A5}")
ok_grp = (A2 == (1, -1) and A3 == (1, 0) and A4 == (0, -1) and A5 is None)
verdict("P3a: A generates Z/5 on 11.a3", ok_grp,
        "2A=(1,-1), 3A=(1,0), 4A=(0,-1), 5A=O")

# --- transported function y_S on 11.a3 (birational map from kappa_exact.py) ---
Xw, Yw = sp.symbols('Xw Yw')
xw = (Xw + Yw) / (Xw - 1)                    # inverse map x(X,Y)
yS = Xw - 1 + xw - xw**2                     # y_S = (u - x^2 - 1)/2, u = F - x^2 + 2x, F = 2X-1
log("transported y_S = X - 1 + x - x^2,  x = (X+Y)/(X-1)")

# certify transport: substitute forward map (x,y) -> (X,Y) and reduce mod S_0
xx, yy = sp.symbols('xx yy')
uu = 2*yy + xx**2 + 1
FF = uu + xx**2 - 2*xx
GG = xx * (FF - 1)
Xf = (FF + 1) / 2
Yf = (GG - FF - 1) / 2
xwf = (Xf + Yf) / (Xf - 1)
ySf = sp.together(Xf - 1 + xwf - xwf**2)
num = sp.together(ySf - yy).as_numer_denom()[0]
rem = sp.rem(sp.Poly(sp.expand(num), yy), sp.Poly(yy**2 + (xx**2 + 1)*yy + xx**3, yy))
transport_ok = sp.expand(rem.as_expr()) == 0
log("y_S(X(x,y),Y(x,y)) - y  reduces to 0 mod S_0:", transport_ok)

# values at A and 4A (must be -1, i.e. y does NOT vanish there)
val_A  = sp.Rational(0 - 1 + (0 + 0) / (0 - 1) - ((0 + 0) / (0 - 1))**2)
val_4A = sp.Rational(0 - 1 + (0 - 1) / (0 - 1) - ((0 - 1) / (0 - 1))**2)
log(f"y_S(A) = {val_A},  y_S(4A) = {val_4A}   (expect both -1)")

# --- local expansions on 11.a3 (uniformizer t; ord of rational exprs) ---
def ord_expr(expr):
    n_, d_ = sp.fraction(sp.together(expr))
    return ord_poly(n_) - ord_poly(d_)

# at 2A = (1,-1): uniformizer t = X-1 (2Y+1 = -1 != 0).  Y = -1 + s, s^2 - s = t(1+t)^2
s = sp.Integer(0)
for _ in range(NSER):
    s = trunc_poly(s**2 - t * (1 + t)**2, NSER)
Xl = 1 + t; Yl = -1 + s
xl_2A = sp.together((Xl + Yl) / (Xl - 1))
ySl_2A = sp.together(Xl - 1 + xl_2A - xl_2A**2)
ord_x_2A  = ord_expr((Xl + Yl) / (Xl - 1))
ord_yS_2A = ord_expr(ySl_2A)
log(f"at 2A: x = {sp.series(xl_2A, t, 0, 5).removeO()} ...  ord = {ord_x_2A} (expect 1)")
log(f"       y_S = {sp.series(ySl_2A, t, 0, 6).removeO()} ...  ord = {ord_yS_2A} (expect 3)")

# at 3A = (1,0): uniformizer t = X-1 (2Y+1 = 1 != 0).  Y -> 0:  Y = t(1+t)^2 - Y^2
r = sp.Integer(0)
for _ in range(NSER):
    r = trunc_poly(t * (1 + t)**2 - r**2, NSER)
Xl = 1 + t; Yl = r
xl_3A = sp.together((Xl + Yl) / (Xl - 1))
ySl_3A = sp.together(Xl - 1 + xl_3A - xl_3A**2)
ord_x_3A  = ord_expr((Xl + Yl) / (Xl - 1))
ord_yS_3A = ord_expr(ySl_3A)
log(f"at 3A: ord(x) = {ord_x_3A} (expect -1);  ord(y_S) = {ord_yS_3A} (expect -2)")

# at O: Tate uniformizer t = -X/Y, w = -1/Y:  X = t/w, Y = -1/w,
#   w = t^3 - t^2 w + w^2
w = sp.Integer(0)
for _ in range(NSER):
    w = trunc_poly(t**3 - t**2 * w + w**2, NSER)
Xl = sp.together(t / w); Yl = sp.together(-1 / w)
xl_O = sp.together((Xl + Yl) / (Xl - 1))
ySl_O = sp.together(Xl - 1 + xl_O - xl_O**2)
ord_x_O  = ord_expr(xl_O)
ord_yS_O = ord_expr(ySl_O)
log(f"at O:  ord(X) = {ord_expr(t / w)} (expect -2);  ord(x) = {ord_x_O} (expect -1);"
    f"  ord(y_S) = {ord_yS_O} (expect -1)")

# at A = (0,0): uniformizer X (2Y+1 = 1).  Y = X^3 - X^2 - Y^2
rA = sp.Integer(0)
for _ in range(NSER):
    rA = trunc_poly(t**3 - t**2 - rA**2, NSER)
Xl = t; Yl = rA
xl_A = sp.together((Xl + Yl) / (Xl - 1))
ord_x_A = ord_expr(xl_A)
log(f"at A:  ord(x) = {ord_x_A} (expect 1)")

ok_yS = (val_A == -1 and val_4A == -1 and ord_yS_2A == 3
         and ord_yS_3A == -2 and ord_yS_O == -1 and transport_ok)
verdict("P3b: div(y) = 3[2A]-2[3A]-[O] on 11.a3", ok_yS,
        f"triple zero at 2A, double pole 3A, simple pole O, y(A)=y(4A)=-1")
ok_xE = (ord_x_A == 1 and ord_x_2A == 1 and ord_x_3A == -1 and ord_x_O == -1)
verdict("P3c: div(x) = [A]+[2A]-[O]-[3A] on 11.a3", ok_xE)

# --- Abel's criterion (group law, exact) ---
S_ = emul(3, A2)          # 3*(2A)
T_ = emul(2, A3)          # 2*(3A)
abel_y = eadd(S_, eneg(T_)) is None
Sx = eadd(A, A2)          # A + 2A
abel_x = eadd(Sx, eneg(A3)) is None
log(f"Abel y: 3*(2A) - 2*(3A) = {S_} - {T_} = O ?  {abel_y}")
log(f"Abel x: A + 2A - 3A = {Sx} - {A3} = O ?  {abel_x}")
verdict("P3d: Abel criterion for both divisors", abel_y and abel_x)

# bonus: image of O1 = [0:1:0] under the birational map (exact series at O1).
# At O1: u/x^2 -> -1, i.e. u = -x^2 s(1/x), s(t) = sqrt(1 - 4t + 2t^2 + t^4), t = 1/x.
u_at_O1 = sp.sqrt(1 - 4*t + 2*t**2 + t**4).series(t, 0, 8).removeO()
u_neg = sp.expand(-u_at_O1 / t**2)          # u = -x^2 s(1/x), x = 1/t
F_at_O1 = sp.expand(u_neg + 1/t**2 - 2/t)   # F = u + x^2 - 2x  (poles cancel exactly)
G_at_O1 = sp.expand((F_at_O1 - 1) / t)      # G = x(F-1) = (F-1)/t
X_at_O1 = sp.expand((F_at_O1 + 1) / 2)
Y_at_O1 = sp.expand((G_at_O1 - F_at_O1 - 1) / 2)
def const_term(expr):
    return sp.Poly(sp.expand(expr), t).coeff_monomial(t**0)
img_O1 = (sp.Rational(const_term(X_at_O1)), sp.Rational(const_term(Y_at_O1)))
log(f"bonus: image of O1 under birational map = {img_O1}  (expect 3A = (1,0))")
bonus_ok = (img_O1 == (1, 0))
verdict("P3e: O1 maps to 3A = (1,0) under the birational map", bonus_ok)

# =====================================================================
log("")
log("=" * 72)
log("PART 4: corrected diamond product (x)<>(y) in Z[E]^-, 5A = O")
log("=" * 72)
from collections import defaultdict
divx = {1: 1, 2: 1, 0: -1, 3: -1}     # [A]+[2A]-[O]-[3A]
divy = {2: 3, 3: -2, 0: -1}           # 3[2A]-2[3A]-[O]   (CORRECTED)
raw = defaultdict(int)
for P, m in divx.items():
    for Q, n in divy.items():
        raw[(P - Q) % 5] += m * n
log("raw diamond (coeff of (P-Q)):", dict(sorted(raw.items())))
# odd convention: (3A) = -(2A), (4A) = -(A)
red = defaultdict(int)
for P, v in raw.items():
    if P == 0:
        red[0] += v
    elif P in (1, 2):
        red[P] += v
    else:
        red[5 - P] -= v
red = dict(sorted(red.items()))
log("reduced in Z[E]^- (odd convention):", red)
expected = {0: 6, 1: -5, 2: 5}
ok_diamond = red == expected
verdict("P4: (x)<>(y) = 6(O) - 5(A) + 5(2A)", ok_diamond)

# numeric cross-check of D_E value (not part of the exact certificate)
DEP0 = 0.1911937370843316957549544343121738161012   # D_E(A)
val = (-5 + 5 * 1.5) * DEP0                          # exotic D_E(2A) = 3/2 D_E(A)
pib11 = 3.1415926535897932385 * 0.15214714172591805  # pi * b_11
log(f"numeric cross-check: D_E(diamond) = (5/2) D_E(A) = {val:.16f}")
log(f"                     pi*b_11                  = {pib11:.16f}"
    f"   rel diff = {abs(val - pib11) / pib11:.2e}")

# =====================================================================
log("")
log("=" * 72)
log(f"OVERALL: {'ALL CHECKS PASS' if all(PASS) else '*** SOME CHECKS FAILED ***'}")
log("=" * 72)
OUT.close()
sys.exit(0 if all(PASS) else 1)
