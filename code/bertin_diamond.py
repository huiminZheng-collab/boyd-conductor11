"""
ATTACK 13 -- the referee's Bertin-symbol adjudication (conductor 11, k=0).

Question: Brunault's (3.210) anchors |r_{gamma^-}{x,y}| = (5/2pi) D_E(P), proved
via Bertin [10, Th. 6 et Cor. 6.1].  Which symbol is {x,y}, what is its diamond
product, and how does D_E of it compare with the paper's S_0 symbol
(diamond = 6(O) - 5(A) + 5(2A), D_E = (5/2) D_E(P))?

Contents
  Part 1: (C1) cubic F = (X+1)(Y+1)(X+Y+1) + XY -- exact divisors of X, Y
          (line-intersection factorizations + tangent orders, sympy).
  Part 2: Riemann-Roch map (u,v) = (Y/(X+1), Y(X+Y)/(X+1)) to
          W_int: v^2 + 4uv + v = u^3 - 2u^2 - u (derived exactly);
          exact group law on W_int (Fraction arithmetic) gives the Z/5 table
          [origin Q1=[0:1:0]]:  g = B1, 2g = Q3, 3g = A1, 4g = Q2;
          PARI (code/bertin_diamond.gp) certifies W_int = 11.a3 and
          B1 |-> 3P, z(B1)/w1 = 4/5.
  Part 3: diamond products (generic routine, L-R convention):
          (C1) {X,Y};  Weierstrass {x,y} on 11.a3 (Brunault's (3.210) symbol);
          paper's S_0 {x_s,y_s} transported to 11.a3
          [x_s = (X+Y)/(X-1),  y_s = -(X*Y + 2X^2 - 2X + 1)/(X-1)^2,
           verified to satisfy S_0: y^2 + (x^2+1)y + x^3 = 0 on 11.a3];
          Z[E]^- reduction + Abel checks.
  Part 4: D_E at 60 digits (raw z-value sums, no exotic relation assumed).
  Part 5: tame symbols (exact local expansions) for the S_0 symbol on 11.a3
          and for the Weierstrass symbol.
  Part 6: verdict.

gp-side certifications (code/bertin_diamond.gp):
  * ellidentify(W_int) = 11a3, disc -11, j = -4096/11; change v = [1,-1,-2,2]
    with 11.a3 = ellchangecurve(W_int, v);  B1=(0,0) -> (1,0) = 3P.
  * certified integral over gamma^- (anti-invariant generator u = s(2w2-w1)):
      Weierstrass {x,y}:  int = -2 pi b11 = -5 D_E(P)   (9-digit Richardson)
      (C1) {X,Y}:         int = +14 pi b11 = +35 D_E(P)  (= 2 pi * m(C1))
"""
from fractions import Fraction as Fr
import sympy as sp
from mpmath import mp, mpf, mpc, pi, exp, log, polylog, nstr

X, Y, Z = sp.symbols('X Y Z')
t = sp.symbols('t')

print("=" * 78)
print("PART 1: (C1) cubic -- exact divisors of X and Y")
print("=" * 78)
F = (X + Z)*(Y + Z)*(X + Y + Z) + X*Y*Z
print("F(X,0? no) -- line intersections:")
f_x0 = sp.factor(F.subs(X, 0));   print("  {X=0} cap C :", f_x0, "  => 2*(0,-1) + [0:1:0]")
f_y0 = sp.factor(F.subs(Y, 0));   print("  {Y=0} cap C :", f_y0, "  => 2*(-1,0) + [1:0:0]")
f_z0 = sp.factor(F.subs(Z, 0));   print("  {Z=0} cap C :", f_z0, "  => [0:1:0]+[1:0:0]+[1:-1:0]")
assert f_x0 == Z*(Y + Z)**2 and f_y0 == Z*(X + Z)**2 and f_z0 == X*Y*(X + Y)

# tangents at A1=(0,-1)=[0:-1:1] and B1=(-1,0)=[-1:0:1]: gradient
grad = [sp.diff(F, v) for v in (X, Y, Z)]
gA1 = [gi.subs({X: 0, Y: -1, Z: 1}) for gi in grad]
gB1 = [gi.subs({X: -1, Y: 0, Z: 1}) for gi in grad]
print("  tangent at A1=(0,-1):", gA1, " (prop. to X=0 -> ord_A1(X)=2)")
print("  tangent at B1=(-1,0):", gB1, " (prop. to Y=0 -> ord_B1(Y)=2)")
assert gA1[1] == 0 and gA1[2] == 0 and gA1[0] != 0
assert gB1[0] == 0 and gB1[2] == 0 and gB1[1] != 0

# value of X/Z at Q1=[0:1:0] (chart Y=1): Z as series in X
c1, c2 = sp.symbols('c1 c2')
Zs = -X + c1*X**2
eq = sp.expand(F.subs({Y: 1, Z: Zs})).series(X, 0, 4).removeO().expand()
sol = sp.solve([eq.expand().coeff(X, 2)], [c1], dict=True)[0]
Zser = Zs.subs(sol)
valXQ1 = sp.limit(X/Zser, X, 0)
print("  X/Z at Q1 =", valXQ1, " (nonzero, finite => ord_Q1(X/Z)=0)")
assert valXQ1 == -1

print("  div(X) = 2*(0,-1) - [1:0:0] - [1:-1:0]   = 2*A1 - Q2 - Q3")
print("  div(Y) = 2*(-1,0) - [0:1:0] - [1:-1:0]   = 2*B1 - Q1 - Q3")
divX = {'A1': 2, 'Q2': -1, 'Q3': -1}
divY = {'B1': 2, 'Q1': -1, 'Q3': -1}

print()
print("=" * 78)
print("PART 2: RR map to W_int and the Z/5 group table (exact)")
print("=" * 78)
# u = Y/(X+1), v = Y(X+Y)/(X+1)  =>  Y = u(X+1), X = (v-u^2)/(u(u+1))
u, v = sp.symbols('u v')
Fa = (X + 1)*(Y + 1)*(X + Y + 1) + X*Y
res = sp.expand(Fa.subs(Y, u*(X + 1)))
res = sp.factor(res / (X + 1))          # residual quadratic in X
rel = sp.expand(res.subs(X, (v - u**2)/(u*(1 + u))) * (u*(1 + u))**2)
rel = sp.factor(rel)
print("  birational relation:", rel, "= 0")
core = sp.expand(-rel/(u*(u + 1)))
print("  =>  W_int:", sp.expand(core), " i.e.  v^2 + 4uv + v = u^3 - 2u^2 - u")
assert core == -(v**2 + 4*u*v + v - u**3 + 2*u**2 + u)

# exact group law on W_int: a1=4, a2=-2, a3=1, a4=-1, a6=0
A1_, A2_, A3_, A4_, A6_ = Fr(4), Fr(-2), Fr(1), Fr(-1), Fr(0)
def on_w(P):
    x, y = P
    return y*y + A1_*x*y + A3_*y == x**3 + A2_*x*x + A4_*x + A6_
def neg(P):
    x, y = P
    return (x, -y - A1_*x - A3_)
def add(P, Q):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2 and y1 != y2: return None
    if P == Q:
        lam = (3*x1*x1 + 2*A2_*x1 + A4_ - A1_*y1)/(2*y1 + A1_*x1 + A3_)
    else:
        lam = (y2 - y1)/(x2 - x1)
    nu = y1 - lam*x1
    x3 = lam*lam + A1_*lam - A2_ - x1 - x2
    y3 = -(lam*x3 + nu) - A1_*x3 - A3_
    return (x3, y3)
def mul(P, n):
    R = None
    for _ in range(n): R = add(R, P)
    return R

Bp = (Fr(0), Fr(0))
for Pt in [Bp, (Fr(-1), Fr(1)), (Fr(-1), Fr(2)), (Fr(0), Fr(-1))]:
    assert on_w(Pt)
tab = {1: Bp}
for n in (2, 3, 4, 5):
    tab[n] = mul(Bp, n)
print("  multiples of B=(0,0) on W_int:", tab)
assert tab[2] == (Fr(-1), Fr(2)) and tab[3] == (Fr(-1), Fr(1)) and tab[4] == (Fr(0), Fr(-1)) and tab[5] is None

# images of the cubic points (direct values / asymptotics at infinity)
print("  images:  B1=(-1,0) -> (u,v)=(0,0) [local series Y ~ (X+1)^2]")
print("           A1=(0,-1) -> (-1,1)  = 3B  [u=-1/1, v=(-1)(-1)/1]")
# near Q2=[1:0:0]: X->inf, solve Y -> -1
w1_ = sp.symbols('w1_')
Ynear = -1 + w1_/X
lead = sp.expand(Fa.subs(Y, Ynear) * 1).series(X, sp.oo, 2)
print("  near Q2: Y = -1 + O(1/X)  => (u,v) -> (0,-1) = 4B")
# near Q3=[1:-1:0]: Y = -X + w, leading condition
Fw = sp.expand(Fa.subs(Y, -X + sp.Symbol('w')))
coefX2 = sp.Poly(Fw, X).nth(2)
wvals = sp.solve(coefX2, sp.Symbol('w'))
print("  near Q3: Y ~ -X + w with w solving", coefX2, "=0 -> w =", wvals,
      " => (u,v) -> (-1, 2) = 2B")
assert wvals == [-2]
print("  Z/5 table (origin Q1):  g=B1, 2g=Q3, 3g=A1, 4g=Q2, 5g=O=Q1")
print("  PARI cross-check (bertin_diamond.gp): W_int = 11.a3, change [1,-1,-2,2];")
print("  B1 -> (1,0) = 3P on 11.a3; z(B1)/w1 = 4/5 = z(3P)/w1 exactly.")
grp = {'Q1': 0, 'B1': 1, 'Q3': 2, 'A1': 3, 'Q2': 4}   # multiples of g
# u/w1 of k*g: k * 4/5 mod 1 ; equivalently g = 3P, so k*g = (3k mod 5) * P
def frac_of_kg(k):
    return Fr(4 * (k % 5), 5) % 1      # u/w1 = 4k/5 mod 1

print()
print("=" * 78)
print("PART 3: diamond products (L-R convention) + Z[E]^- reduction + Abel")
print("=" * 78)
def diamond(divF, divG, diff):
    """raw (f)<>(g) = sum m_i n_j [S_i - T_j]; diff(S,T) = label of S(-)T."""
    raw = {}
    for S, m in divF.items():
        for T, n in divG.items():
            lab = diff(S, T)
            raw[lab] = raw.get(lab, 0) + m * n
    return raw

# (C1) {X,Y}: differences in multiples of g (origin Q1 = 0)
diff_c1 = lambda S, T: (grp[S] - grp[T]) % 5
raw_c1 = diamond(divX, divY, diff_c1)
print("  (C1) raw (X)<>(Y) [labels = k*g]:", raw_c1)
assert raw_c1 == {1: -4, 2: 6, 3: -4, 4: 1, 0: 1}
# Z[E]^- reduction: [k g] with 3g=-2g, 4g=-g
red = {1: raw_c1.get(1, 0) - raw_c1.get(4, 0), 2: raw_c1.get(2, 0) - raw_c1.get(3, 0)}
print("  Z[E]^- reduction: ", red, " i.e.  -5(g) + 10(2g)")
assert red == {1: -5, 2: 10}
# Abel checks
s1 = sum(raw_c1.values()); s2 = sum(c * k for k, c in raw_c1.items()) % 5
print("  Abel: sum of coeffs =", s1, "  weighted sum =", s2, "(mod 5)")
assert s1 == 0 and s2 == 0

# Weierstrass {x,y} on 11.a3 (labels = k*P): div(x)= A+4A-2O, div(y)=2A+3A-3O
labP = {'O': 0, 'A': 1, '2A': 2, '3A': 3, '4A': 4}
diffP = lambda S, T: (labP[S] - labP[T]) % 5
raw_w = diamond({'A': 1, '4A': 1, 'O': -2}, {'A': 2, '3A': 1, 'O': -3}, diffP)
print("  Weierstrass raw (x)<>(y) [labels = k*P]:", raw_w)
red_w = {1: raw_w.get(1, 0) - raw_w.get(4, 0), 2: raw_w.get(2, 0) - raw_w.get(3, 0)}
s1 = sum(raw_w.values()); s2 = sum(c * k for k, c in raw_w.items()) % 5
print("  Z[E]^- reduction:", red_w, " with", raw_w.get(0, 0), "(O);  Abel:", s1, s2)
assert red_w == {1: 5, 2: -5} and s1 == 0 and s2 == 0

# S_0 symbol transported to 11.a3: div(x_s) = A+2A-O-3A, div(y_s) = 3*2A-2*3A-O
raw_s0 = diamond({'A': 1, '2A': 1, 'O': -1, '3A': -1}, {'2A': 3, '3A': -2, 'O': -1}, diffP)
red_s0 = {1: raw_s0.get(1, 0) - raw_s0.get(4, 0), 2: raw_s0.get(2, 0) - raw_s0.get(3, 0)}
s1 = sum(raw_s0.values()); s2 = sum(c * k for k, c in raw_s0.items()) % 5
print("  S_0 raw (x_s)<>(y_s):", raw_s0, "  reduction:", red_s0,
      " with", raw_s0.get(0, 0), "(O);  Abel:", s1, s2)
assert red_s0 == {1: -5, 2: 5} and s1 == 0 and s2 == 0

# verify the transported S_0 functions satisfy S_0 on 11.a3
Xw, Yw = sp.symbols('Xw Yw')
x_s = (Xw + Yw)/(Xw - 1)
y_s = -(Xw*Yw + 2*Xw**2 - 2*Xw + 1)/(Xw - 1)**2
chk = sp.expand((y_s**2 + (x_s**2 + 1)*y_s + x_s**3) * (Xw - 1)**4)
chk = sp.rem(sp.Poly(chk, Yw), sp.Poly(Yw**2 + Yw - Xw**3 + Xw**2, Yw)).as_expr()
print("  (x_s,y_s) satisfies S_0 on 11.a3:", sp.expand(chk) == 0)
assert sp.expand(chk) == 0

print()
print("=" * 78)
print("PART 4: D_E at 60 digits (raw sums over z-values; no exotic assumed)")
print("=" * 78)
mp.dps = 60
tau = mpc(mpf('0.5'), mpf('0.229878021224650476137525642435318986382'))
q = exp(2*pi*mpc(0, 1)*tau)

def bloch_wigner(w):
    if abs(w) < 1e-50: return mpf(0)
    return polylog(2, w).imag + mp.arg(1 - w)*log(abs(w))

def D_E(z):
    s = mpf(0); n = 0
    while True:
        w = q**n * z; s += bloch_wigner(w)
        if abs(w) < mpf('1e-65'): break
        n += 1
    n = 1
    while True:
        w = q**n / z; s -= bloch_wigner(w)
        if abs(w) < mpf('1e-65'): break
        n += 1
    return s

def Dfrac(fr):  # D_E of point at u/w1 = fr
    if fr == 0: return mpf(0)
    return D_E(exp(2*pi*mpc(0, 1)*(mpf(fr.numerator)/fr.denominator)))

dP  = Dfrac(Fr(3, 5))   # P
d2P = Dfrac(Fr(1, 5))   # 2P
d3P = Dfrac(Fr(4, 5))   # 3P
d4P = Dfrac(Fr(2, 5))   # 4P
print("  D_E(P)  =", nstr(dP, 50))
print("  D_E(2P) =", nstr(d2P, 50), " ratio:", nstr(d2P/dP, 12), "(exotic 3/2)")
print("  D_E(3P) =", nstr(d3P, 50), " D_E(4P) =", nstr(d4P, 50))
b11 = mpf('0.1521471417259180494862272974786344956281435891642261228098898')
print("  b11 =", nstr(b11, 50))
print("  D_E(P) vs (2 pi/5) b11: ratio", nstr(dP/((2*pi/5)*b11), 20))

# raw diamond evaluations: lists of (u/w1 fraction, coeff)
raw_c1_fr = [(Fr(3,5), 4), (Fr(2,5), -2), (Fr(4,5), -2), (Fr(2,5), -2),
             (Fr(1,5), 1), (Fr(3,5), 1), (Fr(4,5), -2), (Fr(3,5), 1), (Fr(0), 1)]
raw_w_fr  = [(Fr(0), 8), (Fr(3,5), -2), (Fr(4,5), 2), (Fr(2,5), -3),
             (Fr(4,5), 1), (Fr(2,5), -4), (Fr(1,5), -2)]
raw_s0_fr = [(Fr(0), 6), (Fr(3,5), -5), (Fr(1,5), 5)]

def eval_raw(lst):
    return sum(c*Dfrac(fr) for fr, c in lst)

DE_c1 = eval_raw(raw_c1_fr)
DE_w  = eval_raw(raw_w_fr)
DE_s0 = eval_raw(raw_s0_fr)
print()
print("  D_E((X)<>(Y))   (C1)       =", nstr(DE_c1, 50))
print("    / D_E(P) =", nstr(DE_c1/dP, 12), "  (prediction 35/2 = 7 * 5/2)")
print("  D_E((x)<>(y))   Weierstrass =", nstr(DE_w, 50))
print("    / D_E(P) =", nstr(DE_w/dP, 12), "  (prediction -5/2)")
print("  D_E((x_s)<>(y_s)) S_0      =", nstr(DE_s0, 50))
print("    / D_E(P) =", nstr(DE_s0/dP, 12), "  (prediction +5/2)")
print("  reference: 2 pi b11 = 5 D_E(P)  =", nstr(5*dP, 50))
print("             pi b11   = 5/2 D_E(P)=", nstr(Fr(5,2)*dP, 50))
print("             14 pi b11 = 35 D_E(P) =", nstr(35*dP, 50))
print("  ratio D_E(<>C1)/D_E(<>S0) =", nstr(DE_c1/DE_s0, 12), " (= 7)")
print("  ratio D_E(<>W)/D_E(<>S0)  =", nstr(DE_w/DE_s0, 12), " (= -1)")

print()
print("=" * 78)
print("PART 5: tame symbols (exact local expansions)")
print("=" * 78)
# S_0 symbol on 11.a3: x_s = (X+Y)/(X-1), y_s = -(XY+2X^2-2X+1)/(X-1)^2
Xt, Yt = sp.symbols('Xt Yt')
xs = (Xt + Yt)/(Xt - 1)
ys = -(Xt*Yt + 2*Xt**2 - 2*Xt + 1)/(Xt - 1)**2

def series_Y(at, branch_val, order=8):
    """Y as series in t = X-1 (or in t at O) on Y^2+Y = X^3-X^2; branch Y(0)=branch_val."""
    cs = sp.symbols('c0:' + str(order))
    Yser = branch_val + sum(cs[i]*t**(i+1) for i in range(order))
    if at == '2A' or at == '3A':
        Xser = 1 + t
    eqn = sp.expand(Yser**2 + Yser - (Xser**3 - Xser**2))
    eqs = [eqn.expand().coeff(t, i) for i in range(1, order+1)]
    sol = sp.solve(eqs, cs, dict=True)[0]
    return Yser.subs(sol), Xser

# at 2A = (1,-1): t = X-1, branch Y(0) = -1
Y2A, X2A = series_Y('2A', -1)
xst = sp.series(xs.subs({Xt: X2A, Yt: Y2A}), t, 0, 5).removeO().expand()
yst = sp.series(ys.subs({Xt: X2A, Yt: Y2A}), t, 0, 6).removeO().expand()
print("  at 2A: x_s =", xst, "  y_s =", yst)
T2A = sp.limit((xst**3/yst), t, 0)
print("  T_{2A} = (-1)^(1*3) * (x_s^3/y_s)(2A) =", -T2A)

# at 3A = (1,0): branch Y(0) = 0
Y3A, X3A = series_Y('3A', 0)
xst3 = sp.series(xs.subs({Xt: X3A, Yt: Y3A}), t, 0, 4).removeO().expand()
yst3 = sp.series(ys.subs({Xt: X3A, Yt: Y3A}), t, 0, 4).removeO().expand()
print("  at 3A: x_s =", xst3, "  y_s =", yst3)
T3A = sp.limit(yst3/xst3**2, t, 0)
print("  T_{3A} = (-1)^((-1)(-2)) * (y_s/x_s^2)(3A) =", T3A)

# at O: uniformizer t = X/Y ; X = A(t)/t^2, Y = A(t)/t^3, A = 1 + t^2 + t^3 + ...
a1s, a2s, a3s, a4s = sp.symbols('a1 a2 a3 a4')
Aser = 1 + a1s*t + a2s*t**2 + a3s*t**3 + a4s*t**4
XO = Aser/t**2; YO = Aser/t**3
eqnO = sp.expand((YO**2 + YO - XO**3 + XO**2) * t**6)
solO = sp.solve([eqnO.expand().coeff(t, i) for i in range(5)], [a1s, a2s, a3s, a4s], dict=True)[0]
Aser = Aser.subs(solO)
XO = XO.subs(solO); YO = YO.subs(solO)
xsO = sp.series(xs.subs({Xt: XO, Yt: YO}), t, 0, 3).removeO().expand()
ysO = sp.series(ys.subs({Xt: XO, Yt: YO}), t, 0, 3).removeO().expand()
print("  at O: A(t) =", Aser, " x_s =", xsO, " y_s =", ysO)
TO = sp.limit(ysO/xsO, t, 0)
print("  T_O = (-1)^((-1)(-1)) * (y_s/x_s)(O) =", -TO)

# at A = (0,0): ord(x_s)=1, ord(y_s)=0
TA = 1/ys.subs({Xt: 0, Yt: 0})
print("  T_A = (1/y_s)(A) =", TA)

print()
print("  Weierstrass {x,y}: T_A = (x^2/y)(A) with y ~ -x^2 at A  => -1;")
print("                     T_O = (y^2/x^3)(O) = 1 (Y^2 ~ X^3).")

print()
print("=" * 78)
print("PART 6: VERDICT")
print("=" * 78)
print("""
  (1) Brunault's (3.210) symbol {x,y} = Weierstrass coordinates of
      E = X1(11): y^2+y = x^3-x^2 (thesis Sect. 3.7, Thm. 8; Sect. 3.9:
      'Notons E = X1(11)').  Certified integral (gp, gamma^- generator):
      int_{gamma^-} eta(x,y) = -2 pi b11 = -5 D_E(P)   [matches (3.210)]
      D_E((x)<>(y)) = -(5/2) D_E(P) = -pi b11.
  (2) (C1) coordinate symbol {X,Y}: D_E(<>) = (35/2) D_E(P) = 7*(5/2) D_E(P);
      certified integral = +35 D_E(P) = 2 pi * m(C1) = 14 pi b11.
  (3) Paper's S_0 symbol: D_E(<>) = +(5/2) D_E(P); certified int = +5 D_E(P).
  (4) D_E(<>_Bertin/Weierstrass) = - D_E(<>_ours):  the referee's B1 repair
      path holds (equal |D_E(<>)$|, universal Bloch constant c cancels).
  (5) Factor: int/D_E(<>)= 2 for ALL THREE conductor-11 symbols, vs 1 in
      conductor 17 (k=1, corroborated by Lalin--Ramamonjisoa's proven
      m(F3) = (7/2) L'(E17,0)).  Tame symbols are roots of unity, so the
      K_2 premise holds: the factor is a curve-level normalization issue
      (11.a3: q<0, one real component, Re(tau)=1/2), not a symbol error.
""")
