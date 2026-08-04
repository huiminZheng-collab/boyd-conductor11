"""
Endpoint torsion check, pure-Python exact arithmetic (no sympy).

Quartic u^2 = f(x) = x^4 - 4x^3 + 2x^2 + 1, group law via monic parabolas,
neutral O = infinity_- (see torsion.py for the derivation).

Fields:
  QS2  = Q(sqrt2):  (a, b) = a + b*sqrt2
  QZ8  = Q(zeta8):  (a,b,c,d) = a + b*z + c*z^2 + d*z^3,  z^4 = -1
         (contains i = z^2 and sqrt2 = z + z^7 = z - z^3... check: z+z^7 = z - z^3
          since z^7 = z^4*z^3 = -z^3; and (z - z^3)^2 = z^2 - 2z^4 + z^6 = z^2 + 2 - z^2*... )
Points:
  P_pi   = (-1, 2*sqrt2)              in Q(sqrt2)
  P_pi/2 = (i, sqrt2*(1+i)) = (z^2, 2*z)  in Q(zeta8)
"""
from fractions import Fraction as F

# ---------- Q(sqrt2) ----------
def s2_add(p, q): return (p[0]+q[0], p[1]+q[1])
def s2_neg(p): return (-p[0], -p[1])
def s2_sub(p, q): return (p[0]-q[0], p[1]-q[1])
def s2_mul(p, q): return (p[0]*q[0]+2*p[1]*q[1], p[0]*q[1]+p[1]*q[0])
def s2_inv(p):
    d = p[0]**2 - 2*p[1]**2
    return (p[0]/d, -p[1]/d)
def s2_div(p, q): return s2_mul(p, s2_inv(q))
S2_0 = (F(0), F(0)); S2_1 = (F(1), F(0))

# ---------- Q(zeta8), z^4 = -1 ----------
def z8_add(p, q): return tuple(p[i]+q[i] for i in range(4))
def z8_neg(p): return tuple(-v for v in p)
def z8_sub(p, q): return tuple(p[i]-q[i] for i in range(4))
def z8_mul(p, q):
    r = [F(0)]*7
    for i in range(4):
        for j in range(4):
            r[i+j] += p[i]*q[j]
    for k in range(6, 3, -1):
        r[k-4] -= r[k]   # z^k = -z^{k-4}
    return tuple(r[:4])
def z8_inv(p):
    # solve p * v = 1 by 4x4 linear system over Q
    M = []
    for j in range(4):          # column j: p * z^j
        col = z8_mul(p, tuple(F(1) if i==j else F(0) for i in range(4)))
        M.append(col)
    # M v = (1,0,0,0): rows
    A = [[M[j][i] for j in range(4)] + [F(1) if i==0 else F(0)] for i in range(4)]
    n = 4
    for c in range(n):
        piv = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[piv] = A[piv], A[c]
        pv = A[c][c]
        A[c] = [v/pv for v in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f_ = A[r][c]
                A[r] = [a - f_*b for a, b in zip(A[r], A[c])]
    return tuple(A[i][4] for i in range(4))
def z8_div(p, q): return z8_mul(p, z8_inv(q))
Z8_0 = (F(0),)*4

# ---------- generic group law ----------
def make_ops(add_, neg_, sub, mul, div, zero):
    return add_, neg_, sub, mul, div, zero

def ec_add(P, Q, ops):
    add_, neg_, sub, mul, div, zero = ops
    if P == Q:
        r1, s1 = P
        # b = f'(r1)/(2 s1) - 2 r1 ; f'(x) = 4x^3 - 12x^2 + 4x
        r1_2 = mul(r1, r1); r1_3 = mul(r1_2, r1)
        num = sub(sub(mul(r1_3, _c(4, zero)), mul(r1_2, _c(12, zero))), mul(r1, _c(-4, zero)) if False else mul(r1_3, _c(0, zero)))
        # do it plainly:
        num = add_(add_(mul(r1_3, _c(4, zero)), neg_(mul(r1_2, _c(12, zero)))), mul(r1, _c(4, zero)))
        b = sub(div(num, mul(s1, _c(2, zero))), mul(r1, _c(2, zero)))
        c = sub(s1, add_(mul(r1, r1), mul(b, r1)))
    else:
        r1, s1 = P; r2, s2 = Q
        b = div(sub(sub(s1, s2), sub(mul(r1, r1), mul(r2, r2))), sub(r1, r2))
        c = sub(s1, add_(mul(r1, r1), mul(b, r1)))
    # f - g^2, g = x^2 + b x + c : coeffs of x^3, x^2:
    # g^2 = x^4 + 2b x^3 + (b^2+2c) x^2 + 2bc x + c^2
    a3 = sub(_c(-4, zero), mul(b, _c(2, zero)))
    if a3 == zero:
        return "O"   # third intersection at infinity => sum = O
    a2 = sub(_c(2, zero), add_(mul(b, b), mul(c, _c(2, zero))))
    x3 = sub(neg_(div(a2, a3)), add_(P[0], Q[0]))
    y3 = add_(add_(mul(x3, x3), mul(b, x3)), c)
    return (x3, neg_(y3))

def _c(n, zero):
    """constant integer n embedded in the field (detect size)."""
    if len(zero) == 2:
        return (F(n), F(0))
    return (F(n), F(0), F(0), F(0))

def check(P, ops, maxn, name):
    print(f"{name} = {P}")
    cur = P
    for n in range(2, maxn+1):
        cur = ec_add(cur, P, ops)
        if cur == "O":
            print(f"  {n}P = O  ==> ORDER DIVIDES {n}")
            return n
        print(f"  {n}P = {cur}")
    print(f"  not torsion up to {maxn}")
    return None

S2 = make_ops(s2_add, s2_neg, s2_sub, s2_mul, s2_div, S2_0)
Z8 = make_ops(z8_add, z8_neg, z8_sub, z8_mul, z8_div, Z8_0)

P_pi = ((F(-1), F(0)), (F(0), F(2)))          # (-1, 2 sqrt2)
check(P_pi, S2, 20, "P_pi")

z = (F(0), F(1), F(0), F(0))
P_half = (z8_mul(z, z), z8_mul(z, (F(2), F(0), F(0), F(0))))   # (z^2, 2z) = (i, 2 zeta8)
check(P_half, Z8, 20, "P_{pi/2}")
