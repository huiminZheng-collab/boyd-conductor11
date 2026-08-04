"""
The decisive check for the BMZ proof route (cf. Samart 2023 §3 "path becomes closed").

Split-path boundary divisor:  D = 2[P_{pi/2}] - [P_0] - [P_pi]
on E: u^2 = x^4-4x^3+2x^2+1,  P_0=(1,0) is 5-torsion (=-A).
Conjecture (C3) via regulator/BMZ needs  T := 2*P_{pi/2} - P_pi  to be TORSION in Pic^0(E)
(then a multiple of the path is closed between cusps, and BMZ gives the L-value).

Everything computed in Q(zeta8):  i = z^2, sqrt2 = z - z^3  (z^4=-1).
"""
from fractions import Fraction as F

def z8_add(p, q): return tuple(p[i]+q[i] for i in range(4))
def z8_neg(p): return tuple(-v for v in p)
def z8_sub(p, q): return tuple(p[i]-q[i] for i in range(4))
def z8_mul(p, q):
    r = [F(0)]*7
    for i in range(4):
        for j in range(4):
            r[i+j] += p[i]*q[j]
    for k in range(6, 3, -1):
        r[k-4] -= r[k]
    return tuple(r[:4])
def z8_inv(p):
    M = []
    for j in range(4):
        M.append(z8_mul(p, tuple(F(1) if i==j else F(0) for i in range(4))))
    A = [[M[j][i] for j in range(4)] + [F(1) if i==0 else F(0)] for i in range(4)]
    for c in range(4):
        piv = next(r for r in range(c, 4) if A[r][c] != 0)
        A[c], A[piv] = A[piv], A[c]
        pv = A[c][c]; A[c] = [v/pv for v in A[c]]
        for r in range(4):
            if r != c and A[r][c] != 0:
                f_ = A[r][c]
                A[r] = [a - f_*b for a, b in zip(A[r], A[c])]
    return tuple(A[i][4] for i in range(4))
def z8_div(p, q): return z8_mul(p, z8_inv(q))
Z0 = (F(0),)*4
def C(n): return (F(n), F(0), F(0), F(0))

def ec_add(P, Q):
    if P == "O": return Q
    if Q == "O": return P
    if P == Q:
        r1, s1 = P
        r1_2 = z8_mul(r1, r1); r1_3 = z8_mul(r1_2, r1)
        num = z8_add(z8_add(z8_mul(r1_3, C(4)), z8_neg(z8_mul(r1_2, C(12)))), z8_mul(r1, C(4)))
        b = z8_sub(z8_div(num, z8_mul(s1, C(2))), z8_mul(r1, C(2)))
        c = z8_sub(s1, z8_add(z8_mul(r1, r1), z8_mul(b, r1)))
    else:
        r1, s1 = P; r2, s2 = Q
        b = z8_div(z8_sub(z8_sub(s1, s2), z8_sub(z8_mul(r1, r1), z8_mul(r2, r2))), z8_sub(r1, r2))
        c = z8_sub(s1, z8_add(z8_mul(r1, r1), z8_mul(b, r1)))
    a3 = z8_sub(C(-4), z8_mul(b, C(2)))
    if a3 == Z0:
        return "O"
    a2 = z8_sub(C(2), z8_add(z8_mul(b, b), z8_mul(c, C(2))))
    x3 = z8_sub(z8_neg(z8_div(a2, a3)), z8_add(P[0], Q[0]))
    y3 = z8_add(z8_add(z8_mul(x3, x3), z8_mul(b, x3)), c)
    return (x3, z8_neg(y3))

def ec_neg(R):
    r, s = R
    c = z8_sub(s, z8_sub(z8_mul(r, r), z8_mul(r, C(2))))   # c = s - r^2 + 2r
    onepc = z8_add(C(1), c)
    rp = z8_sub(z8_div(z8_mul(c, C(2)), onepc), r)
    gp = z8_add(z8_sub(z8_mul(rp, rp), z8_mul(rp, C(2))), c)
    return (rp, gp)

def ec_sub(P, Q): return ec_add(P, ec_neg(Q))

z = (F(0), F(1), F(0), F(0))
P_half = (z8_mul(z, z), z8_mul(z, C(2)))                       # (i, 2 zeta8)
sqrt2 = (F(0), F(1), F(0), F(-1))                              # z - z^3
P_pi = (C(-1), z8_mul(sqrt2, C(2)))                            # (-1, 2 sqrt2)

# sanity: points on curve? u^2 = f(x)  (spot-check via known construction, skip)

T = ec_sub(ec_add(P_half, P_half), P_pi)
print("T = 2P_{pi/2} - P_pi =", T)
cur = T
for n in range(2, 31):
    cur = ec_add(cur, T)
    if cur == "O":
        print(f"{n}T = O  ==> T IS TORSION, order divides {n}")
        break
    if n <= 6 or n in (10, 15, 20, 25, 30):
        sx, sy = cur
        h = max(max(len(str(v.numerator)), len(str(v.denominator))) for v in sx if v != 0)
        print(f"  {n}T: x-height digits ~ {h}")
else:
    print("T not torsion up to 30")
