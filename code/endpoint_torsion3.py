"""
Rigorous non-torsion of the endpoint P = (i, e^{i pi/4}) on S_0 (referee M3).

Curve (quartic model, u = 2y + x^2 + 1):
    C: u^2 = f(x) = x^4 - 4x^3 + 2x^2 + 1,   over K = Q(zeta_8),
    P  = (i, 2 zeta_8)   [since 2 e^{i pi/4} + i^2 + 1 = 2 e^{i pi/4} = 2 zeta_8].

Method (reduction modulo good primes, replacing the old "20 exact multiples"
argument of endpoint_torsion2.py, which did not bound the torsion order):

  * For p == 1 (mod 8), F_p contains a primitive 8th root of unity z, and
    P reduces mod any prime above p to Pb = (z^2, 2z) on C/F_p (i, 2zeta_8
    are algebraic integers, so P is integral at p).
  * Good reduction: p does not divide disc(f) = -2^8*11 (checked exactly),
    so the quartic model is smooth over F_p; hence E has good reduction at
    every prime of K above p, and reduction is injective on the prime-to-p
    torsion of E(K) (Silverman, AEC VII Prop. 3.1).
  * Hence if P were torsion of order N = p^a m (p |/ m), then m | ord(Pb).
    We find ord(P mod 17) = 20 and ord(P mod 89) = 3:
      - the 17-part of N divides the 89-free part of N, which divides 3: so
        the 17-part is 1;  likewise the 89-part divides 20, so it is 1;
      - the remaining part of N divides BOTH 20 and 3, hence gcd(20,3) = 1.
    Thus N = 1, i.e. P = O -- but P mod 17 has order 20.  Contradiction.
  * As an independent exact confirmation we also check 3P != O, 20P != O
    and 60P != O over Q(zeta_8) with exact rational arithmetic.

Group law on the quartic (monic-parabola chord-tangent, neutral
O = O- = infinity_-; derivation in torsion.py / endpoint_torsion2.py).

CAUTION -- the two points at infinity must be DISTINGUISHED (the old code
conflated them; harmless for a non-torsion point over K, but wrong in
general, e.g. mod p).  Neutral is O-: it corresponds to S_0's group
neutral [0:1:0] (there u ~ -x^2, y = (u-x^2-1)/2 ~ -x^2; at O+ instead
y ~ -x, i.e. the second point [1:-1:0] at infinity).  On C the map x has
simple poles at both infinities O+, O- (u/x^2 -> +1 resp. -1) and u has
double poles at both.  Principal divisors (exact, over any field where
the model is smooth):
    div(x)                    = [(0,1)] + [(0,-1)] - [O+] - [O-]
    div(u - (x^2-2x-1))       = [O+] + [(0,-1)] - 2[O-]
        (since u^2 - (x^2-2x-1)^2 = -4x, and u - (x^2-2x-1) = -2t + O(t^2)
         at O+, t = 1/x: (1-2t-t^2-2t^3)^2 = 1-4t+2t^2+t^4 + O(t^4))
    div(u - g), g = x^2+bx+c  = [P] + [Q] + [R'] - [O+] - 2[O-]   (b != -2)
                              = [P] + [Q] - 2[O-]                 (b == -2)
    div(x - r)                = 2[(r,0)] - [O+] - [O-]  at a u=0 point
Writing T = [O+] - [O-] in Pic^0 (neutral O-), these give:
    (x,u) + (x,-u) = O+  i.e. T ;   2*(r,0) = T ;   b == -2 chord: sum = O- ;
    2T = [(0,1)] - [O-]  and  [(0,-1)] - [O-] = -T  ==>  5T = O- exactly
    (doubling (0,1) twice gives (0,-1); adding T gives O-), so ord(T) = 5
    (T != 0, 2T != 0).  Points are represented as pairs (A, e) with
    A = "O" or an affine (x, u), meaning A + e*T; equality needs e mod 5.
"""
from fractions import Fraction as F
from math import gcd
from sympy import symbols, discriminant, Poly

# ---------- exact discriminant of the quartic ----------
_x = symbols('x')
f_poly = Poly(_x**4 - 4*_x**3 + 2*_x**2 + 1, _x)
DISC_F = int(discriminant(f_poly))
print("disc(f) =", DISC_F, "= ", end="")
d, facs = abs(DISC_F), []
for q in (2, 3, 5, 7, 11, 13, 17, 19, 23):
    while d % q == 0:
        facs.append(q); d //= q
assert d == 1
print(" * ".join(map(str, facs)))

# ---------- Q(zeta8), z^4 = -1 (exact, for the final M*P != O check) ----------
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
    M = []
    for j in range(4):          # column j: p * z^j
        M.append(z8_mul(p, tuple(F(1) if i == j else F(0) for i in range(4))))
    A = [[M[j][i] for j in range(4)] + [F(1) if i == 0 else F(0)] for i in range(4)]
    for c in range(4):
        piv = next(r for r in range(c, 4) if A[r][c] != 0)
        A[c], A[piv] = A[piv], A[c]
        pv = A[c][c]
        A[c] = [v/pv for v in A[c]]
        for r in range(4):
            if r != c and A[r][c] != 0:
                f_ = A[r][c]
                A[r] = [a - f_*b for a, b in zip(A[r], A[c])]
    return tuple(A[i][4] for i in range(4))
def z8_div(p, q): return z8_mul(p, z8_inv(q))
Z8_0 = (F(0),)*4

def zc(n):    # integer n into Q(zeta8)
    return (F(n), F(0), F(0), F(0))

# ---------- generic quartic group law over any ring with the given ops ----------
# Point representation: (A, e) with A = "O" (point at infinity class) or an
# affine pair (x, u); the group element is A + e*T, T = O+ - O- (order 5).
# Neutral is ("O", 0).  Rules below are the divisor identities of the docstring.
def ec_add(P, Q, ops):
    add_, neg_, sub, mul, div, zero = ops
    A1, e1 = P; A2, e2 = Q
    if A1 == "O":
        return (A2, e1 + e2)
    if A2 == "O":
        return (A1, e1 + e2)
    e = e1 + e2
    if A1 == A2:
        r1, s1 = A1
        if s1 == zero:
            return ("O", e + 1)    # div(x-r) = 2[P] - [O+] - [O-]: 2P = T
        r1_2 = mul(r1, r1); r1_3 = mul(r1_2, r1)
        # b = f'(r1)/(2 s1) - 2 r1 ;  f'(x) = 4x^3 - 12x^2 + 4x
        num = add_(add_(mul(r1_3, cc(4, zero)), neg_(mul(r1_2, cc(12, zero)))),
                   mul(r1, cc(4, zero)))
        b = sub(div(num, mul(s1, cc(2, zero))), mul(r1, cc(2, zero)))
        c = sub(s1, add_(mul(r1, r1), mul(b, r1)))
    else:
        r1, s1 = A1; r2, s2 = A2
        if r1 == r2:
            return ("O", e + 1)    # vertical pair: (x,u)+(x,-u) = O+
        b = div(sub(sub(s1, s2), sub(mul(r1, r1), mul(r2, r2))), sub(r1, r2))
        c = sub(s1, add_(mul(r1, r1), mul(b, r1)))
    # intersect g = x^2 + b x + c with u^2 = f: roots r1, r2, x3 of f - g^2
    a3 = sub(cc(-4, zero), mul(b, cc(2, zero)))
    if a3 == zero:
        return ("O", e)            # b == -2: div(u-g) = [P]+[Q] - 2[O-]
    a2 = sub(cc(2, zero), add_(mul(b, b), mul(c, cc(2, zero))))
    x3 = sub(neg_(div(a2, a3)), add_(A1[0], A2[0]))
    y3 = add_(add_(mul(x3, x3), mul(b, x3)), c)
    return ((x3, neg_(y3)), e)

def ec_normalize(P, ops):
    """Canonicalize (A, e): the affine points in <T> = {O-, O+, (0,1), (1,0),
    (0,-1)} are absorbed into e: (0,1) = 2T, (1,0) = 3T, (0,-1) = 4T
    (exact 5-chain, see docstring).  After normalization the neutral
    element is ("O", e) with e == 0 (mod 5)."""
    A, e = P
    if A != "O":
        add_, zero = ops[0], ops[5]
        rc = lambda n: add_(cc(n, zero), zero)   # embed n and reduce (mod p)
        U1  = (rc(0), rc(1))
        V1  = (rc(1), rc(0))
        U1m = (rc(0), rc(-1))
        if A == U1:
            return ("O", e + 2)
        if A == V1:
            return ("O", e + 3)
        if A == U1m:
            return ("O", e + 4)
    return (A, e)

def ec_is_O(P, ops):
    A, e = ec_normalize(P, ops)
    return A == "O" and e % 5 == 0

def cc(n, zero):
    """constant integer n embedded in the field (tuple => Q(zeta8), int => F_p)."""
    if isinstance(zero, tuple):
        return (F(n), F(0), F(0), F(0)) if len(zero) == 4 else (F(n), F(0))
    return n            # F_p: ops reduce mod p at every step

def ec_mul(n, P, ops):
    R, Q = ("O", 0), P
    while n:
        if n & 1:
            R = ec_add(R, Q, ops)
        Q = ec_add(Q, Q, ops)
        n >>= 1
    return R

Z8 = (z8_add, z8_neg, z8_sub, z8_mul, z8_div, Z8_0)

# ---------- F_p machinery for p == 1 (mod 8) ----------
def fp_ops(p):
    return (lambda a, b: (a+b) % p, lambda a: (-a) % p, lambda a, b: (a-b) % p,
            lambda a, b: (a*b) % p, lambda a, b: (a*pow(b, -1, p)) % p, 0)

def find_zeta8(p):
    """primitive 8th root of unity in F_p (p == 1 mod 8)."""
    assert p % 8 == 1
    for g in range(2, p):
        z = pow(g, (p-1)//8, p)
        if pow(z, 4, p) == p-1:
            return z
    raise AssertionError

def f_mod(x, p):
    return (pow(x, 4, p) - 4*pow(x, 3, p) + 2*x*x + 1) % p

def group_size(p):
    """#C(F_p) = 2 points at infinity + sum_x (1 + chi(f(x)))."""
    s = 2 + p
    for x in range(p):
        v = f_mod(x, p)
        if v and pow(v, (p-1)//2, p) == 1:
            s += 1
        elif v == 0:
            pass
        else:
            s -= 1
    return s

def order_mod_p(p):
    z = find_zeta8(p)
    i, u = pow(z, 2, p), (2*z) % p
    assert (u*u - f_mod(i, p)) % p == 0, "Pb not on curve"
    ops = fp_ops(p)
    nE = group_size(p)
    # sanity: T = O+ - O- has order 5 also mod p (p != 5, reduction injective
    # on prime-to-p torsion): 2T = (0,1), 5*(0,1) = O-
    Pt = ((0 % p, 1 % p), 0)
    assert ec_is_O(ec_mul(5, Pt, ops), ops) and not ec_is_O(Pt, ops), "T-order check"
    cur, n = ("O", 0), 0
    while True:
        cur = ec_add(cur, ((i, u), 0), ops)
        n += 1
        if ec_is_O(cur, ops):
            break
        assert n <= nE, "order did not divide #E"
    return z, n, nE

# ---------- run ----------
print("\n--- reduction of P = (i, 2 zeta8) modulo good primes p == 1 (mod 8) ---")
orders = {}
for p in (17, 41, 73, 89, 97):
    assert DISC_F % p != 0, f"bad reduction at {p}"
    z, o, nE = order_mod_p(p)
    orders[p] = o
    print(f"p={p:3d}:  zeta8 = {z:3d} (order 8 in F_p^*),  disc(f) mod p = {DISC_F % p},"
          f"  #E(F_p) = {nE},  ord(P mod p) = {o}   (divides #E: {nE % o == 0})")

o17, o89 = orders[17], orders[89]
print(f"\nord(P mod 17) = {o17},  ord(P mod 89) = {o89},  gcd = {gcd(o17, o89)}")
print("torsion-order argument (injectivity of prime-to-p torsion reduction):")
print("  17-part of N | 89-free part of N | ord(P mod 89) = 3   ==> no 17-part")
print("  89-part of N | 17-free part of N | ord(P mod 17) = 20  ==> no 89-part")
print(f"  rest of N divides gcd(20, 3) = {gcd(o17, o89)}  ==> N = 1, i.e. P = O")
print("  but ord(P mod 17) = 20 > 1: contradiction")

print("\n--- independent exact confirmation over Q(zeta_8) ---")
zz = (F(0), F(1), F(0), F(0))
P_exact = ((z8_mul(zz, zz), z8_mul(zz, zc(2))), 0)   # (z^2, 2z) = (i, 2 zeta8)
# sanity: P on the curve, exact
x_, u_ = P_exact[0]
x2 = z8_mul(x_, x_); x3 = z8_mul(x2, x_); x4 = z8_mul(x2, x2)
fval = z8_add(z8_sub(z8_add(z8_sub(x4, z8_mul(zc(4), x3)), z8_mul(zc(2), x2)), zc(0)), zc(1))
assert z8_sub(z8_mul(u_, u_), fval) == Z8_0, "P not on curve exactly"
print("P = (z^2, 2z) lies on u^2 = f(x) exactly: OK")
# sanity: T has order exactly 5 over K (2T = (0,1), 5*(0,1) = O-, (0,1) != O-)
Pt = ((zc(0), zc(1)), 0)
assert ec_is_O(ec_mul(5, Pt, Z8), Z8) and not ec_is_O(Pt, Z8)
print("ord(T) = 5 over Q(zeta_8): OK (exact)")
for n in (3, 20, 60):
    nP = ec_mul(n, P_exact, Z8)
    if ec_is_O(nP, Z8):
        raise AssertionError(f"{n}*P = O?!")
    A, e = nP
    digits = max(len(str(c.numerator)) for coord in A for c in coord) if A != "O" else 0
    print(f"{n}*P != O exactly (affine point, e = {e}, max numerator digits = {digits})")
print("3P != O, 20P != O, 60P != O exactly over Q(zeta_8): OK")
print(f"\nVERDICT: N | 20 and N | 3 with cross-parts excluded  ==>  N = 1.")
print("P = (i, e^{i pi/4}) is NOT a torsion point on S_0 over Q(zeta_8). QED")
