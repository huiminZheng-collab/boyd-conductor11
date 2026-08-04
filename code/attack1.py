"""
Boyd conductor-11 attack — Step 2: Mahler measures and the split integral.

Targets
-------
1. [proven, Brunault 2005/06]  m((1+x)(1+y)(1+x+y)+xy) = 7 b_11
2. [proven, Brunault 2006]     m(y^2+(x^2+2x-1)y+x^3) = 5 b_11
3. [OPEN, Boyd 1998 (2-33), Samart 2023 eq. (4.1)]  S_0 = y^2+(x^2+1)y+x^3:
       I_split := (1/pi)∫_0^{pi/2} log|y_-(e^{iθ})| dθ
                - (1/pi)∫_{pi/2}^{pi} log|y_-(e^{iθ})| dθ   =?  -b_11
   where y_-(x) = -(x^2+1)(1/2 - sqrt(1/4 - x^3/(x^2+1)^2))  (principal sqrt).
4. m(S_0) itself: Boyd says "seemingly not r*b_11" — reproduce and PSLQ.

Method: reduce to 1D by Jensen's formula. For P = A(x) y^2 + B(x) y + C(x),
   m(P) = (1/2π) ∫_0^{2π} [ log|A(e^{iθ})| + log+|y_+| + log+|y_-| ] dθ.
"""
from mpmath import mp, mpf, mpc, pi, log, sqrt, exp, fabs, nstr, pslq
import sys, time

sys.path.insert(0, ".")
from b11 import b11

mp.dps = 80

def roots_y(A, B, C, x):
    """Roots in y of A(x)y^2+B(x)y+C(x) at x (|x|=1)."""
    A_, B_, C_ = A(x), B(x), C(x)
    disc = sqrt(B_**2 - 4*A_*C_)
    return (-B_ + disc)/(2*A_), (-B_ - disc)/(2*A_)

def logp(v):
    """log^+|v| = max(log|v|,0)."""
    lv = fabs(v)
    return log(lv) if lv > 1 else mpf(0)

def mahler_quady(A, B, C, verbose=False):
    f = lambda th: (log(fabs(A(exp(2*mpc(0,1)*pi*th))))
                    + sum(logp(y) for y in roots_y(A, B, C, exp(2*mpc(0,1)*pi*th))))
    return mp.quad(f, [0, mpf('0.25'), mpf('0.5'), mpf('0.75'), 1])

# --- polynomials (as coefficient functions of x) ---
one = lambda x: mpc(1)
# P_B1 = (1+x)(1+y)(1+x+y)+xy = (1+x) y^2 + (x^2+4x+2) y + (1+x)^2
A1 = lambda x: 1+x;  B1 = lambda x: x**2+4*x+2;  C1 = lambda x: (1+x)**2
# P_B2 = y^2+(x^2+2x-1)y+x^3
A2 = one;    B2 = lambda x: x**2+2*x-1;  C2 = lambda x: x**3
# S_0  = y^2+(x^2+1)y+x^3
A3 = one;    B3 = lambda x: x**2+1;      C3 = lambda x: x**3

def yminus_S0(x):
    return -(x**2+1)*(mpf('0.5') - sqrt(mpf('0.25') - x**3/(x**2+1)**2))

def main():
    t0 = time.time()
    b = b11(80)
    print("b_11      =", nstr(b, 60)); print()

    print("[1] m((1+x)(1+y)(1+x+y)+xy) vs 7 b_11  [Brunault, proven]")
    m1 = mahler_quady(A1, B1, C1)
    print("  m   =", nstr(m1, 60)); print("  7b  =", nstr(7*b, 60))
    print("  |m-7b| =", nstr(fabs(m1-7*b), 8), f"  ({time.time()-t0:.0f}s)"); print()

    print("[2] m(y^2+(x^2+2x-1)y+x^3) vs 5 b_11  [Brunault, proven]")
    t0 = time.time()
    m2 = mahler_quady(A2, B2, C2)
    print("  m   =", nstr(m2, 60)); print("  5b  =", nstr(5*b, 60))
    print("  |m-5b| =", nstr(fabs(m2-5*b), 8), f"  ({time.time()-t0:.0f}s)"); print()

    print("[3] OPEN: split integral for S_0 = y^2+(x^2+1)y+x^3 vs -b_11")
    t0 = time.time()
    g = lambda th: log(fabs(yminus_S0(exp(mpc(0,1)*th))))
    I1 = mp.quad(g, [0, pi/2]) / pi
    I2 = mp.quad(g, [pi/2, pi]) / pi
    Isplit = I1 - I2
    print("  I1 (0..pi/2)      =", nstr(I1, 60))
    print("  I2 (pi/2..pi)     =", nstr(I2, 60))
    print("  I_split = I1 - I2 =", nstr(Isplit, 60))
    print("  -b_11             =", nstr(-b, 60))
    print("  |I_split+b_11|    =", nstr(fabs(Isplit+b), 8), f"  ({time.time()-t0:.0f}s)"); print()

    print("[4] m(S_0) itself (Boyd: seemingly not r*b_11)")
    t0 = time.time()
    mS0 = mahler_quady(A3, B3, C3)
    print("  m(S_0) =", nstr(mS0, 60))
    print("  m(S_0)/b_11 =", nstr(mS0/b, 30), f"  ({time.time()-t0:.0f}s)")
    r = pslq([mS0, b], tol=mpf(10)**(-70), maxcoeff=10**8)
    print("  PSLQ(m(S0), b_11), maxcoeff=1e8:", r, "-> nil = NOT a rational multiple" if r is None else "-> RELATION!")

if __name__ == "__main__":
    main()
