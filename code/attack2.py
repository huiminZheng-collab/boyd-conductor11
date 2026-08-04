"""
Boyd conductor-11 attack — Step 2: sharpen and dig deeper.

Findings from attack1 (80 dps):
  * m((1+x)(1+y)(1+x+y)+xy) = 7 b_11   to 52 digits   [Brunault, proven - reproduced]
  * m(y^2+(x^2+2x-1)y+x^3)   = 5 b_11  to 52 digits   [Brunault, proven - reproduced]
  * OPEN CONJECTURE: I_split := I1 - I2 = +b_11 to 52 digits (!!)
    (Samart 2023 eq. (4.1) writes -L'(E,0); sign depends on which root is called y_-)
  * m(S_0) = 0.4056029559150104... (Boyd's 0.4056029), NOT a rational multiple of b_11
  * STRUCTURE: I1 + I2 = -m(S_0) exactly (all 60 digits) =>
        I1 = (b_11 - m(S_0))/2,   I2 = -(b_11 + m(S_0))/2.

This script:
  A. re-verifies I_split = b_11 at 150 dps (~110 digit agreement);
  B. verifies the structural identity I1+I2 = -m(S_0) at 150 dps;
  C. second model y^2+y+x^3+x^2 (Boyd slides): signed integral (1/pi)∫_0^pi log|y_2| =? b_11,
     and its Mahler measure vs m(S_0);
  D. PSLQ hunts for m(S_0) among standard constants (expected negative -> documents that
     m(S_0) needs elliptic dilogarithms, not elementary constants).
"""
from mpmath import (mp, mpf, mpc, pi, log, sqrt, exp, fabs, nstr, pslq,
                    zeta, euler, catalan, phi)
import sys, time

sys.path.insert(0, ".")
from b11 import b11

DPS = 150
mp.dps = DPS

def logp(v):
    lv = fabs(v)
    return log(lv) if lv > 1 else mpf(0)

def mahler_quady(A, B, C):
    def f(th):
        x = exp(2*mpc(0,1)*pi*th)
        A_, B_, C_ = A(x), B(x), C(x)
        disc = sqrt(B_**2 - 4*A_*C_)
        y1 = (-B_ + disc)/(2*A_); y2 = (-B_ - disc)/(2*A_)
        return log(fabs(A_)) + logp(y1) + logp(y2)
    return mp.quad(f, [0, mpf('0.25'), mpf('0.5'), mpf('0.75'), 1])

def yminus_S0(x):
    return -(x**2+1)*(mpf('0.5') - sqrt(mpf('0.25') - x**3/(x**2+1)**2))

def main():
    t0 = time.time()
    b = b11(DPS)
    print(f"=== precision: {DPS} dps ===")
    print("b_11 =", nstr(b, 120)); print()

    print("[A] OPEN conjecture I_split = b_11 at 150 dps")
    g = lambda th: log(fabs(yminus_S0(exp(mpc(0,1)*th))))
    I1 = mp.quad(g, [0, pi/2]) / pi
    I2 = mp.quad(g, [pi/2, pi]) / pi
    Isplit = I1 - I2
    print("  I_split =", nstr(Isplit, 120))
    print("  b_11    =", nstr(b, 120))
    print("  |I_split - b_11| =", nstr(fabs(Isplit-b), 6), f"  ({time.time()-t0:.0f}s)"); print()

    print("[B] structure: I1+I2 =? -m(S_0)")
    mS0 = mahler_quady(lambda x: mpc(1), lambda x: x**2+1, lambda x: x**3)
    print("  I1+I2   =", nstr(I1+I2, 120))
    print("  -m(S_0) =", nstr(-mS0, 120))
    print("  |diff|  =", nstr(fabs(I1+I2+mS0), 6))
    print("  m(S_0)  =", nstr(mS0, 120), f"  ({time.time()-t0:.0f}s)"); print()

    print("[C] second model P' = y^2+y+x^3+x^2 (Boyd slides p.28)")
    # roots y = (-1 +- sqrt(1-4x^3-4x^2))/2 ; y_2 := branch with |y_2|>=1 on the arc,
    # slides: (1/pi)∫|log|y2|| = 0.4056029 (NOT b_11) but (1/pi)∫ log|y2| = b_11 (50 dp)
    mP = mahler_quady(lambda x: mpc(1), lambda x: mpc(1), lambda x: x**3+x**2)
    print("  m(P')   =", nstr(mP, 120))
    print("  |m(P') - m(S_0)| =", nstr(fabs(mP-mS0), 6))
    def ysmall(th):  # branch with |y|<=1 (principal sqrt gives which? test both)
        x = exp(mpc(0,1)*th)
        d = sqrt(1-4*x**3-4*x**2)
        ya = (-1+d)/2; yb = (-1-d)/2
        return ya if fabs(ya) <= fabs(yb) else yb
    sg = mp.quad(lambda th: log(fabs(ysmall(th))), [0, pi]) / pi
    print("  (1/pi)∫_0^pi log|y_small| =", nstr(sg, 120))
    print("  vs -b_11, |diff| =", nstr(fabs(sg+b), 6), f"  ({time.time()-t0:.0f}s)"); print()

    print("[D] PSLQ: m(S_0) among elementary constants (maxcoeff=1e10)")
    # Smyth's m(1+x+y) = (3√3/4π) L(χ_{-3},2)
    def chi3(n):
        r = n % 3
        return 1 if r == 1 else (-1 if r == 2 else 0)
    Lchi3_2 = mp.nsum(lambda n: chi3(n)/n**2, [1, mp.inf])
    smyth = 3*sqrt(3)/(4*pi) * Lchi3_2
    vec = [mS0, b, log(2), log(3), catalan, smyth]
    names = ["m(S_0)", "b_11", "log2", "log3", "Catalan", "m(1+x+y)"]
    tol = mpf(10)**(-(DPS-40))
    r = pslq(vec, tol=tol, maxcoeff=10**10, maxsteps=10000)
    print("  basis:", names)
    print("  PSLQ result:", r)
    if r is None:
        print("  -> NO relation with coefficients <= 1e10 among these constants")
    print(f"  ({time.time()-t0:.0f}s)")

if __name__ == "__main__":
    main()
