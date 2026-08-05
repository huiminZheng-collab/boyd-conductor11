"""
Attack 13 (review strengthening): (C3) from 149 digits to ~300 digits.

I_split = (1/pi) ∫_0^{pi/2} log|y_-(e^{iθ})| dθ - (1/pi) ∫_{pi/2}^pi log|y_-(e^{iθ})| dθ
compared against b_11 = L'(E11,0) from the eta-product Mellin series (b11.py),
both at 800 dps. Empirical quadrature error scales ~ 10^{-0.48*dps},
so 800 dps should give ~ 350 certified digits of agreement.
Cross-check: PARI/GP lfun(E11, 0 derivative) at 320 digits (see notes/attack13-c3-300.txt).
"""
from mpmath import mp, mpf, mpc, pi, log, sqrt, exp, fabs, nstr
import sys, time
sys.path.insert(0, "code")
from b11 import b11

DPS = 800
mp.dps = DPS

def yminus_S0(x):
    return -(x**2+1)*(mpf('0.5') - sqrt(mpf('0.25') - x**3/(x**2+1)**2))

def main():
    t0 = time.time()
    b = b11(DPS)
    print(f"b_11 ({DPS} dps series, nmax auto):")
    print(nstr(b, 310))
    g = lambda th: log(fabs(yminus_S0(exp(mpc(0,1)*th))))
    I1 = mp.quad(g, [0, pi/2]) / pi
    print(f"I1 done ({time.time()-t0:.0f}s)")
    I2 = mp.quad(g, [pi/2, pi]) / pi
    Is = I1 - I2
    print("I_split:")
    print(nstr(Is, 310))
    err = fabs(Is - b)
    print(f"|I_split - b_11| = {nstr(err, 6)}   ({time.time()-t0:.0f}s total)")
    print(f"agreement in digits: {nstr(-mp.log10(err), 8)}")

if __name__ == "__main__":
    main()
