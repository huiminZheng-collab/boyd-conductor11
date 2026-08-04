"""
Attack 3: push precision to 300 dps; fix the second-model branch tracking;
verify |y_-(e^{iθ})| <= 1 on [0,pi] (makes I1+I2 = -m(S_0) a *theorem*, see report).
"""
from mpmath import mp, mpf, mpc, pi, log, sqrt, exp, fabs, nstr
import sys, time
sys.path.insert(0, ".")
from b11 import b11

DPS = 300
mp.dps = DPS

def yminus_S0(x):
    return -(x**2+1)*(mpf('0.5') - sqrt(mpf('0.25') - x**3/(x**2+1)**2))

def main():
    t0 = time.time()
    b = b11(DPS)
    g = lambda th: log(fabs(yminus_S0(exp(mpc(0,1)*th))))
    I1 = mp.quad(g, [0, pi/2]) / pi
    I2 = mp.quad(g, [pi/2, pi]) / pi
    Is = I1 - I2
    print(f"=== {DPS} dps ===")
    print("I_split =", nstr(Is, 260))
    print("b_11    =", nstr(b, 260))
    print("|I_split - b_11| =", nstr(fabs(Is-b), 6), f"({time.time()-t0:.0f}s)")

    # max |y_-| on [0,pi] (theorem support for I1+I2 = -m(S_0))
    import math
    mx = 0.0; arg = None
    for k in range(20001):
        th = math.pi * k / 20000
        v = abs(complex(yminus_S0(exp(mpc(0,1)*th))))
        if v > mx: mx, arg = v, th
    print(f"max |y_-(e^{{iθ}})| on [0,pi] = {mx:.15f} at θ={arg:.6f}  (<=1 confirms |y_-|<=1 everywhere)")

    # [C] redo: P' = y^2+y+x^3+x^2, track continuous branch y_2 that is >1 near θ=0
    # and ->0 at θ=π. Signed integral (1/pi)∫ log|y_2| =? b_11 ; m(P') =? (1/pi)∫|log|y_2||
    mp.dps = 100
    b100 = b11(100)
    def rootsP(x):
        d = sqrt(1-4*x**3-4*x**2)
        return (-1+d)/2, (-1-d)/2
    N = 4000
    dth = pi / N
    # initial: both roots have |y|=sqrt2 at θ=0; pick either (they are conjugates, same |.|)
    prev = rootsP(exp(mpc(0,1)*dth*mpf('0.5')))[0]
    sg = mpf(0); ab = mpf(0); cross = []
    inbig = True
    for k in range(1, N+1):
        th = dth*(k - mpf('0.5'))
        r = rootsP(exp(mpc(0,1)*th))
        # continuous choice: closest to prev
        y = r[0] if fabs(r[0]-prev) <= fabs(r[1]-prev) else r[1]
        prev = y
        lv = log(fabs(y))
        sg += lv; ab += fabs(lv)
        if (lv > 0) != inbig:
            cross.append(float(th)); inbig = not inbig
    sg *= dth/pi; ab *= dth/pi
    print("\n[C] P' = y^2+y+x^3+x^2, continuous branch y_2:")
    print("  signed (1/pi)∫log|y_2| =", nstr(sg, 60))
    print("  b_11                   =", nstr(b100, 60))
    print("  |signed - b_11| =", nstr(fabs(sg-b100), 6))
    print("  (1/pi)∫|log|y_2||    =", nstr(ab, 60))
    print("  |y_2(pi)| =", nstr(fabs(prev), 6), " crossings of |y_2|=1 at θ ≈", cross)

if __name__ == "__main__":
    main()
