"""
CERTIFICATE numerics for the cycle lemma, S_1 (k=1, conductor 17).
Mirrors n1_certify.py (k=0, conductor 11) step by step.

Chain algebra (same as k=0, see notes/proof-k1.md):
  gamma~ = +y_big[-c,c] - y_big(outside), c = 2*pi/3  (fold: arccos(-k/2), k=1)
  Jumps at +-c EXACT: at x = exp(+-2 pi i/3) = omega, omega^2 one has
  B = x^2+x+1 = 0, so S_1 = y^2 + x^3 = y^2 + 1 = 0  =>  y^2 = -1, y = +- i.
  P = (omega, i), P~ = (omega^2, -i) = conjugate, -P = (omega, -i).
  beta0 = alpha1 + alpha2 (small-branch closing arcs, same as k=0).
  C' = gamma~ + beta0: closed, anti-invariant, integral.

KEY DIFFERENCE from k=0: minimal model [1,-1,1,-1,0] has disc = +17 > 0
(two real components), conj(w2) = -w2 exactly, so the primitive
anti-invariant period is w_anti = w2 itself (k=0: 2*w2 - w1 = 2i Im w2).
We print the ratio against both w2 and 2*w2.

D(z) = (z^2+z+1)^2 - 4 z^3 = z^4 - 2 z^3 + 3 z^2 + 2 z + 1;  D(1) = 5 != 0,
so NO 1/sqrt(theta) singularity at theta=0 (unlike k=0).

w2, b_17 from PARI (code/k1_pari.gp, 80 digits).
"""
from mpmath import mp, mpf, mpc, pi, sqrt, exp, nstr, quad, fabs, log

for dps in [60, 80]:
    mp.dps = dps
    def roots(th):
        x = exp(mpc(0, 1)*th)
        B = x**2 + x + 1
        d = sqrt(B**2 - 4*x**3)
        a, b = (-B + d)/2, (-B - d)/2
        return (a, b) if fabs(a) >= fabs(b) else (b, a)
    def per(th, br):
        x = exp(mpc(0, 1)*th)
        y = roots(th)[0 if br == 'b' else 1]
        return mpc(0, 1)*x / (2*y + x**2 + x + 1)
    c = 2*pi/3
    inner_b = quad(lambda t: per(t,'b'), [-c,0]) + quad(lambda t: per(t,'b'), [0,c])
    outer_b = quad(lambda t: per(t,'b'), [-pi,-c]) + quad(lambda t: per(t,'b'), [c,pi])
    I_signed = inner_b - outer_b
    A_in  = quad(lambda t: per(t,'s'), [-c,0]) + quad(lambda t: per(t,'s'), [0,c])
    A_out = quad(lambda t: per(t,'s'), [-pi,-c]) + quad(lambda t: per(t,'s'), [c,pi])
    # PARI minimal model [1,-1,1,-1,0]: w1 real, w2 pure imaginary, disc=+17
    w2 = mpc(0, -2.7457391180897536720341878803801090713441744772379115706927434014592934271565948)
    w_anti = w2                      # primitive anti-invariant period (disc>0)
    P = I_signed + A_out - A_in
    print(f"--- dps={dps} ---")
    print("I_signed          =", nstr(I_signed, 40))
    print("I_signed/w_anti   =", nstr(I_signed/w_anti, 40))
    print("period(C')        =", nstr(P, 40))
    print("period(C')/w_anti =", nstr(P/w_anti, 40))
    print("period(C')/(2 w2) =", nstr(P/(2*w2), 40))
    for n in range(-4, 5):
        if n and fabs(P/w_anti - n) < mpf('0.1'):
            print(f"  candidate integer n = {n}, |P/w_anti - n| =",
                  nstr(fabs(P/w_anti - n), 6))
    # regulator side: int_gamma~ eta = 2(J1 - J2), exact algebra as k=0
    J1 = quad(lambda t: log(fabs(roots(t)[1])), [0, c])
    J2 = quad(lambda t: log(fabs(roots(t)[1])), [c, pi])
    b17 = mpf('0.29935558688291539005379974769003145062973595594493634132226823157334179816452589')
    int_eta = 2*(J1 - J2)
    print("int_gamma~ eta =", nstr(int_eta, 40))
    print("  /(2 pi b_17) =", nstr(int_eta/(2*pi*b17), 40))
    print("  ntilde(1) = -(J1-J2)/pi =", nstr(-(J1-J2)/pi, 40), " vs b_17")
    print("  (structural identity ntilde(1) = -(1/2pi) int_gamma~ eta, k=1 sign)")

# --- exact jump values at the corners (mpmath 60) ---
mp.dps = 60
c = 2*pi/3
omega = exp(mpc(0,1)*c)
def yvals(th):
    x = exp(mpc(0,1)*th)
    B = x**2 + x + 1
    d = sqrt(B**2 - 4*x**3)
    a, b = (-B + d)/2, (-B - d)/2
    return (a, b) if fabs(a) >= fabs(b) else (b, a)
print("\ncorner checks (c = 2 pi/3, omega = exp(2 pi i/3)):")
print("omega^2+omega+1 =", nstr(omega**2+omega+1, 5), " omega^3 =", nstr(omega**3, 5))
for side, th in [("c^-", c-1e-40), ("c^+", c+1e-40)]:
    yb, ys = yvals(mpf(th) if side=="c^-" else mpf(th))
    print(f"y_big({side}) = {nstr(yb,20)}  y_small({side}) = {nstr(ys,20)}")
print("exact candidates: +i, -i; |y_big*y_small| at c:",
      nstr(fabs(yvals(c)[0]*yvals(c)[1]), 10))
