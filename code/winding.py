"""
Winding number of the signed cycle gamma~ in H_1(E, Z)^-  (S_0 curve, k=0).

Method: period pairing with the invariant differential w = dx/u on the
quartic model u^2 = g(x) = x^4 - 4x^3 + 2x^2 + 1 (u = 2*y_big + x^2 + 1).

Findings (60 dps):
  * the geometric path |x|=1 carrying y_big is DISCONTINUOUS at the torus
    crossings th = +-pi/2 (jumps between the two crossing points) -- the
    plain loop integral I_loop = -0.47447...i is NOT a period of any
    integral cycle (ratio 0.16262... vs w_anti, irrational);
  * the SIGNED cycle (Samart Lemma 9 weights +1 on (-c,c), -1 outside)
    has period I_signed = -2.917633233876990458...i = w_anti EXACTLY,
    i.e. gamma~ = generator of H_1(E, Z)^-  (winding n = 1).

Cross-check: PARI E.omega for 11.a3 [0,-1,1,0,0] gives
w_anti = 2i*Im(w2) = -2.91763323387699045866177922580735051...i
and w1 = 6.34604652139776710844... = our oval period Om_re (model constant 1).
"""
from mpmath import mp, mpf, mpc, pi, sqrt, exp, nstr, quad

mp.dps = 60

def ybig(th):
    x = exp(mpc(0, 1)*th)
    B = x**2 + 1
    d = sqrt(B**2 - 4*x**3)
    a, b = (-B + d)/2, (-B - d)/2
    return a if abs(a) >= abs(b) else b

def integrand(th):
    x = exp(mpc(0, 1)*th)
    u = 2*ybig(th) + x**2 + 1
    return mpc(0, 1)*x / u

if __name__ == "__main__":
    c = pi/2
    print("jump test at crossings (y_big left/right limits):")
    for side, th in [("c-", c-1e-8), ("c+", c+1e-8), ("-c+", -c+1e-8), ("-c-", -c-1e-8)]:
        print(f"  y_big({side}) =", nstr(ybig(th), 20))
    inner = quad(integrand, [-c, 0]) + quad(integrand, [0, c])
    outer = quad(integrand, [-pi, -c]) + quad(integrand, [c, pi])
    I_loop, I_signed = inner + outer, inner - outer
    w_anti = mpc(0, -2.91763323387699045866177922580735051431848685790533032293924)
    print("\nI_loop   =", nstr(I_loop, 30), " (not a cycle period)")
    print("I_signed =", nstr(I_signed, 30))
    print("w_anti   =", nstr(w_anti, 30), " (PARI 2i*Im(w2), 11.a3)")
    print("I_signed/w_anti =", nstr(I_signed/w_anti, 25), "  =>  winding n = 1")
