"""
D_E numerics for the (C3) proof assembly (no PARI needed here).

E = 11.a3, lattice basis w1 (real), tau = 0.5 + 0.229878021224650476... i
=> q = e^{2 pi i tau} is REAL NEGATIVE.
P = (0,0) sits at u/w1 = 3/5, 2P = (1,-1) at 1/5  (z-values on unit circle).

D_E(S) = sum_{n in Z} D(q^n z_S),  D = Bloch-Wigner,
sum over negative n via D(w) = -D(1/w).

Predictions (as of attack 13, see notes/attack13-bertin-diamond.txt):
  (a) exotic relation: D_E(2P) = (3/2) D_E(P)      [Bertin, cited in Brunault (3.151);
      confirmed numerically to 60 digits -- an early draft of this file said "expect 2",
      which was wrong]
  (b) diamond value:   D_E((x)⋄(y)) = (5/2) D_E(P) = pi b_11, while the
      certified regulator integral is 2 pi b_11 (the factor-2 phenomenon of
      paper Remark rem:diamondk0; curve-level, symbol-independent).
      Earlier version "(x).(y) = 5(A)-5(2A) => -5 D_E(P) = 2 pi b_11" mixed
      the raw and Z[E]--reduced expansions; the reduced diamond is
      6(O)-5(A)+5(2A).
"""
from mpmath import mp, mpf, mpc, pi, exp, log, polylog, nstr

mp.dps = 60

tau = mpc(mpf('0.5'), mpf('0.229878021224650476137525642435318986382'))
q = exp(2*pi*mpc(0,1)*tau)
zP  = exp(2*pi*mpc(0,1)*mpf(3)/5)
z2P = exp(2*pi*mpc(0,1)*mpf(1)/5)

def D(w):
    if w == 0:
        return mpf(0)
    return (polylog(2, w).imag + mpc(0,1).phase*0) + 0  # placeholder, replaced below

def bloch_wigner(w):
    """D(w) = Im Li_2(w) + arg(1-w) log|w|,  real-valued."""
    if abs(w) < 1e-50:
        return mpf(0)
    return polylog(2, w).imag + mp.arg(1-w)*log(abs(w))

def D_E(z):
    s = mpf(0)
    # n >= 0: D(q^n z)
    n = 0
    while True:
        w = q**n * z
        t = bloch_wigner(w)
        s += t
        if abs(w) < mpf('1e-65'):
            break
        n += 1
    # n < 0: D(q^n z) = -D(q^{-n} / z)
    n = 1
    while True:
        w = q**n / z
        t = bloch_wigner(w)
        s -= t
        if abs(w) < mpf('1e-65'):
            break
        n += 1
    return s

DEP  = D_E(zP)
DE2P = D_E(z2P)
print("q =", nstr(q, 20))
print("D_E(P)  =", nstr(DEP, 40))
print("D_E(2P) =", nstr(DE2P, 40))
print("(a) D_E(2P)/D_E(P) =", nstr(DE2P/DEP, 25), "  (exotic: 3/2)")
b11 = mpf('0.1521471417259180494862272974786344956281')
print("(b) (5/2) D_E(P) =", nstr(mpf(5)/2*DEP, 40))
print("    pi b_11     =", nstr(pi*b11, 40), "  <- D_E of the reduced diamond")
print("    2 pi b_11   =", nstr(2*pi*b11, 40), "  <- certified regulator integral (factor 2)")
L2 = mpf('0.5460480362150135183341266604334443385907')
print("(c) D_E(P)/L(E,2) =", nstr(DEP/L2, 25), "  vs 11/(10 pi) =", nstr(11/(10*pi), 25))
