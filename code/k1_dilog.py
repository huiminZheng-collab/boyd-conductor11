"""
D_E numerics for the k=1 (conductor 17) regulator assembly.
Mirrors dilog.py (k=0).

Minimal model E: y^2 + x y + y = x^3 - x^2 - x  ([1,-1,1,-1,0], disc = +17).
disc > 0: w1 real, w2 pure imaginary; use tau = -w2/w1 = 0.8873941733... i
(Im > 0), q = e^{2 pi i tau} = 0.003788966335379...  REAL POSITIVE.

Torsion Z/4, generator P = (0,0):  2P = (1,-1) [2-torsion], 3P = (0,-1).
PARI ellpointtoz (code/k1_zvals.gp):
  z(P)/w1  = 1/4 + tau/2   =>  z_q = e^{2 pi i (1/4 + tau/2)} = i * q^{1/2}
  z(2P)/w1 = 1/2           =>  z_q = -1

Diamond product (exact algebra on the S_1 cubic model, Abel-verified,
notes/proof-k1.md):   (x) <> (y) = 6(O) + 4(A) - 6(2A)   in Z[E]^-
=>  D_E(<>) = 4 D_E(A) - 6 D_E(2A).

Predictions to check:
  (a) D_E(2A) = 0            (2-torsion, q real => all terms vanish)
  (b) Bloch synthesis (factor-1 normalization of Lalin--Ramamonjisoa Thm. 6):
      int_{gamma^-} eta = +-D_E(<>) = +-4 D_E(A)
      vs measured int_{gamma^-} eta = -2 pi b_17  (k1_certify.py)
      => expect D_E(A) = +- pi b_17 / 2 = +-0.47022...
  (c) Brunault-style constant: D_E(A)/L(E,2) =? rational/pi
      with L(E,2) = 0.69518146300948461345173162509571731230...
"""
from mpmath import mp, mpf, mpc, pi, exp, log, polylog, nstr, sqrt

mp.dps = 60

tau = mpc(0, mpf('0.88739417337317839326832472783180088710153577029494590500547192561034443944405467'))
q = exp(2*pi*mpc(0,1)*tau)
zA  = exp(2*pi*mpc(0,1)*(mpf(1)/4 + tau/2))     # i * sqrt(q)
z2A = mpf(-1)                                     # 2-torsion

def bloch_wigner(w):
    """D(w) = Im Li_2(w) + arg(1-w) log|w|,  real-valued."""
    if abs(w) < 1e-50:
        return mpf(0)
    return polylog(2, w).imag + mp.arg(1-w)*log(abs(w))

def D_E(z):
    s = mpf(0)
    n = 0
    while True:
        w = q**n * z
        s += bloch_wigner(w)
        if abs(w) < mpf('1e-65'):
            break
        n += 1
    n = 1
    while True:
        w = q**n / z
        s -= bloch_wigner(w)
        if abs(w) < mpf('1e-65'):
            break
        n += 1
    return s

DEA  = D_E(zA)
DE2A = D_E(z2A)
b17 = mpf('0.29935558688291539005379974769003145062973595594493634132226823157334179816452589')
L2  = mpf('0.69518146300948461345173162509571731229515793682980952137951318196803875189356525')

print("q =", nstr(q, 25))
print("(a) D_E(2A) =", nstr(DE2A, 20), "   (expect 0)")
print("    D_E(A)  =", nstr(DEA, 40))
print("(b) 4 D_E(A) =", nstr(4*DEA, 40))
print("    2 pi b_17 =", nstr(2*pi*b17, 40))
print("    4 D_E(A) / (2 pi b_17) =", nstr(4*DEA/(2*pi*b17), 25))
print("    D_E(A) / (pi b_17) =", nstr(DEA/(pi*b17), 25))
print("(c) D_E(A)/L(E,2) =", nstr(DEA/L2, 25))
print("    vs 17/(8 pi) =", nstr(17/(8*pi), 25), "  17/(16 pi) =", nstr(17/(16*pi), 25))
print("    pi D_E(A)/L(E,2) =", nstr(pi*DEA/L2, 25))
# also: the full diamond combination and Bloch form
Dd = 4*DEA - 6*DE2A
print("    D_E(<>) = 4 D_E(A) - 6 D_E(2A) =", nstr(Dd, 40))
print("    D_E(<>) / (2 pi b_17) =", nstr(Dd/(2*pi*b17), 25), "  (factor-1 Bloch: expect +-1)")
