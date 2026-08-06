"""attack16 step 3: numerical computation of I_c = int_0^{c/11} 2 pi i f_11 dz,
to calibrate PARI's mseval normalization.  f_11 = eta(tau)^2 eta(11 tau)^2.
Method: regularized vertical integrals with Fricke acceleration (exact formulas,
q-series with exponentially small tails)."""
from mpmath import mp, mpf, mpc, pi, exp, sqrt
mp.dps = 60

# q-expansion of f_11 via eta products: f = q * prod (1-q^n)^2 (1-q^(11n))^2
NMAX = 3000
import math
def f11_coeffs(nmax):
    # log f / q = 2 sum log(1-q^n) + 2 sum log(1-q^(11n))
    a = [mpf(0)]*(nmax+1)
    a[1] = mpf(1)  # f = q * prod ...; compute via series multiplication
    # prod_{n>=1} (1 - q^n)^2 (1 - q^(11 n))^2 up to q^nmax
    P = [mpf(0)]*(nmax+1); P[0] = mpf(1)
    for n in range(1, nmax+1):
        # multiply by (1-q^n)^2
        for _ in range(2):
            for k in range(nmax, n-1, -1):
                P[k] -= P[k-n]
        if 11*n <= nmax:
            # multiply by (1-q^(11n))^2
            for _ in range(2):
                for k in range(nmax, 11*n-1, -1):
                    P[k] -= P[k-11*n]
    # f = q * P: a[n] = P[n-1]
    return [mpf(0)] + P

a = f11_coeffs(NMAX)
print("a[1..12] =", [int(a[n]) for n in range(1,13)])

def S(x):
    """sum a_n/n x^n, |x|<1"""
    s = mpf(0)
    for n in range(1, NMAX+1):
        s += a[n]*x**n/n
    return s

epi = lambda t: exp(-2*pi*t)   # e^{-2 pi t}
zeta11 = exp(2*pi*mpc(0,1)/11)

# Fricke eigenvalue: determine numerically: f(-1/(11 tau)) = eps * 11 tau^2 f(tau)
tau = mpc(0, 1.3)
def fval(tau):
    q = exp(2*pi*mpc(0,1)*tau)
    return sum(a[n]*q**n for n in range(1, 400))
eps = fval(-1/(11*tau)) / (11*tau**2*fval(tau))
print("Fricke eigenvalue eps =", eps)   # expect +-1
eps = int(round(eps.real))

# reg int_0^{i oo} omega = eps * S(e^{-2pi/11}) - S(e^{-2pi})
J0 = eps*S(epi(mpf(1)/11)) - S(epi(1))
print("reg int_0^{i inf} =", J0)

# L(f,1) via lfun-independent formula: L(f,1) = sum a_n/n e^{-2pi n/sqrt(11)} * (1+eps)?  just print J0
# I_c for c = k/11:
# reg int_{k/11}^{i inf} = reg int_0^{i inf} + S(e^{2 pi i (i/121 - d/11)}) - S(e^{2 pi i k/11} e^{-2 pi})
# where d = k^{-1} mod 11  (gamma_k = [[k, b],[11, d]], k d - 11 b = 1)
def Ic(k):
    d = pow(k, -1, 11)
    t1 = mpc(-d, 0)/11 + mpc(0,1)/121
    x1 = exp(2*pi*mpc(0,1)*t1)          # e^{2 pi i (i/121 - d/11)}
    x2 = zeta11**k * epi(1)             # e^{2 pi i k/11} e^{-2 pi}
    regk = J0 + S(x1) - S(x2)
    return J0 - regk

w1 = mpf("6.3460465213977671084439730837727365260974612091530433398356846061520416930638886")
w2 = mpc("3.1730232606988835542219865418863682630487306045765216699178423030760208465319443",
         "-1.4588166169384952293308896129036752571592434289526651614696187624505378966090287")
wim = w1 - 2*w2
phi = {1:(4,25,0,1), 2:(3,50,1,2), 3:(-7,50,1,2), 4:(-6,25,0,1), 5:(-6,25,0,1),
       6:(-7,50,-1,2), 7:(3,50,-1,2), 8:(4,25,0,1), 9:(-1,25,0,1), 10:(0,1,0,1)}
for k in range(1,11):
    v = Ic(k)
    p, q, r, s = phi[k]
    cand = (mpf(p)/q)*w1 + (mpf(r)/s)*wim
    print(f"k={k:2d}  I_c = {v}")
    print(f"      cand w1*phi+ + wim*phi- = {cand}   diff = {v-cand}")
