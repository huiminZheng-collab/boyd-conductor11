# Abel values D_k = int_0^{k/11} omega on X_1(11) (omega = 2 pi i f11(tau) dtau)
# Piecewise convergent path, no regularization ambiguity:
#   cusp 0   -> i*H0        : Fricke conjugate  tau' = -1/(11 tau)   (exact series)
#   i*H0     -> i*H         : direct q-series quadrature
#   i*H      -> k/11 + i*H  : direct q-series quadrature
#   k/11+i*H -> k/11+i*h1   : direct q-series quadrature
#   k/11+i*h1 -> k/11       : Gamma_0(11) conjugate gamma_k=[[k,b],[11,d]] (exact series)
# Cross-check via S(k)-S(0) (vertical paths through i*inf), and L(f,1).
from mpmath import mp, mpf, mpc, pi, exp, nstr, quad

mp.dps = 60
NMAX = 2500

P = [0]*(NMAX+1); P[0] = 1
for n in range(1, NMAX+1):
    for _ in range(2):
        for k in range(NMAX, n-1, -1):
            P[k] -= P[k-n]
    if 11*n <= NMAX:
        for _ in range(2):
            for k in range(NMAX, 11*n-1, -1):
                P[k] -= P[k-11*n]
a = [0] + P   # a[n] = coefficient of q^n in f11
assert a[1:13] == [1,-2,-1,2,1,2,-2,0,-2,-2,1,-2], a[1:13]
print("coefficients ok")

TWO_PI = 2*pi
TOL = mpf(10)**(-65)

def fval(re_, y):
    q = exp(TWO_PI*mpc(0,1)*mpc(re_, y))
    s = mpc(0); qn = q
    for n in range(1, NMAX+1):
        s += a[n]*qn
        qn *= q
        if n > 80 and abs(qn) < TOL:
            break
    return s

def S_upper(c, y):
    # integral_{c+iy}^{c+i*inf} omega = -sum a_n e^{2 pi i n c} e^{-2 pi n y}/n
    zz = exp(TWO_PI*mpc(0,1)*c)
    s = mpc(0)
    zn = zz
    for n in range(1, NMAX+1):
        t = a[n]*zn*exp(-TWO_PI*n*y)/n
        s += t
        zn *= zz
        if n > 80 and abs(t) < TOL:
            break
    return -s

H  = mpf('0.75')
H0 = mpf('0.35')   # split for cusp-0 vertical
h1 = mpf('0.35')   # split for cusp-k/11 vertical

def inv_mod(k): return pow(k, -1, 11)

seg1 = S_upper(0, 1/(11*H0))          # int_0^{i H0}
def seg2():                            # int_{i H0}^{i H}
    return quad(lambda y: -TWO_PI*fval(0, y), [H0, H])
def seg3(k):                           # int_{i H}^{k/11 + i H}
    return quad(lambda t: TWO_PI*mpc(0,1)*fval(t, H), [0, mpf(k)/11])
def seg4(k):                           # int_{k/11+i H}^{k/11+i h1}
    return quad(lambda y: TWO_PI*fval(mpf(k)/11, y), [h1, H])
def seg5(k):                           # int_{k/11+i h1}^{k/11}
    d = inv_mod(k)
    return S_upper(mpf(-d)/11, 1/(121*h1))

s2 = seg2()
L_check = seg1 + s2 + S_upper(0, H)    # = int_0^{i inf} = -L(f,1)
print("int_0^iinf omega =", nstr(L_check, 50))
print("  -L(f,1) (PARI) = -0.2538418608559106843377589233509094610439...")

w1  = mpf("6.3460465213977671084439730837727365260974612091530")
wim = mpf("2.9176332338769904586617792258073505143184868579053")  # w1-2w2 = i*wim

D = {}
for k in range(1, 6):
    Dk = seg1 + s2 + seg3(k) + seg4(k) + seg5(k)
    # cross-check: D_k = int_0^iinf - int_{k/11}^iinf
    d = inv_mod(k)
    int_k_inf = -S_upper(mpf(-d)/11, 1/(121*h1)) - seg4(k) + S_upper(mpf(k)/11, H)
    Dk_b = L_check - int_k_inf
    D[k] = Dk
    print("k=%d  D_k = %s" % (k, nstr(Dk, 45)))
    print("      cross = %s   diff = %s" % (nstr(Dk_b, 45), nstr(abs(Dk-Dk_b), 8)))

print()
mmap = {0:0, 1:2, 2:4, 3:1, 4:3}   # (3m/5 mod 1 = r/5) -> m
for k in range(2, 6):
    V = D[k] - D[1]                # = int_{1/11}^{k/11} = z(pi(k/11)) since pi(1/11)=O
    beta = V.imag / wim
    alpha = V.real / w1
    b2 = round(beta*2)             # beta should be in (1/2) Z
    rb = mpf(b2)/2
    resid = (alpha + rb) % 1       # alpha + beta == 3m/5 (mod 1)
    r5 = round(resid*5) % 5
    m = mmap[r5]
    print("k=%d  alpha=%s  beta=%s  (beta*2=%d)  resid=%s  -> m_k = %d"
          % (k, nstr(alpha, 20), nstr(beta, 20), b2, nstr(resid, 20), m))
