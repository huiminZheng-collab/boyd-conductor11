# MAIN EVALUATION: int_{gamma^-} eta(x,y) via Brunault Thm 1.
#   gamma^- = {0,3/11} - {0,8/11}  (anti-invariant cycle on X_1(11), period wim)
#   Manin-symbol decomposition (continued fractions):
#     {0,3/11} = [[1,0],[3,1]] - [[1,1],[3,4]] + [[3,1],[11,4]]
#     {0,8/11} = [[1,0],[1,1]] - [[1,2],[1,3]] + [[3,2],[4,3]] - [[3,8],[4,11]]
#   x o pi = - G4 G5 / G2^2,  y o pi = G1 G5^3 / (G2^3 G3)   (G_a = prod_b g_{a,b})
#   => row-uniform symbols rho_a = {2:-2,4:1,5:1}, sigma_c = {1:1,2:-3,3:-1,5:3}
#   int_{(alpha)} eta = sum rho_a sigma_c sum_{b,d} [ L((a',d'),(b',-c')) + L((a',-d'),(b',c')) ]
#   with (a',b') = (a,b)alpha, (c',d') = (c,d)alpha,
#   L(P,Q) = Lambda(e_P e_Q, 0),  W_121(e_P e_Q) = -(1/121) f_P f_Q.
from mpmath import mp, mpf, mpc, pi, exp, sqrt, nstr, e1 as E1
from fractions import Fraction as Fr

mp.dps = 50
N = 11; M = 121; NMAX = 200

def e_coeffs(a, b):
    a %= N; b %= N
    c = [0]*(NMAX+1)
    if a == 0 and b == 0: c[0] = Fr(0)
    elif a == 0: c[0] = Fr(1, 2) - Fr(b, N)
    elif b == 0: c[0] = Fr(1, 2) - Fr(a, N)
    else: c[0] = Fr(0)
    for n in range(1, NMAX+1):
        s = 0
        for m in range(1, n+1):
            if n % m == 0:
                nn = n//m
                if m % N == a and nn % N == b: s += 1
                if m % N == (-a) % N and nn % N == (-b) % N: s -= 1
        c[n] = s
    return c

zeta = exp(2*pi*mpc(0, 1)/N)
def f_coeffs(a, b):
    a %= N; b %= N
    z = zeta
    if a == 0 and b == 0: c0 = mpc(0)
    elif a == 0: c0 = (1+z**b)/(2*(1-z**b))
    elif b == 0: c0 = (1+z**a)/(2*(1-z**a))
    else: c0 = ((1+z**a)/(1-z**a) + (1+z**b)/(1-z**b))/2
    c = [c0]
    for n in range(1, NMAX+1):
        s = mpc(0)
        for m in range(1, n+1):
            if n % m == 0:
                nn = n//m
                s += z**(a*m + b*nn) - z**(-(a*m + b*nn))
        c.append(s)
    return c

print("building e/f coefficient tables...", flush=True)
EC = {(a, b): e_coeffs(a, b) for a in range(N) for b in range(N)}
FC = {(a, b): f_coeffs(a, b) for a in range(N) for b in range(N)}

def conv(A, B):
    C = [0]*(NMAX+1)
    for i in range(NMAX+1):
        Ai = A[i]
        if Ai == 0: continue
        for j in range(NMAX+1-i):
            if B[j] != 0:
                C[i+j] += Ai*B[j]
    return C

sq = sqrt(M); cc = 1/sq
def Lambda0(A, B):
    a0 = A[0]; b0 = B[0]
    S1 = mpc(0)
    for n in range(1, NMAX+1):
        if A[n] != 0:
            S1 += A[n]*E1(2*pi*n*cc)
    S2 = mpc(0)
    for n in range(1, NMAX+1):
        if B[n] != 0:
            S2 += B[n]*exp(-2*pi*n*cc)*(sq/(2*pi*n) + M/(2*pi*n)**2)
    return S1 + S2 - b0/2

pairLam = {}
def Lam_pair(P, Q):
    key = (P, Q) if P <= Q else (Q, P)
    v = pairLam.get(key)
    if v is None:
        A = conv(EC[P[0] % N][P[1] % N] if False else EC[(P[0] % N, P[1] % N)], EC[(Q[0] % N, Q[1] % N)])
        Bf = conv(FC[(P[0] % N, P[1] % N)], FC[(Q[0] % N, Q[1] % N)])
        B = [-x/121 for x in Bf]
        v = Lambda0(A, B)
        pairLam[key] = v
    return v

rho = {2: -2, 4: 1, 5: 1}
sig = {1: 1, 2: -3, 3: -1, 5: 3}
symbols = [(1, (1, 0, 3, 1)), (-1, (1, 1, 3, 4)), (1, (3, 1, 11, 4)),
           (-1, (1, 0, 1, 1)), (1, (1, 2, 1, 3)), (-1, (3, 2, 4, 3)),
           (1, (3, 8, 4, 11))]

total = mpc(0)
for sgn, (r, s, t, u) in symbols:
    Iv = mpc(0)
    for a, ra in rho.items():
        for c, sc in sig.items():
            for b in range(11):
                ap = (a*r + b*t) % N
                bpp = (a*s + b*u) % N
                for d in range(11):
                    cp = (c*r + d*t) % N
                    dp = (c*s + d*u) % N
                    Iv += ra*sc*(Lam_pair((ap, dp), (bpp, (-cp) % N))
                                 + Lam_pair((ap, (-dp) % N), (bpp, cp)))
    print("symbol (%+d) [[%d,%d],[%d,%d]]: I = %s   (pairs cached: %d)"
          % (sgn, r, s, t, u, nstr(Iv, 40), len(pairLam)), flush=True)
    total += sgn*Iv

print("="*70)
print("int_{gamma^-} eta(x,y) =", nstr(total, 46))
print("+2 pi b11 = +0.95596868542165847877477217156086908050293759...")
print("-2 pi b11 = -0.95596868542165847877477217156086908050293759...")
print("ratio to 2 pi b11:", nstr(total/0.95596868542165847877477217156086908050293759, 30))
