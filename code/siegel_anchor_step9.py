# EXACT structure of F_total (rational q-coefficients) behind int_{gamma^-} eta(x,y).
# Same Manin-symbol sum as step8, but with exact Fraction arithmetic.
# Then: analyze the cusp-form part (f11 coefficient) and Eisenstein decomposition.
from fractions import Fraction as Fr

N = 11
NMAX = 250

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

print("e tables...", flush=True)
EC = {(a, b): e_coeffs(a, b) for a in range(N) for b in range(N)}

def conv(P, Q):
    # returns Fraction array; series parts are integral, constants rational
    A = EC[P]; B = EC[Q]
    C = [Fr(0)]*(NMAX+1)
    C[0] = A[0]*B[0]
    for n in range(1, NMAX+1):
        t = 0
        for i in range(1, n):
            t += A[i]*B[n-i]
        C[n] = A[0]*B[n] + A[n]*B[0] + Fr(t)
    return C

pairCache = {}
def pair(P, Q):
    key = (P, Q) if P <= Q else (Q, P)
    v = pairCache.get(key)
    if v is None:
        v = conv(*key)
        pairCache[key] = v
    return v

rho = {2: -2, 4: 1, 5: 1}
sig = {1: 1, 2: -3, 3: -1, 5: 3}
symbols = [(1, (1, 0, 3, 1)), (-1, (1, 1, 3, 4)), (1, (3, 1, 11, 4)),
           (-1, (1, 0, 1, 1)), (1, (1, 2, 1, 3)), (-1, (3, 2, 4, 3)),
           (1, (3, 8, 4, 11))]

F = [Fr(0)]*(NMAX+1)
per_sym = []
for sgn, (r, s, t, u) in symbols:
    Fs = [Fr(0)]*(NMAX+1)
    for a, ra in rho.items():
        for c, sc in sig.items():
            for b in range(N):
                ap = (a*r + b*t) % N
                bpp = (a*s + b*u) % N
                for d in range(N):
                    cp = (c*r + d*t) % N
                    dp = (c*s + d*u) % N
                    w = ra*sc
                    P1 = pair((ap, dp), (bpp, (-cp) % N))
                    P2 = pair((ap, (-dp) % N), (bpp, cp))
                    for n in range(NMAX+1):
                        Fs[n] += w*(P1[n] + P2[n])
    per_sym.append((sgn, (r, s, t, u), Fs))
    for n in range(NMAX+1):
        F[n] += sgn*Fs[n]
    print("symbol (%+d) [[%d,%d],[%d,%d]] done, F_sym[0..6] = %s"
          % (sgn, r, s, t, u, [str(x) for x in Fs[:7]]), flush=True)

print("="*70)
print("F_total coefficients [0..20]:")
for n in range(21):
    print("  q^%d: %s" % (n, F[n]))

# f11 coefficients
Pf = [0]*(NMAX+1); Pf[0] = 1
for n in range(1, NMAX+1):
    for _ in range(2):
        for k in range(NMAX, n-1, -1):
            Pf[k] -= Pf[k-n]
    if 11*n <= NMAX:
        for _ in range(2):
            for k in range(NMAX, 11*n-1, -1):
                Pf[k] -= Pf[k-11*n]
f11 = [0] + Pf   # f11[n]

# E2,11(tau) = E2(tau) - 11 E2(11 tau), E2 = 1 - 24 sum sigma1 q^n
def sigma1(n):
    return sum(m for m in range(1, n+1) if n % m == 0)
E211 = [Fr(0)]*(NMAX+1)
E211[0] = Fr(-10)
for n in range(1, NMAX+1):
    v = -24*sigma1(n)
    if n % 11 == 0:
        v += 264*sigma1(n//11)
    E211[n] = Fr(v)

print("="*70)
# try F = lam f11 + mu E211  (M_2(Gamma_0(11)) is 2-dimensional)
lam = Fr(F[1], f11[1])
mu = Fr(F[0] - lam*f11[0] if False else F[0], E211[0])  # f11[0]=0
ok = all(F[n] == lam*f11[n] + mu*E211[n] for n in range(NMAX+1))
print("M_2(Gamma_0(11)) test: lam =", lam, " mu =", mu, " exact match:", ok)
if not ok:
    bad = [n for n in range(NMAX+1) if F[n] != lam*f11[n] + mu*E211[n]]
    print("  first mismatches:", bad[:8])
    for n in bad[:4]:
        print("   n=%d F=%s f11=%s E211=%s" % (n, F[n], f11[n], E211[n]))

import json
with open("code/siegel_anchor_step9_out.json", "w") as f:
    json.dump({"F": [str(x) for x in F],
               "per_sym": [[sgn, list(al), [str(x) for x in Fs]] for sgn, al, Fs in per_sym]}, f)
print("written code/siegel_anchor_step9_out.json")
