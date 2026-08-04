"""
b-values for the S_k family by POINT COUNTING (no LMFDB needed).

E_k: Y^2 Z + (X^2 + k X Z + Z^2) Y + X^3 = 0.
a_p = p + 1 - #E(F_p); two rational points at infinity ([0:1:0], [1:-1:0]).
Bad primes (p | N, multiplicative): a_p = p - #smooth(F_p).
a_n multiplicative; then for w = +1:
  b(N) = L'(E,0) = Lambda(f,2) = sum a_n [ e^{-t}(1/t + 1/t^2) + E_1(t) ], t = 2 pi n / sqrt(N).

Pipeline validated on k=0 against b11.py (eta-product coefficients).
"""
from mpmath import mp, mpf, exp, sqrt, pi, nstr
import sys
sys.path.insert(0, ".")
from b11 import b11

mp.dps = 50

def primes_upto(n):
    sieve = [True]*(n+1)
    ps = []
    for i in range(2, n+1):
        if sieve[i]:
            ps.append(i)
            for j in range(i*i, n+1, i):
                sieve[j] = False
    return ps

def ap_good(k, p):
    """a_p at good prime via affine count + 2 points at infinity."""
    s = 0
    for x in range(p):
        B = (x*x + k*x + 1) % p
        c = x*x*x % p
        if p == 2:                # quadratic formula invalid in char 2: count directly
            s += sum(1 for y in range(2) if (y*y + B*y + c) % 2 == 0)
            continue
        d = (B*B - 4*c) % p
        if d == 0:
            s += 1
        elif pow(d, (p-1)//2, p) == 1:
            s += 2
    return p + 1 - 2 - s          # #E = 2 (infinity) + s

def ap_bad(k, p):
    """a_p at bad prime: #E_ns(F_p) = p - a_p (multiplicative reduction)."""
    smooth = 0
    for x in range(p):
        for y in range(p):
            f  = (y*y + (x*x+k*x+1)*y + x*x*x) % p
            if f == 0:
                fx = (2*x*y + k*y + 3*x*x) % p
                fy = (2*y + x*x + k*x + 1) % p
                if fx or fy:
                    smooth += 1
    # infinity: [0:1:0] and [1:-1:0] are both always smooth
    # (partial wrt Z resp. Y is nonzero there)
    return p - (smooth + 2)

def ap(k, p, N):
    return ap_good(k, p) if N % p else ap_bad(k, p)

def an_list(k, N, nmax):
    ps = primes_upto(nmax)
    a = [0]*(nmax+1)
    a[1] = 1
    for p in ps:
        v = ap(k, p, N)
        if N % p == 0:                    # bad prime: a_{p^e} = v^e
            pe, e = p, 1
            while pe <= nmax:
                a[pe] = v**e
                e += 1; pe *= p
        else:                             # good: a_{p^e} = v a_{p^e-1} - p a_{p^e-2}
            a[p] = v
            f0, f1, pe = 1, v, p*p
            while pe <= nmax:
                f2 = v*f1 - p*f0
                a[pe] = f2
                f0, f1, pe = f1, f2, pe*p
    for n in range(2, nmax+1):            # multiplicative fill, ascending
        if a[n] == 0:
            m, pe = n, 1
            for p in ps:
                if m % p == 0:
                    while m % p == 0:
                        m //= p; pe *= p
                    break
            if pe > 1 and a[m]:
                a[n] = a[pe]*a[m]
    return a

def bvalue(a, N, nmax):
    t = 2*pi/sqrt(N)
    s = mpf(0)
    for n in range(1, nmax+1):
        if a[n]:
            tn = t*n
            s += a[n]*(exp(-tn)*(1/tn + 1/tn**2) + mp.expint(1, tn))
    return s

if __name__ == "__main__":
    nmax = 700
    # ---- validate pipeline on k=0, N=11 (against eta-product b_11)
    a = an_list(0, 11, nmax)
    b0 = bvalue(a, 11, nmax)
    ref = b11(50)
    print("pipeline k=0: b =", nstr(b0, 30))
    print("eta-product : b =", nstr(ref, 30))
    print("  |diff| =", nstr(abs(b0-ref), 6))
    print("  a_p check p=2,3,5,7,11:", [ap(0, p, 11) for p in [2,3,5,7,11]],
          " (expect -2,-1,1,0,-1 for 11.a3)")
    # ---- family targets
    for k, N in [(1, 17), (-1, 53)]:
        a = an_list(k, N, nmax)
        b = bvalue(a, N, nmax)
        print(f"\nk={k}, N={N}: a_p (p=2,3,5,7,11,13,{N}) =",
              [ap(k, p, N) for p in [2,3,5,7,11,13,N]])
        print(f"  b_{N} =", nstr(b, 30))
