# Sanity checks for Brunault Thm 1 normalization BEFORE the conductor-11 main run.
#   Part A: pin alpha0(e_{a,b}) numerically via Lemma 13:  e_{a,b}(-1/(N^2 tau)) = -tau f_{a,b}(tau)
#   Part B: W-split Lambda(F,0); check Lambda(f11,0) = L'(11a1,0) = b11
#   Part C: Lambda(e_{0,d}e_{b,0} + e_{0,-d}e_{b,0}, 0) = 0 (eta=0 on imaginary axis when a=0)
#   Part D: conductor-14 example of Brunault sec 5.1: total = 4 L'(14a4,0)
#   Part E: direct imaginary-axis integration of int_0^iinf eta(g_u,g_v) vs Thm 1 (mixed u,v)
from mpmath import mp, mpf, mpc, pi, exp, log, sqrt, nstr, quad, e1 as E1

mp.dps = 50
I = mpc(0, 1)

def zN(N): return exp(2*pi*I/N)

# ---------- f_{a,b} (Def 12) ----------
def al0_f(a, b, N):
    a %= N; b %= N
    z = zN(N)
    if a == 0 and b == 0: return mpc(0)
    if a == 0: return (1+z**b)/(2*(1-z**b))
    if b == 0: return (1+z**a)/(2*(1-z**a))
    return ((1+z**a)/(1-z**a) + (1+z**b)/(1-z**b))/2

def f_ab(a, b, N, tau):
    # al0 + sum_{m,n>=1} (z^{am+bn} - z^{-(am+bn)}) q^{mn}
    a %= N; b %= N
    z = zN(N)
    q = exp(2*pi*I*tau)
    s = al0_f(a, b, N) + 0
    qm = q
    m = 1
    while abs(qm) > mpf(10)**(-52):
        x1 = z**b * qm; x2 = z**(-b) * qm
        s += z**(a*m) * x1/(1-x1) - z**(-a*m) * x2/(1-x2)
        m += 1; qm *= q
    return s

# ---------- e_{a,b} (Def 10); alpha0 pinned in Part A ----------
AL0E = {}   # (N, a%N, b%N) -> value, filled by Part A / formula

def al0_e(a, b, N):
    a %= N; b %= N
    if a == 0 and b == 0: return mpf(0)
    if a == 0: return mpf(1)/2 - mpf(b)/N
    if b == 0: return mpf(1)/2 - mpf(a)/N
    return mpf(0)

def e_sums(a, b, N, tau):
    # the two double sums (without constant term)
    a %= N; b %= N
    q = exp(2*pi*I*tau)
    s = mpc(0)
    for sgn, aa, bb in ((1, a, b), (-1, (-a) % N, (-b) % N)):
        m0 = aa if aa else N
        n0 = bb if bb else N
        m = m0
        while True:
            qm = q**m
            if abs(qm) < mpf(10)**(-52): break
            # sum_{n equiv bb} q^{m n} = q^{m n0}/(1 - q^{m N})
            qmN = q**(m*N)
            s += sgn * q**(m*n0) / (1 - qmN)
            m += N
    return s

def e_ab(a, b, N, tau):
    return al0_e(a, b, N) + e_sums(a, b, N, tau)

print("="*70)
print("Part A: pin alpha0(e) via Lemma 13: e(-1/(N^2 tau)) = -tau f(tau)")
for N in (3, 4, 5):
    for (a, b) in ((0, 1), (0, 2), (1, 0), (2, 0), (1, 1), (1, 2), (2, 1)):
        if (a % N == 0) and (b % N == 0): continue
        vals = []
        for t in (mpf('0.6'), mpf('0.83')):
            tau = I*t
            taup = -1/(N*N*tau)
            rhs = -tau * f_ab(a, b, N, tau)
            vals.append(rhs - e_sums(a, b, N, taup))
        v = vals[0]
        print("N=%d (a,b)=(%d,%d): alpha0 = %s  (stability %s)" %
              (N, a, b, nstr(v, 30), nstr(abs(vals[0]-vals[1]), 5)))
print("formula check: al0_e = 1/2 - {b/N} (a=0), = 1/2 - {a/N} (b=0), = 0 (both nonzero)")
for N in (3, 4, 5):
    for (a, b) in ((0, 1), (0, 2), (1, 0), (2, 0), (1, 1), (1, 2), (2, 1)):
        if (a % N == 0) and (b % N == 0): continue
        t = mpf('0.71'); tau = I*t
        v = -tau*f_ab(a, b, N, tau) - e_sums(a, b, N, -1/(N*N*tau))
        err = abs(v - al0_e(a, b, N))
        assert err < mpf(10)**(-40), (N, a, b, err)
print("alpha0 formula CONFIRMED to 40 digits")

# ---------- coefficient arrays ----------
def e_coeffs(a, b, N, nmax):
    c = [mpf(0)]*(nmax+1)
    c[0] = al0_e(a, b, N)
    a %= N; b %= N
    for n in range(1, nmax+1):
        s = 0
        for m in range(1, n+1):
            if n % m: continue
            nn = n//m
            if m % N == a and nn % N == b: s += 1
            if m % N == (-a) % N and nn % N == (-b) % N: s -= 1
        c[n] = mpf(s)
    return c

def f_coeffs(a, b, N, nmax):
    z = zN(N)
    c = [mpc(0)]*(nmax+1)
    c[0] = al0_f(a, b, N)
    a %= N; b %= N
    for n in range(1, nmax+1):
        s = mpc(0)
        for m in range(1, n+1):
            if n % m: continue
            nn = n//m
            s += z**(a*m + b*nn) - z**(-(a*m + b*nn))
        c[n] = s
    return c

def conv(c1, c2, nmax):
    c = [0]*(nmax+1)
    for i in range(nmax+1):
        if c1[i] == 0: continue
        for j in range(nmax+1-i):
            if c2[j] != 0:
                c[i+j] += c1[i]*c2[j]
    return c

# ---------- Lambda(F, 0) via W-split (M = level, W_M F = B series) ----------
def Lambda0(A, B, M, nmax=None):
    # A, B: coefficient arrays (A = F, B = W_M F); weight 2
    nmax = nmax or (len(A)-1)
    a0 = A[0]; b0 = B[0]
    sq = sqrt(M)
    cc = 1/sq
    S1 = mpc(0)
    for n in range(1, nmax+1):
        if A[n] == 0: continue
        S1 += A[n] * E1(2*pi*n*cc)
        if n > 20 and abs(A[n]*E1(2*pi*n*cc)) < mpf(10)**(-55): break
    S2 = mpc(0)
    for n in range(1, nmax+1):
        if B[n] == 0: continue
        t = B[n] * exp(-2*pi*n*cc) * (sq/(2*pi*n) + M/(2*pi*n)**2)
        S2 += t
        if n > 20 and abs(t) < mpf(10)**(-55): break
    return S1 + S2 - b0/2

print("="*70)
print("Part B: Lambda(f11, 0) = L'(11a1,0) = b11 = 0.15214714172591804948622729747863449563")
nmax = 400
# f11 coefficients (level 11 newform), W_11 f11 = -f11
P = [0]*(nmax+1); P[0] = 1
for n in range(1, nmax+1):
    for _ in range(2):
        for k in range(nmax, n-1, -1):
            P[k] -= P[k-n]
    if 11*n <= nmax:
        for _ in range(2):
            for k in range(nmax, 11*n-1, -1):
                P[k] -= P[k-11*n]
A11 = [mpf(0)] + [mpf(x) for x in P]
B11 = [-x for x in A11]   # W_11 f11 = -f11 (Fricke eigenvalue -1)
val = Lambda0(A11, B11, 11)
print("Lambda(f11,0) =", nstr(val, 45))
print("reference b11 = 0.15214714172591804948622729747863449562814...")

print("="*70)
print("Part C: Lambda(e_{0,d} e_{b,0} + e_{0,-d} e_{b,0}, 0) should be 0 (eta=0 on iR+)")
for N, b, d in ((7, 2, 3), (11, 3, 5)):
    M = N*N
    c1 = e_coeffs(0, d, N, nmax); c2 = e_coeffs(b, 0, N, nmax)
    c3 = e_coeffs(0, (-d) % N, N, nmax)
    A = [c1[i]+c3[i] for i in range(nmax+1)]
    A = conv(A, c2, nmax)
    # W_{N^2}(e_{a,b}) = -(i/N) f_{a,b}; product of two: (-i/N)^2 f f = -1/N^2 f f
    f1 = f_coeffs(0, d, N, nmax); f2 = f_coeffs(b, 0, N, nmax)
    f3 = f_coeffs(0, (-d) % N, N, nmax)
    B = conv([f1[i]+f3[i] for i in range(nmax+1)], f2, nmax)
    B = [-x/N**2 for x in B]
    val = Lambda0(A, B, M)
    print("N=%d b=%d d=%d: value = %s" % (N, b, d, nstr(val, 20)))

print("="*70)
print("Part D: conductor-14 example, target 4 L'(14a4,0) = 0.90992489204940443157992925809053")
N = 14; M = N*N
nmax = 500
nu = {(0,5):1, (0,6):1, (0,1):-1, (0,2):-1}
nv = {(0,3):1, (0,5):1, (0,6):2, (0,1):-2, (0,2):-1, (0,4):-1}
symbols = [(-1, ((2,1),(7,4))), (-1, ((1,0),(4,1))), (1, ((1,0),(-4,1))), (1, ((-2,1),(7,-4)))]
# W-coefficient arrays are shared: build per transformed pair on the fly
total = mpc(0)
for sgn, al in symbols:
    ((r,s),(t,u_)) = al
    def act(ab):
        a, b = ab
        return ((a*r + b*t) % N, (a*s + b*u_) % N)
    usym = [(act(ab), n) for ab, n in nu.items()]
    vsym = [(act(ab), n) for ab, n in nv.items()]
    for (a, b), n1 in usym:
        for (c, d), n2 in vsym:
            # Thm 1: F = e_{a,d} e_{b,-c} + e_{a,-d} e_{b,c}
            A1 = conv(e_coeffs(a, d, N, nmax), e_coeffs(b, (-c) % N, N, nmax), nmax)
            A2 = conv(e_coeffs(a, (-d) % N, N, nmax), e_coeffs(b, c, N, nmax), nmax)
            A = [A1[i]+A2[i] for i in range(nmax+1)]
            B1 = conv(f_coeffs(a, d, N, nmax), f_coeffs(b, (-c) % N, N, nmax), nmax)
            B2 = conv(f_coeffs(a, (-d) % N, N, nmax), f_coeffs(b, c, N, nmax), nmax)
            B = [-(B1[i]+B2[i])/N**2 for i in range(nmax+1)]
            lam = Lambda0(A, B, M)
            total += sgn * n1 * n2 * lam
print("total =", nstr(total, 45))
print("target= 0.90992489204940443157992925809053394097")
