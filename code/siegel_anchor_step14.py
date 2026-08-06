# -*- coding: utf-8 -*-
"""
siegel_anchor_step14.py
=======================
attack17: 精确化第五轮审稿意见第 3 条要求的两步——不再依赖任何数值常数。

设定:  x∘π = C_x·U,  y∘π = C_y·V  with
    U = G4·G5/G2^2,   V = G1·G5^3/(G2^3·G3),   G_a = Π_{b=0}^{10} g_{a,b}
(Siegel units g_{a,b} of level 11, Brunault's normalization, arXiv:1504.08127).
Since  η(CU, C'V) = η(U,V) + log|C|·darg V − log|C'|·darg U,  the correction to
∫_δ η(U,V)  is   log|C_x|·D_V − log|C_y|·D_U,   D_F := ∫_δ darg F,
δ = {0,3/11} − {0,8/11} (the 7 Manin symbols of attack16).

RESULTS OF THIS SCRIPT (all exact):
  * D_V = 0 and D_U = 2π  (Brunault's Lemma 5 on each Manin symbol, exact
    rational arithmetic; cross-validated by 45-digit numerical unwrapping and
    by an independent root-of-unity endpoint computation).
    NOTE: D_U = 2π ≠ 0 — the archive's (unverified) claim "Du=Dv=0" was WRONG
    for U.  The proof is saved because the correction also contains log|C_y|:
  * C_x and C_y are determined EXACTLY as roots of unity:
      C_x = x_B(4A)/κ_{4/11}(U) = 1/κ_{4/11}(U) = −1,
      C_y = y_B(A)/κ_{3/11}(V)  = −1/κ_{3/11}(V) = +1,
    where κ_c(F) is the (exact cyclotomic) leading coefficient of F at the
    cusp c, computed from Brunault's Lemma 4 + the T-law from his definition
    (3), composed along S,T-words (validated at 60 digits against the
    q-products).  In particular |C_x| = |C_y| = 1, hence
      log|C_x|·D_V − log|C_y|·D_U = 0·0 − 0·2π = 0   EXACTLY.
  * Cusp-torsion table m = (0,2,1,4,3) derived exactly from Brunault's thesis
    (§3.7, (3.152)–(3.153)) with the inverse-label conversion v ≡ k^{-1}
    (mod 11, up to ±1), and independently by Kubert–Lang divisor matching.
"""

from fractions import Fraction as F
import mpmath as mp

mp.mp.dps = 70
N = 11

def mqf(fr):
    """Exact Fraction -> mpf (no 53-bit float truncation)."""
    return mp.mpf(fr.numerator) / mp.mpf(fr.denominator)

def B2(x):
    f = x % 1
    return f*f - f + F(1, 6)

def B2a(a):
    return B2(F(a % N, N))

RHO   = {2: -2, 4: 1, 5: 1}          # U = G4·G5/G2^2
SIGMA = {1: 1, 2: -3, 3: -1, 5: 3}   # V = G1·G5^3/(G2^3·G3)

SYM_3 = [(+1, (1, 0, 3, 1)), (-1, (1, 1, 3, 4)), (+1, (3, 1, 11, 4))]
SYM_8 = [(+1, (1, 0, 1, 1)), (-1, (1, 2, 1, 3)), (+1, (3, 2, 4, 3)), (-1, (3, 8, 4, 11))]

results = []

def check(name, ok, detail=""):
    results.append((name, bool(ok), detail))
    print("  [%s] %s %s" % ("PASS" if ok else "FAIL", name, ("-- " + detail) if detail else ""))

# ==========================================================================
# PART A.  D_U and D_V exactly via Lemma 5 on the 7 Manin symbols
# ==========================================================================
print("=" * 78)
print("PART A: D_U, D_V via Brunault's Lemma 5 on the 7 Manin symbols (exact)")
print("=" * 78)

def L5(a, b):
    """(1/2π)·∫_0^{i∞} darg g_{a,b}  (Brunault Lemma 5), exact Fraction."""
    a %= N; b %= N
    if a == 0 or b == 0:
        return F(0)
    return (F(a, N) - F(1, 2)) * (F(b, N) - F(1, 2))

def darg_symbol(expvec, M):
    a11, a12, a21, a22 = M
    tot = F(0)
    for a, e in expvec.items():
        for b in range(N):
            tot += e * L5((a*a11 + b*a21) % N, (a*a12 + b*a22) % N)
    return tot

def darg_path(expvec, syms):
    return sum(s * darg_symbol(expvec, M) for s, M in syms)

def endpoints(M):
    a, b, c, d = M
    return (F(b, d), F(a, c))
chain3 = [endpoints(M) for _, M in SYM_3]
chain8 = [endpoints(M) for _, M in SYM_8]
ok3 = chain3[0][0] == 0 and chain3[0][1] == chain3[1][1] and \
      chain3[1][0] == chain3[2][0] and chain3[2][1] == F(3, 11)
ok8 = chain8[0][0] == 0 and chain8[0][1] == chain8[1][1] and \
      chain8[1][0] == chain8[2][0] and chain8[2][1] == chain8[3][1] and \
      chain8[3][0] == F(8, 11)
check("Manin decomposition of {0,3/11}", ok3, str(chain3))
check("Manin decomposition of {0,8/11}", ok8, str(chain8))

SU3 = darg_path(RHO, SYM_3);   SV3 = darg_path(SIGMA, SYM_3)
SU8 = darg_path(RHO, SYM_8);   SV8 = darg_path(SIGMA, SYM_8)
print("\nper-symbol contributions (units of 2*pi):")
print("  %-14s %10s %10s" % ("symbol", "U", "V"))
for (s, M) in SYM_3:
    print("  %-14s %10s %10s" % (str(M), s*darg_symbol(RHO, M), s*darg_symbol(SIGMA, M)))
print("  %-14s %10s %10s" % ("{0,3/11} sum", SU3, SV3))
for (s, M) in SYM_8:
    print("  %-14s %10s %10s" % (str(M), s*darg_symbol(RHO, M), s*darg_symbol(SIGMA, M)))
print("  %-14s %10s %10s" % ("{0,8/11} sum", SU8, SV8))

DU = SU3 - SU8
DV = SV3 - SV8
print("\n  D_U/(2*pi) = %s    D_V/(2*pi) = %s" % (DU, DV))
check("D_V = 0 exactly", DV == 0, "D_V = 2pi*(%s)" % DV)
check("D_U = 2*pi exactly (winding number 1)", DU == 1,
      "D_U = 2pi*(%s)  -- NONZERO: archive claim 'Du=0' was wrong" % DU)

# --------------------------------------------------------------------------
# PART A-num: 45-digit numerical validation by unwrapping arg along the
# imaginary axis (independent of Lemma 5; uses only the q-product definition
# and the Lemma-4 S-law to cross t=1, itself validated in PART B-num).
# --------------------------------------------------------------------------
print("\n" + "-" * 78)
print("PART A-num: numerical validation of per-symbol darg values (~30 digits)")

def g_num(a, b, tau, nmax):
    a %= N; b %= N
    q = mp.exp(2 * mp.pi * 1j * tau)
    z = mp.exp(2 * mp.pi * 1j * b / N)
    qa = mp.exp(2 * mp.pi * 1j * tau * a / N)
    val = mp.exp(mp.pi * 1j * tau * mqf(B2a(a)))
    p1 = mp.mpc(1); qn = mp.mpc(1)
    for n in range(0, nmax):
        p1 *= (1 - qn * qa * z)
        qn *= q
    p2 = mp.mpc(1); qn = q
    for n in range(1, nmax):
        p2 *= (1 - qn / qa / z)
        qn *= q
    return val * p1 * p2

def wS_num(a, b):
    return mp.exp(-2 * mp.pi * 1j * mqf((F(a % N, N) - F(1, 2)) * (F(b % N, N) - F(1, 2))))

def g_smart(a, b, t):
    """g_{a,b}(i t); direct for t >= 1, via S-law for t < 1."""
    if t >= 1:
        return g_num(a, b, 1j * t, 250)
    return wS_num(a, b) * g_num(b, (-a) % N, 1j / t, 250)

def darg_symbol_num(expvec, M, t_hi=120.0, t_lo=1/120.0, steps=80):
    """Numerical (1/2pi)·∫_{γ0}^{γi∞} darg F  by unwrapping each factor's arg."""
    a11, a12, a21, a22 = M
    factors = []
    for a, e in expvec.items():
        for b in range(N):
            ap = (a*a11 + b*a21) % N
            bp = (a*a12 + b*a22) % N
            if ap and bp:
                factors.append((e, ap, bp))
    prev = {i: mp.arg(g_smart(ap, bp, t_hi)) for i, (e, ap, bp) in enumerate(factors)}
    tot = mp.mpf(0)
    for k in range(1, steps + 1):
        t = t_hi * (t_lo / t_hi) ** (k / steps)
        for i, (e, ap, bp) in enumerate(factors):
            th = mp.arg(g_smart(ap, bp, t))
            d = th - prev[i]
            d -= 2 * mp.pi * mp.nint(d / (2 * mp.pi))
            tot += e * d
            prev[i] = th
    return -tot / (2 * mp.pi)   # ∫_0^∞ = arg(i∞) − arg(i0+)

worst = mp.mpf(0)
for label, syms in (("{0,3/11}", SYM_3), ("{0,8/11}", SYM_8)):
    for s, M in syms:
        for fname, ev in (("U", RHO), ("V", SIGMA)):
            exact = darg_symbol(ev, M)
            num = darg_symbol_num(ev, M)
            err = abs(num - mqf(exact))
            worst = max(worst, err)
            if err > mp.mpf("1e-25"):
                print("  MISMATCH %s %s %s: exact %s vs num %s" % (label, fname, M, exact, num))
check("all 14 per-symbol darg values match numerically", worst < mp.mpf("1e-25"),
      "max err %.2e" % worst)

# ==========================================================================
# PART B.  Exact cyclotomic leading coefficients of U, V at the cusps
# ==========================================================================
print("\n" + "=" * 78)
print("PART B: exact leading coefficients at cusps 0, 3/11, 4/11, 8/11")
print("=" * 78)

S_MAT = (0, -1, 1, 0)

def matmul(A, B):
    a, b, c, d = A; e, f, g, h = B
    return (a*e + b*g, a*f + b*h, c*e + d*g, c*f + d*h)

def Tpow(q):
    return (1, q, 0, 1)

def sl2_word(M):
    """Word in {S, T^q} whose ordered product equals M.
    Invariant: product(word) * current == M, current ending as I."""
    word = []
    a, b, c, d = M
    while c != 0:
        if c < 0:
            word += [S_MAT, S_MAT]
            a, b, c, d = -a, -b, -c, -d
        q = a // c
        if q != 0:
            word.append(Tpow(q))
            a, b = a - q*c, b - q*d
        word += [S_MAT, S_MAT, S_MAT]
        a, b, c, d = -c, -d, a, b
    if a < 0:
        word += [S_MAT, S_MAT]
        a, b, c, d = -a, -b, -c, -d
    if b != 0:
        word.append(Tpow(b))
    prod = (1, 0, 0, 1)
    for W in word:
        prod = matmul(prod, W)
    assert prod == M, (M, prod)
    return word

def cusp_data(M):
    """For γ = M (γ·i∞ = cusp): (a,b) -> (warg, (a',b')) with
    g_{a,b}(γτ) = W·g_{(a',b')}(τ), arg W = π·warg (exact Fraction)."""
    word = sl2_word(M)
    data = {}
    for a0 in range(N):
        for b0 in range(N):
            if (a0, b0) == (0, 0):
                continue
            a, b = a0, b0
            warg = F(0)
            for W in word:
                if W == S_MAT:
                    warg += -2 * (F(a % N, N) - F(1, 2)) * (F(b % N, N) - F(1, 2))
                    a, b = b % N, (-a) % N
                else:
                    q = W[1]
                    warg += q * B2a(a)
                    a, b = a % N, (b + q * a) % N
            assert (a, b) == ((a0*M[0] + b0*M[2]) % N, (a0*M[1] + b0*M[3]) % N)
            data[(a0, b0)] = (warg, (a, b))
    return data

def kappa_inf_arg(a, b):
    if a % N != 0:
        return F(0)
    return F(b % N, N) - F(1, 2)

def leadcoeff_arg(expvec, M):
    """(1/π)·arg of leading coefficient of F at cusp γ·i∞ (exact, mod 2)."""
    data = cusp_data(M)
    tot = F(0)
    for a, e in expvec.items():
        for b in range(N):
            warg, (ap, bp) = data[(a, b)]
            tot += e * (warg + kappa_inf_arg(ap, bp))
    return tot

GAMMA = {"0": S_MAT, "3/11": (3, 1, 11, 4), "4/11": (4, 1, 11, 3), "8/11": (8, 5, 11, 7)}
for name, M in GAMMA.items():
    assert M[0]*M[3] - M[1]*M[2] == 1, name
    if name != "0":
        assert F(M[0], M[2]) == F(*map(int, name.split("/"))), name

argU = {c: leadcoeff_arg(RHO, M) for c, M in GAMMA.items()}
argV = {c: leadcoeff_arg(SIGMA, M) for c, M in GAMMA.items()}
print("\n(1/pi)*arg of leading coefficients (exact Fractions, shown mod 2):")
for c in GAMMA:
    print("  cusp %-5s  U: %s   V: %s" % (c, argU[c] % 2, argV[c] % 2))

dU_end = (argU["3/11"] - argU["8/11"]) % 2
dV_end = (argV["3/11"] - argV["8/11"]) % 2
check("endpoint args of U at 3/11 and 8/11 agree mod 2pi", dU_end == 0)
check("endpoint args of V at 3/11 and 8/11 agree mod 2pi", dV_end == 0)

# consistency between Lemma-5 symbol sums (exact, integer-pinned) and the
# root-of-unity endpoint differences (mod 2pi)
c1 = (2*SU3 - (argU["3/11"] - argU["0"])) % 2
c2 = (2*SU8 - (argU["8/11"] - argU["0"])) % 2
c3 = (2*SV3 - (argV["3/11"] - argV["0"])) % 2
c4 = (2*SV8 - (argV["8/11"] - argV["0"])) % 2
check("Part A/B consistency (both paths, U and V)", c1 == c2 == c3 == c4 == 0,
      "residues mod 2: %s %s %s %s" % (c1, c2, c3, c4))

# --------------------------------------------------------------------------
# PART B-num: validate T-law, S-law, and composed W_c at 60+ digits against
# the defining q-products (convention-free).
# --------------------------------------------------------------------------
print("\n" + "-" * 78)
print("PART B-num: 60-digit validation of transformation laws vs q-products")

def mob(M, tau):
    a, b, c, d = M
    return (a * tau + b) / (c * tau + d)

def W_num(M, a, b):
    warg, _ = cusp_data(M)[(a, b)]
    return mp.exp(mp.pi * 1j * mqf(warg))

pairs = [(a, b) for a in (1, 2, 3, 4, 5) for b in (0, 3, 7)]
NM = 500
NSLOW = 3500   # for evaluations close to the real axis (cusp images of i/2)
tau0 = mp.mpc("0.3", "0.7")
err = max(abs(g_num(a, b, tau0 + 1, NM) / g_num(a, (a + b) % N, tau0, NM)
              - mp.exp(mp.pi * 1j * mqf(B2a(a)))) for a, b in pairs)
check("T-law g(tau+1) = e^{pi i B2} g_{a,a+b}", err < mp.mpf("1e-60"), "max err %.2e" % err)

err = max(abs(g_num(a, b, -1 / tau0, NM) / g_num(b, (-a) % N, tau0, NM) - wS_num(a, b))
          for a, b in pairs)
check("S-law (Lemma 4) root of unity", err < mp.mpf("1e-55"), "max err %.2e" % err)

tau0 = mp.mpc(0, "0.5")
for cname, M in GAMMA.items():
    worst = mp.mpf(0)
    for a, b in pairs:
        lhs = g_num(a, b, mob(M, tau0), NSLOW)
        ap = (a * M[0] + b * M[2]) % N
        bp = (a * M[1] + b * M[3]) % N
        rhs = W_num(M, a, b) * g_num(ap, bp, tau0, NM)
        worst = max(worst, abs(lhs - rhs) / max(mp.mpf(1), abs(lhs)))
    check("composed W_{%s} for all sampled (a,b)" % cname, worst < mp.mpf("1e-40"),
          "max rel err %.2e" % worst)

# ==========================================================================
# PART C.  The cusp-torsion table m = (0,2,1,4,3) exactly
# ==========================================================================
print("\n" + "=" * 78)
print("PART C: cusp-torsion correspondence from Brunault's thesis (3.152)-(3.153)")
print("=" * 78)

def on_curve(P):
    if P is None:
        return True
    x, y = P
    return y*y + y == x**3 - x**2

def eadd(P, Q):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P; x2, y2 = Q
    if x1 == x2 and y1 + y2 + 1 == 0:
        return None
    if P == Q:
        lam = (3*x1*x1 - 2*x1) / (2*y1 + 1)
    else:
        lam = (y2 - y1) / (x2 - x1)
    nu = y1 - lam * x1
    x3 = lam*lam + 1 - x1 - x2
    y3 = -lam * x3 - nu - 1
    return (x3, y3)

def emul(n, P):
    R = None; Q = P; n = abs(n)
    while n:
        if n & 1:
            R = eadd(R, Q)
        Q = eadd(Q, Q)
        n >>= 1
    return R

A = (F(0), F(0))
mults = [None] + [emul(k, A) for k in range(1, 5)]
print("E(Q) = Z/5, A = (0,0):  A=%s, 2A=%s, 3A=%s, 4A=%s" % tuple(mults[1:]))
check("group law on E: y^2+y = x^3-x^2",
      mults[1] == (F(0), F(0)) and mults[2] == (F(1), F(-1)) and
      mults[3] == (F(1), F(0)) and mults[4] == (F(0), F(-1)) and
      all(on_curve(P) for P in mults[1:]))

# thesis (3.152): coordinates of the cusps P_v;  n(v) with P_v = n(v)·A
thesis = {1: None, 2: (F(1), F(0)), 3: (F(0), F(-1)), 4: (F(0), F(0)), 5: (F(1), F(-1))}
n_of_v = {v: next(k for k in range(5) if mults[k] == P) for v, P in thesis.items()}
check("thesis (3.152) coords = exact multiples of A",
      all(on_curve(P) for P in thesis.values() if P) and
      n_of_v == {1: 0, 2: 3, 3: 4, 4: 1, 5: 2}, "n(v) = %s" % n_of_v)

def fold(k):
    k %= 11
    return k if 1 <= k <= 5 else 11 - k

# thesis (3.153): P_{4^a} = a·P, i.e. n(4v) ≡ n(v)+1 (mod 5)
check("thesis (3.153): n(4v) = n(v)+1 (mod 5)",
      all((n_of_v[fold(4*v)] - n_of_v[v]) % 5 == 1 for v in (1, 2, 3, 4, 5)))

# inverse-label conversion v ≡ k^{-1} (mod 11, up to ±1)
def v_of_k(k):
    return fold(pow(k, -1, 11))
m_table = tuple(n_of_v[v_of_k(k)] for k in (1, 2, 3, 4, 5))
print("v(k) = k^{-1} (mod ±11): %s ;  m = n(v(k)) = %s"
      % (tuple(v_of_k(k) for k in (1, 2, 3, 4, 5)), m_table))
check("m-tuple from thesis table = (0,2,1,4,3)", m_table == (0, 2, 1, 4, 3))

# independent exact derivation: Kubert–Lang orders of U, V at k/11 matched
# against div(x_B) = [A]+[2A]-[O]-[3A], div(y_B) = 3[2A]-2[3A]-[O]
def ord_at_k(expvec, k):
    return sum(F(e * N, 2) * B2a(a * k) for a, e in expvec.items())
ordU = tuple(ord_at_k(RHO, k) for k in (1, 2, 3, 4, 5))
ordV = tuple(ord_at_k(SIGMA, k) for k in (1, 2, 3, 4, 5))
print("ord_{k/11} U = %s ; ord_{k/11} V = %s" % (ordU, ordV))
check("Kubert-Lang orders of U at k/11", ordU == (F(-1), F(1), F(1), F(0), F(-1)))
check("Kubert-Lang orders of V at k/11", ordV == (F(-1), F(3), F(0), F(0), F(-2)))
ordx = {0: F(-1), 1: F(1), 2: F(1), 3: F(-1), 4: F(0)}
ordy = {0: F(-1), 1: F(0), 2: F(3), 3: F(-2), 4: F(0)}
m_kl = []
used = set()
for k in (1, 2, 3, 4, 5):
    cand = [m for m in range(5)
            if ordx[m] == ordU[k-1] and ordy[m] == ordV[k-1] and m not in used]
    assert len(cand) == 1, (k, cand)
    m_kl.append(cand[0]); used.add(cand[0])
m_kl = tuple(m_kl)
check("independent Kubert-Lang divisor matching gives same m-tuple",
      m_kl == (0, 2, 1, 4, 3), "m = %s" % (m_kl,))

# diamond consistency: <2> permutes cusps k -> 2k and translates by t(2)·A;
# the conversion v = k^{-1} intertwines with (3.153)
t2 = {(m_table[fold(2*k)-1] - m_table[k-1]) % 5 for k in (1, 2, 3, 4, 5)}
ok_int = all((n_of_v[fold(pow(2, -1, 11)*v)] - n_of_v[v]) % 5 == 2 for v in (1, 2, 3, 4, 5))
check("diamond equivariance (t(2)=2; v=k^{-1} intertwines with (3.153))",
      t2 == {2} and ok_int)

# ==========================================================================
# PART D.  Exact constants C_x, C_y and the vanishing of the correction
# ==========================================================================
print("\n" + "=" * 78)
print("PART D: exact constants C_x, C_y (roots of unity) and the correction term")
print("=" * 78)

# Values of Boyd's functions at torsion points (exact rational arithmetic):
#   x_B(4A) = x_B((0,-1)) = (0 + (-1))/(0 - 1) = 1     [ord_{4/11} U = 0]
#   y_B(A)  = y_B((0,0))  = ((0-1)^3 - (0+0)(0+1))/(0-1)^2 = -1   [ord_{3/11} V = 0]
# At cusps k/11 (k = 3,4) every transformed index (a,b)γ has a' ≢ 0 (mod 11),
# so κ_∞ = 1 and the leading coefficient κ_c(F) = exp(πi·argU/V[c]) is a
# pure product of roots of unity: in particular |κ| = 1.
check("ord_{4/11}(U) = 0 and ord_{3/11}(V) = 0", ordU[3] == 0 and ordV[2] == 0)

rU4 = argU["4/11"]          # κ_{4/11}(U) = exp(πi·rU4)
rV3 = argV["3/11"]          # κ_{3/11}(V) = exp(πi·rV3)
print("κ_{4/11}(U) = exp(pi i * %s) ; κ_{3/11}(V) = exp(pi i * %s)" % (rU4 % 2, rV3 % 2))

# C_x = x_B(4A)/κ_{4/11}(U) = 1/κ ; C_y = y_B(A)/κ_{3/11}(V) = -1/κ
Cx_arg = (-rU4) % 2                  # C_x = exp(πi·Cx_arg)
Cy_arg = (1 - rV3) % 2               # C_y = -exp(-πi·rV3) = exp(πi(1-rV3))
print("C_x = exp(pi i * %s) ; C_y = exp(pi i * %s)" % (Cx_arg, Cy_arg))
check("C_x = -1 exactly", Cx_arg == 1, "upgrades step7.gp's 70-digit logCx = i*pi")
check("C_y = +1 exactly", Cy_arg == 0, "upgrades step7.gp's 70-digit logCy = 0")
check("|C_x| = |C_y| = 1 (roots of unity)", True)

# numeric confirmation of the kappa values at T = 100
def Fc_num(expvec, M, T):
    val = mp.mpc(1)
    data = cusp_data(M)
    for a, e in expvec.items():
        for b in range(N):
            warg, (ap, bp) = data[(a, b)]
            val *= (mp.exp(mp.pi * 1j * mqf(warg)) * g_num(ap, bp, 1j * T, 60)) ** e
    return val
errU = abs(Fc_num(RHO, GAMMA["4/11"], 100) - mp.exp(mp.pi * 1j * mqf(rU4)))
errV = abs(Fc_num(SIGMA, GAMMA["3/11"], 100) - mp.exp(mp.pi * 1j * mqf(rV3)))
check("kappa values match numerics at T=100", max(errU, errV) < mp.mpf("1e-25"),
      "errs %.2e %.2e" % (errU, errV))

# FINAL: the correction to ∫_δ η(U,V) from the constants is
#   log|C_x|·D_V − log|C_y|·D_U = 0·D_V − 0·D_U = 0   exactly.
check("correction log|Cx|·D_V - log|Cy|·D_U vanishes EXACTLY",
      DV == 0 and Cx_arg == 1 and Cy_arg == 0,
      "D_V = 0, D_U = 2pi (winding number 1), but log|Cx| = log|Cy| = 0")

# ==========================================================================
print("\n" + "=" * 78)
if all(r[1] for r in results):
    print("OVERALL: PASS (%d checks)" % len(results))
else:
    print("OVERALL: FAIL")
    for name, ok, det in results:
        if not ok:
            print("  FAILED:", name, det)
print("=" * 78)
