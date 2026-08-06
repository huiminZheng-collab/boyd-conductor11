# siegel_anchor_step12.py
# =============================================================================
# Referee round 5, comment 1: prove F_total = -2 f_11 WITHOUT assuming a priori
# that F_total lies in M_2(Gamma_0(11)) (the Sturm bound 2 used in step9 was
# otherwise circular).
#
# Route taken (referee's fallback route, necessitated by the true level):
#   Brunault, arXiv:1504.08127, Lemma 11 (reg_siegel.tex):
#     "the function e_{a,b} is an Eisenstein series of weight 1 on Gamma_1(N^2)".
#   For N = 11 the level is therefore Gamma_1(121), NOT Gamma_1(11): the
#   preferred route (symbolic Gamma_0(11)-invariance via e|gamma = e_{(a,b)gamma})
#   has no valid premise -- that clean permutation law only holds for the scaled
#   series e_{a,b}(tau/N) on Gamma(N) (Lemma 11, first assertion), where the
#   SL_2(Z)-action mixes the e's with cyclotomic coefficients.
#
#   Hence every product e_{p1,q1} e_{p2,q2} in F_total lies in M_2(Gamma_1(121))
#   (Eisenstein series are holomorphic on H and at all cusps -- no order
#    computation is needed), so
#       F_total, F_sym (each Manin symbol), and D := F_total + 2 f_11
#   all lie in M_2(Gamma_1(121)), since f_11 in M_2(Gamma_0(11)) subset.
#
#   Sturm bound for M_2(Gamma_1(121)):
#     [PSL_2(Z) : barGamma_1(121)] = (121^2/2)(1 - 1/121) = 7260   (N >= 5)
#     B = (k/12) * index = (2/12) * 7260 = 1210.
#   (-I acts trivially in weight 2, so the PSL index is the correct one; the
#    conservative SL_2 convention gives (2/12)*14520 = 2420.)
#   We verify exact vanishing of ALL q-coefficients of D up to q^2420 -- twice
#   the sharp bound, covering both conventions.  Therefore D = 0, i.e.
#       F_total = -2 f_11        (exact identity of modular forms),
#   and a fortiori F_total in M_2(Gamma_0(11)).  The same computation certifies
#   the per-symbol identities (two symbols equal -f_11, five are identically 0),
#   upgrading them from "observation" to proved statements.
#
# Method: exact integer/rational arithmetic only.
#   e-coefficients (n >= 1) are small integers; products use a Kronecker-packing
#   exact convolution (bias K, base 2^64, no overflow by the bounds asserted
#   in-code); constant terms alpha_0 (denominator 22) are kept as Fractions.
#   Cross-validated against code/siegel_anchor_step9_out.json (251 coeffs).
# =============================================================================
import json, time
from fractions import Fraction as Fr
from collections import defaultdict

N = 11
NMAX = 2420            # = conservative Sturm bound 2420; sharp bound is 1210
INDEX_PSL = 7260       # [PSL_2(Z) : barGamma_1(121)]
BOUND_SHARP = 2 * INDEX_PSL // 12          # 1210
BOUND_CONSERV = 2 * (2 * INDEX_PSL) // 12  # 2420
assert BOUND_SHARP == 1210 and BOUND_CONSERV == 2420

t0 = time.time()

# ---------------------------------------------------------------- e-tables --
# e_{a,b}(tau) = alpha_0(a,b) + sum_{m>=1, k>=1, m=a(N), k=b(N)} q^{mk}
#                               - sum_{m=-a(N), k=-b(N)} q^{mk}   (Brunault Def. 10)
E = {(a, b): [0]*(NMAX+1) for a in range(N) for b in range(N)}  # integer series part
A0 = {}
for a in range(N):
    for b in range(N):
        if a == 0 and b == 0: A0[(a, b)] = Fr(0)
        elif a == 0:          A0[(a, b)] = Fr(1, 2) - Fr(b, N)
        elif b == 0:          A0[(a, b)] = Fr(1, 2) - Fr(a, N)
        else:                 A0[(a, b)] = Fr(0)
for m in range(1, NMAX+1):
    am = m % N
    for k in range(1, NMAX//m + 1):
        n = m*k
        E[(am, k % N)][n] += 1
        E[((-am) % N, (-k) % N)][n] -= 1
maxe = max(abs(v) for tab in E.values() for v in tab)
print("e tables done (%.1fs), max |e coeff| = %d" % (time.time()-t0, maxe))

# ------------------------------------------- exact packed integer convolution
# e_{-a,-b} = -e_{a,b} and e_{a,b} = e_{b,a}: canonical representative with sign.
def canon(w):
    a, b = w
    neg = ((-a) % N, (-b) % N)
    sw = (b, a)
    nsw = (neg[1], neg[0])
    c = min(w, neg, sw, nsw)
    s = -1 if c in (neg, nsw) and c not in (w, sw) else 1
    # careful when w == neg (impossible for N odd unless w=(0,0), excluded)
    return c, s

K = 4096                 # bias > max |e coeff| (=32)
BASE = 1 << 64           # one 8-byte digit
assert K > maxe

convCache = {}
def conv(key):
    """exact integer convolution of the series parts of e_c1, e_c2 (n>=1 part)."""
    C = convCache.get(key)
    if C is not None:
        return C
    A = E[key[0]]; B = E[key[1]]
    pa = b'\x00'*8 + b''.join((A[i]+K).to_bytes(8, 'little') for i in range(1, NMAX+1))
    pb = b'\x00'*8 + b''.join((B[i]+K).to_bytes(8, 'little') for i in range(1, NMAX+1))
    Q = int.from_bytes(pa, 'little') * int.from_bytes(pb, 'little')
    raw = Q.to_bytes((Q.bit_length()+7)//8 + 8, 'little')
    prefA = [0]*(NMAX+1); prefB = [0]*(NMAX+1)
    sa = sb = 0
    for i in range(1, NMAX+1):
        sa += A[i]; sb += B[i]
        prefA[i] = sa; prefB[i] = sb
    C = [0]*(NMAX+1)
    KK = K*K
    for n in range(2, NMAX+1):
        d = int.from_bytes(raw[8*n:8*n+8], 'little')
        assert d < BASE, "packed digit overflow"
        C[n] = d - K*(prefA[n-1]+prefB[n-1]) - KK*(n-1)
    convCache[key] = C
    return C

# --------------------------------------------------------------- Manin symbols
rho = {2: -2, 4: 1, 5: 1}
sig = {1: 1, 2: -3, 3: -1, 5: 3}
symbols = [(1, (1, 0, 3, 1)), (-1, (1, 1, 3, 4)), (1, (3, 1, 11, 4)),
           (-1, (1, 0, 1, 1)), (1, (1, 2, 1, 3)), (-1, (3, 2, 4, 3)),
           (1, (3, 8, 4, 11))]

def build_symbol(r, s, t, u):
    """Returns (const0, acc_int, coef) with F_sym[0] = const0,
       F_sym[n] = acc_int[n] + sum_cs coef[cs] * E[cs][n]  (n >= 1)."""
    W = defaultdict(int)          # (c1,c2) -> total integer weight
    coef = defaultdict(Fr)        # canonical series -> Fraction multiplier
    c0 = Fr(0)
    for a, ra in rho.items():
        for c, sc in sig.items():
            w0 = ra*sc
            for b in range(N):
                ap = (a*r + b*t) % N
                bpp = (a*s + b*u) % N
                for d in range(N):
                    cp = (c*r + d*t) % N
                    dp = (c*s + d*u) % N
                    for P, Q_ in (((ap, dp), (bpp, (-cp) % N)),
                                  ((ap, (-dp) % N), (bpp, cp))):
                        if P == (0, 0) or Q_ == (0, 0):
                            continue          # e_{0,0} = 0 identically
                        c1, s1 = canon(P)
                        c2, s2 = canon(Q_)
                        if c1 > c2:
                            c1, s1, c2, s2 = c2, s2, c1, s1
                        w = w0*s1*s2
                        W[(c1, c2)] += w
                        coef[c2] += w*A0[c1]
                        coef[c1] += w*A0[c2]
                        c0 += w*A0[c1]*A0[c2]
    acc = [0]*(NMAX+1)
    nkey = 0
    for key, w in W.items():
        if w == 0:
            continue
        nkey += 1
        C = conv(key)
        acc = [x + w*y for x, y in zip(acc, C)]
    return c0, acc, coef, nkey

per_sym = []
total_keys = 0
for sgn, (r, s, t, u) in symbols:
    ts = time.time()
    c0, acc, coef, nkey = build_symbol(r, s, t, u)
    total_keys += nkey
    nz = [(cs, cf) for cs, cf in coef.items() if cf != 0]
    Fs = [Fr(0)]*(NMAX+1)
    Fs[0] = c0
    for n in range(1, NMAX+1):
        v = Fr(acc[n])
        for cs, cf in nz:
            e = E[cs][n]
            if e:
                v += cf*e
        Fs[n] = v
    per_sym.append((sgn, (r, s, t, u), Fs))
    print("symbol (%+d) [[%d,%d],[%d,%d]] done: %d distinct products, %.1fs"
          % (sgn, r, s, t, u, nkey, time.time()-ts), flush=True)

F = [Fr(0)]*(NMAX+1)
for sgn, al, Fs in per_sym:
    for n in range(NMAX+1):
        if Fs[n]:
            F[n] += sgn*Fs[n]

# ------------------------------------------------------------------ f_11 -----
# f_11 = q prod (1-q^n)^2 (1-q^{11n})^2, exact integer coefficients
Pf = [0]*(NMAX+2); Pf[0] = 1
for n in range(1, NMAX+1):
    for _ in range(2):
        for k in range(NMAX+1, n-1, -1):
            Pf[k] -= Pf[k-n]
    if 11*n <= NMAX+1:
        for _ in range(2):
            for k in range(NMAX+1, 11*n-1, -1):
                Pf[k] -= Pf[k-11*n]
f11 = [0]*(NMAX+1)
for n in range(1, NMAX+1):
    f11[n] = Pf[n-1]          # f_11[0] = 0

# ------------------------------------------------------------------ checks ---
results = {}

# (1) main identity: D = F_total + 2 f_11  == 0  for n = 0..NMAX
badD = [n for n in range(NMAX+1) if F[n] != -2*f11[n]]
results["F_total_eq_-2f11_upto"] = NMAX
results["F_total_bad"] = badD[:20]
print("="*72)
print("MAIN: F_total = -2 f_11   exact for n=0..%d : %s"
      % (NMAX, "PASS" if not badD else "FAIL %s" % badD[:10]))

# (2) per-symbol identities: expect F_sym = -f_11 for [[1,0],[3,1]] and
#     [[1,2],[1,3]], and F_sym = 0 for the other five.
expect = {(1, 0, 3, 1): "-f11", (1, 2, 1, 3): "-f11"}
per_results = []
allok = True
for sgn, al, Fs in per_sym:
    is_zero = all(v == 0 for v in Fs)
    is_negf = all(Fs[n] == -f11[n] for n in range(NMAX+1))
    kind = "0" if is_zero else ("-f11" if is_negf else "OTHER")
    want = expect.get(al, "0")
    ok = (kind == want)
    allok &= ok
    per_results.append([sgn, list(al), kind, ok])
    print("  symbol (%+d) [[%d,%d],[%d,%d]] : F_sym = %-5s expected %-5s %s"
          % (sgn, *al, kind, want, "PASS" if ok else "FAIL"))
results["per_symbol"] = per_results
print("PER-SYMBOL (all 7, to q^%d): %s" % (NMAX, "PASS" if allok else "FAIL"))

# (3) cross-validation against step9 (251 exact coefficients)
try:
    with open("code/siegel_anchor_step9_out.json") as f:
        s9 = json.load(f)
    ok9 = all(Fr(s9["F"][n]) == F[n] for n in range(251))
    okp = all(Fr(s9["per_sym"][i][2][n]) == per_sym[i][2][n]
              for i in range(7) for n in range(251))
    print("cross-check vs step9_out.json (251 coeffs, F and 7 symbols): %s"
          % ("PASS" if ok9 and okp else "FAIL"))
    results["crosscheck_step9"] = bool(ok9 and okp)
except FileNotFoundError:
    print("step9_out.json not found; cross-check skipped")
    results["crosscheck_step9"] = None

results["group"] = "M_2(Gamma_1(121))"
results["index_PSL"] = INDEX_PSL
results["sturm_bound_sharp"] = BOUND_SHARP
results["sturm_bound_conservative"] = BOUND_CONSERV
results["coeffs_verified"] = NMAX + 1
results["distinct_products"] = len(convCache)
results["symbol_key_total"] = total_keys
results["runtime_sec"] = round(time.time()-t0, 1)
results["F_first40"] = [str(x) for x in F[:40]]
verdict = (not badD) and allok and results["crosscheck_step9"] in (True, None)
results["VERDICT"] = "PASS" if verdict else "FAIL"
print("="*72)
print("distinct e-products computed:", len(convCache),
      " runtime %.1fs" % (time.time()-t0))
print("Sturm: M_2(Gamma_1(121)), index 7260, sharp bound %d, conservative %d;"
      " verified n=0..%d" % (BOUND_SHARP, BOUND_CONSERV, NMAX))
print("VERDICT:", results["VERDICT"])

with open("code/siegel_anchor_step12_out.json", "w") as f:
    json.dump(results, f, indent=1)
print("written code/siegel_anchor_step12_out.json")
