"""
Sanity check for the closedness argument (Samart Lemma 9 mirrored to S_0).

On E: S_0 = y^2+(x^2+1)y+x^3 = 0, take the CONTINUOUS big-modulus branch
y_big(e^{iθ}) (|y_big| >= 1 everywhere; at toric points |y|=1 the branch is
continuous and may switch algebraic roots). The signed cycle on the FULL circle
    gamma~ = (+1) * arc(-c -> c) + (-1) * arc(c -> 2pi - c),   c = pi/2
should satisfy:
  (i)  endpoint matching: lim_{θ->±pi} y_big = same point (-1, -1-sqrt2)
       => ∂gamma~ = 2(P_c - P_{-c}) = 2(P_c - conj P_c), anti-invariant;
  (ii) the signed integral equals -2 * I_split (by conj symmetry),
       hence = -2 b_11 (using the 149-digit result).
"""
from mpmath import mp, mpf, mpc, pi, sqrt, exp, fabs, log, nstr
import sys
sys.path.insert(0, ".")
from b11 import b11

mp.dps = 60

def roots(x):
    d = sqrt((x**2+1)**2 - 4*x**3)
    return (-(x**2+1) + d)/2, (-(x**2+1) - d)/2

def ybig(th):
    a, b = roots(exp(mpc(0,1)*th))
    return a if fabs(a) >= fabs(b) else b

# (i) endpoint matching at θ = ±π
for side, th in [("pi-", pi - mpf('1e-12')), ("-pi+", -pi + mpf('1e-12'))]:
    print(f"  y_big({side}) =", nstr(ybig(th), 30))
print("  -1 - sqrt(2) =", nstr(-1-sqrt(2), 30))

# (ii) full-circle signed integral vs -2 b_11
c = pi/2
W = lambda th: 1 if fabs(th) < c else -1
f = lambda th: W(th) * log(fabs(ybig(th)))
val = (mp.quad(f, [-pi, -c]) + mp.quad(f, [-c, c]) + mp.quad(f, [c, pi])) / (2*pi)
print("\n  (1/2π) ∫_{full,signed} log|y_big| =", nstr(val, 50))
b = b11(60)
print("  -b_11/π × 2 × π/2 ... direct compare: -2*b_11/(2π)*π = ")
print("  -b_11 =", nstr(-b, 50))
print("  |val - (-b_11)| =", nstr(fabs(val + b), 8))

# interpretation: val should equal -b_11 (since (1/2π)*signed_full = (1/π)(I2-I1)*(1/1)... )
# Actually by conj symmetry signed_full = 2(I2'-I1') where Ii' = integrals of log|y_big| = -Ii,
# so signed_full = -2(I2-I1) = 2 I_split = 2 b_11, and (1/2π)*that = b_11/π * ... see report.
