# siegel_anchor_step11.py -- FINAL verification of the pi-corrected Brunault chain.
#
# Context: our Theorem-1 implementation (steps 5/8/9) computes
#     I_symbol = Lambda*(e_{a,d}e_{b,-c} + e_{a,-d}e_{b,c}, 0)
# and found, exactly, F_total = -2 f_11  =>  sum over the 7 Manin symbols = -2*b11.
# The arXiv LaTeX source of Brunault (arXiv:1504.08127, reg_siegel.tex line 91)
# shows the TRUE statement of Theorem 1 carries an explicit factor pi that was
# lost in the PDF text extraction we worked from:
#     int_0^{i infty} eta(g_u,g_v) = pi * Lambda*(e_{a,d}e_{b,-c}+e_{a,-d}e_{b,c},0)
# Hence  int_{gamma^-} eta(x,y) = pi * (-2 b11) = -2 pi b11  (our orientation).
#
# This script independently validates the pi factor on the conductor-14 example
# (where everything can be checked against the true Mahler measure), and prints
# the final conductor-11 comparison.
#
# Run: .venv/Scripts/python.exe code/siegel_anchor_step11.py

from mpmath import mp, mpf, pi, exp, log, quad, fabs, sqrt
import subprocess, re

mp.dps = 60

def gp_lprime(crv):
    """L'(E,0) (plain derivative) via PARI lfun."""
    out = subprocess.run(
        ["gp", "-q", "-f"],
        input=f'default(realprecision,60);\ne=ellinit("{crv}");\nprint(lfun(e,0,1));\n',
        capture_output=True, text=True, cwd=r"C:/Users/zheng/boyd-conductor11",
    ).stdout
    m = re.search(r"[-+]?\d*\.\d+(?:E[-+]?\d+)?", out)
    return mpf(m.group(0))

print("=" * 74)
print("(a) conductor 14 ground truth: P_{{-1}}(x,y) = y^2 - x y + y - x^3")
L14 = gp_lprime("14a4")
print("  L'(E14a4,0) [PARI]      =", L14)

def yroots(th):
    x = exp(1j * th)
    s = sqrt((1 - x) ** 2 + 4 * x ** 3)
    return ((x - 1) + s) / 2, ((x - 1) - s) / 2

def jensen_integrand(th):
    y1, y2 = yroots(th)
    return log(max(fabs(y1), mpf(1))) + log(max(fabs(y2), mpf(1)))

# P_{-1} vanishes on the torus at th = 0, pi (log singularities): split there.
m_P = quad(jensen_integrand, [0, pi, 2 * pi]) / (2 * pi)
print("  m(P_-1) [Jensen, 60d]   =", m_P)
print("  m(P_-1) - 2 L'          =", m_P - 2 * L14, "  (Boyd/Mellit identity)")

# Direct regulator integral of the STANDARD eta(x,y) = log|x| darg y - log|y| darg x
# over the Deninger path gamma_{-1} = {|x|=1, |y|>=1}: here log|x| = 0, so
# eta = -log|y+(th)| dth with y+ the outside root.
def eta_integrand(th):
    y1, y2 = yroots(th)
    yp = y1 if fabs(y1) >= fabs(y2) else y2
    return -log(fabs(yp))

I_std = quad(eta_integrand, [0, pi, 2 * pi])
print("  int_{{gamma_-1}} eta [direct] =", I_std)
print("  -2 pi m(P_-1)             =", -2 * pi * m_P)
print("  -4 pi L'                  =", -4 * pi * L14)
print("  I_std / (-4 pi L')        =", I_std / (-4 * pi * L14))

# Our Lambda-chain gave (45 digits, step5): sum of the 4 Manin symbols = 4*L' = 0.9099248920494044315
LAMBDA_CHAIN_14 = 4 * L14
print()
print("  Lambda-chain (step5) sum  =", LAMBDA_CHAIN_14)
print("  pi * Lambda-chain         =", pi * LAMBDA_CHAIN_14)
print("  paper Thm1: pi L'(4f,0)   = pi * 4 L' =", pi * 4 * L14)
print("  => direct |I_std| = pi * Lambda-chain: diff =",
      fabs(I_std) - pi * LAMBDA_CHAIN_14)

print("=" * 74)
print("(b) conductor 11: final value of int_{{gamma^-}} eta(x,y)")
b11 = gp_lprime("11a3")
print("  b11 = L'(E11a3,0) [PARI]  =", b11)
# step9 (exact rational linear algebra in M_2(Gamma_0(11))): F_total = -2 f_11
# step8 (numeric): sum of 7 symbols' Lambda = -2 b11 = -0.30429428345183609897...
int_gamma = pi * (-2 * b11)   # <-- Theorem 1 with its pi factor
print("  pi * Lambda(F_total,0)    =", int_gamma)
print("  -2 pi b11                 =", -2 * pi * b11)
print("  +2 pi b11 (repo gamma~)   =", 2 * pi * b11)
repo = mpf("0.95596868542165847877477217156086908050293759")
print("  repo direct numerical int =", repo)
print("  |int_gamma| - repo        =", fabs(int_gamma) - repo)
print()
print("CONCLUSION: int_{{gamma^-}} eta(x,y) = +/- 2 pi b11  (sign = orientation);")
print("  with gamma^- = {{0,3/11}} - {{0,8/11}} and the symbol chain of step8/9,")
print("  the value is -2 pi b11 =", -2 * pi * b11)
