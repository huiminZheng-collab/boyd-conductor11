"""
Fourth wave: test the structural prediction for k <= -4 (no genuine torus
crossings): m(S_k) itself should be a rational multiple of |L'(E_k,0)|.

Result (25 digits exact):
  m(S_-4) = 7/2 |b_37|     (same conductor 37 as S_2 -- RV-type relation)
  m(S_-5) = 1/4 |b_359|
  m(S_-6) = 1/8 |b_997|
b_N values from PARI (kneg.gp: ellfromeqn -> minimal model -> lfun).
"""
from mpmath import mp, mpf, mpc, pi, sqrt, exp, fabs, log, nstr, quad

mp.dps = 50

def make_ybig(k):
    def ybig(th):
        x = exp(mpc(0, 1)*th)
        B = x**2 + k*x + 1
        d = sqrt(B**2 - 4*x**3)
        a, bb = (-B + d)/2, (-B - d)/2
        return a if fabs(a) >= fabs(bb) else bb
    return ybig

b = {-4: (37, mpf('0.3576204661274976499989425692137969438649688554459')),
     -5: (359, mpf('6.1757484565398443745639817442650057108504649126200')),
     -6: (997, mpf('14.0093879721933903080171234124227422121074090767801'))}

for k, (N, bN) in b.items():
    ybig = make_ybig(k)
    g = lambda th: log(fabs(ybig(th)))
    m = (quad(g, [-pi, 0]) + quad(g, [0, pi]))/(2*pi)
    print(f"k={k}:  N={N}   m = {nstr(m, 25)}")
    print(f"        m/|b_N| = {nstr(m/bN, 25)}")
