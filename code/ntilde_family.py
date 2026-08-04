"""
Second wave: corrected Mahler measure ntilde(k) for the family
    S_k = y^2 + (x^2 + k x + 1) y + x^3
along the closed signed cycle (Samart Lemma 9 construction, cf. closedness_check.py).

Crossings = angles where log|y_big| TOUCHES 0 (g >= 0, tangential zeros) ->
detected by golden-section minimization, integrals split there (cusps).

  m(k)      = (1/2pi) int_{-pi}^{pi} log|y_big|
  ntilde(k) = (1/2pi)[ int_{-c}^{c} - int_{|theta|>c} ] log|y_big|   (single pair +-c)
For k = 0 recover ntilde = -b_11 and m = m(S_0) = 0.40560295591501...
"""
from mpmath import (mp, mpf, mpc, pi, sqrt, exp, fabs, log, nstr,
                    quad, pslq)
import sys
sys.path.insert(0, ".")
from b11 import b11

mp.dps = 50

def make_ybig(k):
    def ybig(th):
        x = exp(mpc(0, 1)*th)
        B = x**2 + k*x + 1
        d = sqrt(B**2 - 4*x**3)
        a, bb = (-B + d)/2, (-B - d)/2
        return a if fabs(a) >= fabs(bb) else bb
    return ybig

def golden_min(g, a, b, tol=mpf('1e-35')):
    gr = (sqrt(5) - 1)/2
    c, d = b - gr*(b-a), a + gr*(b-a)
    fc, fd = g(c), g(d)
    while fabs(b-a) > tol:
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - gr*(b-a); fc = g(c)
        else:
            a, c, fc = c, d, fd
            d = a + gr*(b-a); fd = g(d)
    t = (a+b)/2
    return t, g(t)

def find_crossings(g, n=4000):
    """tangential zeros of g on [0, pi]: local minima with value ~ 0."""
    xs = [pi * mpf(j) / n for j in range(0, n+1)]
    vs = [g(t) for t in xs]
    out = []
    if vs[0] < mpf('1e-6'):
        out.append(mpf(0))
    for i in range(1, n):
        if vs[i] <= vs[i-1] and vs[i] <= vs[i+1] and vs[i] < mpf('0.05'):
            t, v = golden_min(g, xs[i-1], xs[i+1])
            if v < mpf('1e-10'):
                out.append(t)
    if vs[-1] < mpf('1e-6'):
        out.append(pi)
    ded = []
    for r in sorted(out):
        if not ded or fabs(r - ded[-1]) > mpf('1e-20'):
            ded.append(r)
    return ded

def meas(k):
    ybig = make_ybig(k)
    g = lambda th: log(fabs(ybig(th)))
    cs = find_crossings(g)          # in [0, pi]
    cs_open = [c for c in cs if c > mpf('1e-30') and pi - c > mpf('1e-30')]
    # split points for integration on [-pi, pi]: ALL tangential zeros incl. 0
    zs = set()
    for c in cs_open:
        zs.add(c); zs.add(-c)
    if cs and cs[0] == 0:
        zs.add(mpf(0))
    pts = [-pi] + sorted(zs) + [pi]
    m = sum(quad(g, [pts[i], pts[i+1]]) for i in range(len(pts)-1)) / (2*pi)
    if len(cs_open) == 0:
        nt, tag = m, "no crossings"
    elif len(cs_open) == 1:
        c = cs_open[0]
        inner = (quad(g, [-c, 0]) + quad(g, [0, c])) if mpf(0) in zs else quad(g, [-c, c])
        nt = (inner - quad(g, [-pi, -c]) - quad(g, [c, pi])) / (2*pi)
        tag = f"c/pi = {nstr(c/pi, 15)}"
    else:
        nt, tag = None, f"MULTIPLE: {[nstr(c/pi, 10) for c in cs_open]}"
    return m, nt, tag, [nstr(c, 12) for c in cs]

print(f"{'k':>3} | {'m(k)':>27} | {'ntilde(k)':>27} | crossings(0,pi)")
b = b11(50)
results = {}
for k in [-3, -2, -1, 0, 1, 2, 3]:
    try:
        m, nt, tag, cs = meas(k)
        results[k] = (m, nt)
        print(f"{k:>3} | {nstr(m, 24):>27} | "
              f"{nstr(nt, 24) if nt is not None else '---':>27} | {tag}")
    except Exception as e:
        print(f"{k:>3} | FAILED: {e}")

print("\nchecks:  m(0) vs m(S0)=0.40560295591501040...,  ntilde(0) vs -b_11")
m0, nt0 = results[0]
print("  m(0)  =", nstr(m0, 30))
print("  nt(0) =", nstr(nt0, 30), "  -b_11 =", nstr(-b, 30))

print("\nPSLQ: ntilde(k) vs [b_11, pi, log 2, log 3, 1]  (tol 1e-18)")
for k, (m, nt) in results.items():
    if nt is None:
        continue
    r = pslq([nt, b, pi, log(2), log(3), mpf(1)], tol=mpf('1e-18'), maxcoeff=10**8)
    print(f"  k={k:>2}: {r}")
