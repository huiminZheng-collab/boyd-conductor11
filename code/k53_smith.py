"""
Exact integer-rank computation for the k=-1 (conductor 53) chain space,
replacing the coefficient-box search of k53_attack.py.

The arcs and their endpoint identifications are taken from the same
one-sided limits as k53_attack.py (numeric keys are used ONLY as labels;
distinct endpoints differ by O(1), so 8-decimal rounding cannot collide).
The boundary matrix of the anti-invariant symmetrization is then an exact
+1/0/-1 integer matrix, and its kernel is computed exactly (sympy
nullspace = Smith-type rational arithmetic, no coefficient bound).
"""
from mpmath import mp, mpf, mpc, pi, sqrt, exp
import sympy as sp

mp.dps = 50
K = -1
EPS = mpf('1e-9')

def roots(th):
    x = exp(mpc(0, 1)*th)
    B = x**2 + K*x + 1
    d = sqrt(B**2 - 4*x**3)
    a, b = (-B + d)/2, (-B - d)/2
    return (a, b) if abs(a) >= abs(b) else (b, a)

def ybr(th, br):
    return roots(th)[0 if br == 'b' else 1]

def key(th, y):
    x = exp(mpc(0, 1)*th)
    return (round(float(x.real), 8), round(float(x.imag), 8),
            round(float(y.real), 8), round(float(y.imag), 8))

def conjpt(p):
    return (p[0], -p[1], p[2], -p[3])

c = pi/3
J = [(mpf(0), c, 'J1'), (c, pi, 'J2')]
arcs = {}
for t0, t1, nm in J:
    for br in 'bs':
        arcs[(nm, br)] = dict(s=key(t0+EPS, ybr(t0+EPS, br)),
                              e=key(t1-EPS, ybr(t1-EPS, br)))

names = [('J1', 'b'), ('J1', 's'), ('J2', 'b'), ('J2', 's')]

# distinct endpoint labels
pts = []
def pid(p):
    if p not in pts:
        pts.append(p)
    return pts.index(p)

# boundary matrix rows = points, cols = arcs (anti-invariant symmetrized)
M = {}
for j, (nm, br) in enumerate(names):
    d = arcs[(nm, br)]
    for pt, w in [(d['e'], 1), (d['s'], -1),
                  (conjpt(d['e']), -1), (conjpt(d['s']), 1)]:
        M[(pid(pt), j)] = M.get((pid(pt), j), 0) + w

rows = len(pts)
mat = sp.Matrix(rows, 4, lambda i, j: M.get((i, j), 0))
print("endpoints (exact labels):")
for i, p in enumerate(pts):
    print(f"  P{i} = {p}")
print("boundary matrix (rows=points, cols=[J1b,J1s,J2b,J2s]):")
print(mat)
ns = mat.nullspace()
print("nullspace dimension over Q:", len(ns))
for v in ns:
    print("  generator:", list(v))
print("rank check: rank =", mat.rank(), " (4 - dim ker =", 4 - len(ns), ")")
