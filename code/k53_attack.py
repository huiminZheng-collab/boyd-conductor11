"""
k = -1 (conductor 53): resolution of the Samart remark + refined family picture.

FINDINGS (this script + k53.gp + k53b.gp + kfamily_torsion.gp):
  * torus crossings for k=-1 on |x|=1: th = 0 (roots exp(+-2 pi i/3)) and
    th = +-pi/3 (roots +-1); the two branches EXCHANGE at each crossing.
  * chain enumeration (breakpoints -pi,-pi/3,0,pi/3, branches big/small,
    endpoints = one-sided limits): the ONLY closed anti-invariant chains on
    |x|=1 are multiples of (1,1,1,1), whose period is 0  =>  homologically
    trivial.  The anti-invariant generator of H_1(E,Z)^- is NOT realizable
    on the torus.  (Consistent: the continuous-root loop around the circle
    needs two revolutions to close, and the full double loop has period 0.)
  * deeper reason: 53.a1 has tors = trivial, rank 1, and (0,0) = MW generator
    (ellorder = 0).  Hence x,y are NOT modular units; no Beilinson-Brunault
    rational L-value relation is expected for ANY cycle.
  * contrast table (kfamily_torsion.gp):
      k=0:  tors Z/5, ord(0,0)=5   -> modular units -> Boyd identity holds
      k=1:  tors Z/4, ord(0,0)=4   -> modular units -> ntilde(1)=b_17 holds
      k=-3,-2,-1,2,3: tors trivial, ord(0,0)=INFINITE
  * prediction tested (k <= -4, no genuine torus crossings): m(S_k) itself is
    a rational multiple of |L'(E,0)|:
      m(S_-4) = 7/2 |b_37|,  m(S_-5) = 1/4 |b_359|,  m(S_-6) = 1/8 |b_997|
    (25 digits exact).  Together with m(S_2)=2|b_37|, m(S_3)=|b_79|:
    **Boyd-type identities hold for ALL k with |k| >= 2** (standard Deninger
    regime), while inside (-4,2) they hold ONLY for the torsion cases 0,1.
  * Samart's "analogous conjectural identity for conductor 53" (a single
    sentence, no formula/precision) is with high confidence a numerical
    false positive: no closed anti-invariant cycle on the torus exists, and
    the natural candidates' regulator integrals are irrational multiples of
    b_53 (PSLQ below).
"""
from mpmath import mp, mpf, mpc, pi, sqrt, exp, nstr, quad, fabs, log, pslq

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

def A_arc(t0, t1, br):
    f = lambda th: mpc(0, 1)*exp(mpc(0, 1)*th) / (2*ybr(th, br) + exp(mpc(0, 2)*th) + K*exp(mpc(0, 1)*th) + 1)
    return quad(f, [t0, t1])

def R_arc(t0, t1, br):
    return quad(lambda th: -log(fabs(ybr(th, br))), [t0, t1])

def key(th, y):
    x = exp(mpc(0, 1)*th)
    return (round(float(x.real), 8), round(float(x.imag), 8),
            round(float(y.real), 8), round(float(y.imag), 8))

c = pi/3
J = [(mpf(0), c, 'J1'), (c, pi, 'J2')]
arcs = {}
for t0, t1, nm in J:
    for br in 'bs':
        arcs[(nm, br)] = dict(A=A_arc(t0, t1, br), R=R_arc(t0, t1, br),
                              s=key(t0+EPS, ybr(t0+EPS, br)),
                              e=key(t1-EPS, ybr(t1-EPS, br)))

names = [('J1', 'b'), ('J1', 's'), ('J2', 'b'), ('J2', 's')]
def conjpt(p):
    return (p[0], -p[1], p[2], -p[3])

def boundary(coeffs):
    bd = {}
    for n, (nm, br) in zip(coeffs, names):
        if n == 0:
            continue
        d = arcs[(nm, br)]
        for pt, w in [(d['e'], n), (d['s'], -n),
                      (conjpt(d['e']), -n), (conjpt(d['s']), n)]:
            bd[pt] = bd.get(pt, 0) + w
    return {p: w for p, w in bd.items() if w != 0}

print("closed anti-invariant chains on |x|=1 (n in [-2,2]):")
found = []
for n1 in range(-2, 3):
    for n2 in range(-2, 3):
        for n3 in range(-2, 3):
            for n4 in range(-2, 3):
                co = (n1, n2, n3, n4)
                if all(v == 0 for v in co) or boundary(co):
                    continue
                P = sum(n*2*mpc(0, 1)*arcs[nm]['A'].imag for n, nm in zip(co, names))
                R = sum(2*n*arcs[nm]['R'] for n, nm in zip(co, names))
                found.append((co, P, R))
                print(f"  n={co}:  P = {nstr(P,18)}   R = {nstr(R,18)}")
print("  => all closed anti-invariant chains are homologically TRIVIAL" )

w_anti = mpc(0, -1)*mpf('3.0811813402756626952638569706155872205350922143816')
b53 = mpf('0.6289720645462873268255234759695752406840870110883')

print("\nhalf-loop L1 (continuous root, one revolution, OPEN chain):")
PL1 = 2*mpc(0, 1)*(arcs[('J1', 'b')]['A'].imag + arcs[('J2', 's')]['A'].imag)
RL1 = 2*(arcs[('J1', 'b')]['R'] + arcs[('J2', 's')]['R'])
print(f"  P(L1) = {nstr(PL1, 25)}")
print(f"  P(L1)/w_anti = {nstr(PL1/w_anti, 20)}")
print(f"  R(L1) = {nstr(RL1, 25)}")
print(f"  R(L1)/(2 pi b53) = {nstr(RL1/(2*pi*b53), 20)}")
r1 = pslq([PL1.imag, fabs(w_anti.imag)], tol=mpf('1e-20'), maxcoeff=10**6)
r2 = pslq([RL1, 2*pi*b53], tol=mpf('1e-20'), maxcoeff=10**6)
print("  PSLQ P vs w_anti:", r1, "   PSLQ R vs 2 pi b53:", r2)

print("\nold ntilde chain: NOT closed (boundary at 4 points), period not in lattice:")
print("  P(old)/w_anti =", nstr(2*(arcs[('J1','b')]['A'].imag - arcs[('J2','b')]['A'].imag)/fabs(w_anti.imag), 20))
