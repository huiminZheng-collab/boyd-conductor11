"""Recon: behavior of D = B^2 - 4x^3 and root branches on |x|=1, S_0 (k=0).

Determines for each integration interval whether the principal sqrt of D
is continuous (D avoids (-inf,0]) and which sqrt sign gives the big/small
root of  y^2 + B y + x^3 = 0,  B = x^2+1.
"""
from flint import acb, arb, ctx, good
import cmath

ctx.prec = 128
ctx.dps = 40

def vals(th):
    x = acb.exp(acb(0, 1) * arb(str(th)))
    B = x**2 + 1
    D = B**2 - 4 * x**3
    s = D.sqrt()          # principal branch
    yp = (-B + s) / 2
    ym = (-B - s) / 2
    return x, B, D, yp, ym

N = 2000
import math
print("theta        arg(D)      |yp|        |ym|       big=+/-")
prev_arg = None
jumps = 0
for i in range(1, N + 1):
    th = i * math.pi / N
    x, B, D, yp, ym = vals(th)
    a = cmath.phase(complex(D.mid()))
    # unwrapped continuity check on the sampled arg
    if prev_arg is not None:
        da = a - prev_arg
        if da > math.pi: a -= 2 * math.pi
        elif da < -math.pi: a += 2 * math.pi
        if abs(a - prev_arg) > 0.2:
            jumps += 1
    prev_arg = a
    big_is_plus = abs(complex(yp.mid())) >= abs(complex(ym.mid()))
    if i % 200 == 0 or i in (1, N):
        print(f"{th:.6f}  {a:+.6f}  {abs(complex(yp.mid())):.8f}  {abs(complex(ym.mid())):.8f}  {'+' if big_is_plus else '-'}")
print("sampled-arg discontinuities >0.2 rad:", jumps)

# same scan on negative side
prev_arg = None; jumps = 0
for i in range(1, N + 1):
    th = -i * math.pi / N
    x, B, D, yp, ym = vals(th)
    a = cmath.phase(complex(D.mid()))
    if prev_arg is not None:
        da = a - prev_arg
        if da > math.pi: a -= 2 * math.pi
        elif da < -math.pi: a += 2 * math.pi
        if abs(a - prev_arg) > 0.2: jumps += 1
    prev_arg = a
print("negative side discontinuities >0.2 rad:", jumps)

# min |D| margin from negative real axis on each interval (sampled)
for name, lo, hi in [("(0,c]", 1e-9, math.pi/2), ("[c,pi]", math.pi/2, math.pi),
                     ("[-c,0)", -math.pi/2, -1e-9), ("[-pi,-c]", -math.pi, -math.pi/2)]:
    worst = 1e9; wth = None
    for i in range(2001):
        th = lo + (hi - lo) * i / 2000
        x, B, D, yp, ym = vals(th)
        z = complex(D.mid())
        # distance to negative real axis (treat near-zero separately)
        d = abs(z.imag) if z.real < 0 else abs(z)
        if d < worst: worst, wth = d, th
    print(f"interval {name}: sampled min dist of D to (-inf,0] = {worst:.6g} at th={wth:.6f}")
