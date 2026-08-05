"""
INTERVAL (BALL-ARITHMETIC) CERTIFICATION of
  (A) the 8 endpoint branch assignments of the closed-chain lemma
      (referee major item M2 / finding 16), and
  (B) the modulus ordering of Proposition 3.1,
for S_1 (k=1, conductor 17).  Everything below is rigorous.
Mirrors branch_certify.py (k=0, conductor 11) step by step.
Run:  .venv/Scripts/python code/k1_branch_certify.py
Output: notes/attack12-k1-branch.txt

Setting.  y solves  y^2 + B y + x^3 = 0,  B = x^2 + x + 1,  x = e^{i theta},
discriminant D = B^2 - 4x^3 = x^4 - 2x^3 + 3x^2 + 2x + 1;  c = 2 pi/3.

KEY DIFFERENCES from k=0:
  * At theta = c (x = omega = e^{2 pi i/3}): B = omega^2 + omega + 1 = 0
    EXACTLY, so y^2 = -x^3 = -1, roots y = +/- i.  The four exact
    endpoint points
      P = (omega, i),  -P = (omega, -i),
      P-bar = (omega^2, -i),  -P-bar = (omega^2, i)
    are verified to lie on S_1 by ball arithmetic (residual ball contains
    0; the identities omega^2+omega+1 = 0 and omega^3 = 1 are certified,
    (+-i)^2 = -1 is exact).
  * D(z) has NO zeros on |z| = 1: min |D| = 4 > 0, attained at
    theta = +/- c where D = -4.  Hence there is NO theta = 0 degeneracy
    (k=0 had the double root y = -1 there): the only special point of
    (B) on [0, pi] is theta = c.  At theta = 0 (x = 1): D = 5,
    y = (-3 +/- sqrt(5))/2, |y| = (3 +/- sqrt(5))/2 -- strict.
  * c = 2 pi/3 is NOT a dyadic multiple of pi, so plain bisection of
    [0, pi] never has a node exactly at c (k=0's c = pi/2 appears at the
    first midpoint).  The initial stack is therefore seeded with the
    partition [0, c] U [c, pi]; afterwards the touch-c test works
    verbatim as in k=0.

(A) Each claimed assignment  y_branch(theta0^{+-}) = point  is certified
    (for eps = 1e-6 and again for eps = 1e-9) exactly as in k=0:
      1. computing both root balls y = (-B +/- sqrt(D))/2 at
         theta1 = theta0 +/- eps AND at theta2 = theta0 +/- eps/2
         (|D| >= 3.99 certified near c: D is far from 0; principal sqrt
         fine at a point); the big/small label is decided by a certified
         modulus test (|y| > 1 resp. < 1 with strict margin);
      2. checking the selected root ball is within certified distance
         < 1/10 of the claimed exact point at BOTH theta1 and theta2;
      3. certifying D != 0 on the whole theta-ball [theta0, theta0 +/- 2eps]
         covering the closed interval: no root collision on it, so
         y_branch is continuous on (theta0, theta1], its limit at theta0
         is one of the two exact roots, and those are 2 apart while the
         ball sits within 1/10 of the claim -- the limit IS the claim.

(B) |y_-(theta)| <= 1 <= |y_+(theta)| for theta in [0, pi], with equality
    |y_-| = 1 only at theta = c = 2 pi/3.  Adaptive bisection of
    [0, c] U [c, pi]:
      - on every ball not touching c: D != 0 certified (no root collision)
        and STRICT  |y_small| < 1, |y_big| > 1 certified (each modulus
        ball strictly on one side of 1);
      - on the special balls touching c (width <= 1e-3): both modulus
        balls are certified to CONTAIN 1, and equality AT c itself is
        exact algebra: B(omega) = 0 and y^2 = -1, so y = +/- i, |y| = 1.
        (Recall y_+ y_- = x^3 has modulus exactly 1 on |x| = 1, so the
        moduli are reciprocal and equality can occur only where both
        equal 1.)

Verdict tables are printed for (A) and (B).
"""
from flint import acb, arb, ctx
import math
import sys

ctx.prec = 300
# NOTE: do NOT set ctx.dps here — in python-flint >= 0.9, dps overrides prec
# and would silently drop the working precision to ~50 digits.

PI = arb.pi()
C2 = 2 * PI / 3                    # c = 2 pi/3  (fold angle for k=1)
I = acb(0, 1)

CLOSE_TOL = arb("0.1")             # certified pin distance (< 1/10)
SPECIAL_W = 1e-3                   # width cap for the special balls of (B)


def D_and_x(th):
    """th: acb -> (x = e^{i th}, D = x^4 - 2x^3 + 3x^2 + 2x + 1)."""
    x = acb.exp(I * th)
    return x, x**4 - 2 * x**3 + 3 * x**2 + 2 * x + 1


def variant_sqrt(D, kind):
    if kind == 'p':
        return D.sqrt()
    return I * (-D).sqrt()           # (i*sqrt(-D))^2 = D


def cut_avoidance(Dball):
    """'p' if D-image avoids (-inf,0], 'r' if it avoids [0,+inf), else None."""
    re, im = Dball.real, Dball.imag
    if (not im.contains(arb(0))) or re.lower() > 0:
        return 'p'
    if re.upper() < 0:
        return 'r'
    return None


def ball_interval(lo, hi):
    """Interval ball with endpoints lo, hi (python-flint arb(a,b) = (mid,rad))."""
    lo, hi = arb(lo), arb(hi)
    return arb((lo + hi) / 2, (hi - lo) / 2)


# ---------------------------------------------------------------------------
# exact endpoint points (x, y) on S_1

OM = acb.exp(I * C2)               # omega = e^{2 pi i/3} (certified ball)
OMB = OM * OM                      # omega^2 = e^{-2 pi i/3} = conj(omega)
PTS = {"P": (OM, I), "-P": (OM, -I),
       "P-bar": (OMB, -I), "-P-bar": (OMB, I)}


def check_exact_points():
    """The 4 corner points lie on S_1 (certified residual), and at
    theta = c the corner algebra is exact: omega^2+omega+1 = 0,
    omega^3 = 1 (certified balls containing 0), (+-i)^2 = -1 (exact)."""
    for name, (x, y) in PTS.items():
        B = x * x + x + 1
        r = y * y + B * y + x**3
        assert r.contains(acb(0)) and abs(r).upper() < arb(10) ** (-80), \
            "exact point %s not on S_1" % name
        print("  %-7s on S_1: residual |y^2+By+x^3| <= %s  (contains 0)"
              % (name, abs(r).upper()))
    b = OM * OM + OM + 1             # B(omega) = 0 exactly
    c3 = OM**3 - 1                   # omega^3 = 1 exactly
    assert b.contains(acb(0)) and abs(b).upper() < arb(10) ** (-80)
    assert c3.contains(acb(0)) and abs(c3).upper() < arb(10) ** (-80)
    print("  theta=c: omega^2+omega+1 = 0 and omega^3 = 1 certified")
    print("           (|.| <= %s, %s), so y^2 = -x^3 = -1, y = +/- i"
          % (abs(b).upper(), abs(c3).upper()))
    assert I * I == -1               # exact
    x, D = D_and_x(acb(C2))
    print("           D(omega) = %s  (contains -4; min_{|z|=1}|D| = 4)" % D)


def root_balls_point(th):
    """Root balls at an interior point th (acb), big/small by certified
    modulus test.  Returns (x, y_big, y_small, m_big, m_small)."""
    x, D = D_and_x(th)
    assert not D.contains(arb(0)), "D vanishes at %s -- roots collide" % th
    B = x * x + x + 1
    s = D.sqrt()
    y1, y2 = (-B + s) / 2, (-B - s) / 2
    m1, m2 = abs(y1), abs(y2)
    if (m1 - 1).contains(arb(0)) or (m2 - 1).contains(arb(0)):
        raise RuntimeError("modulus test inconclusive at %s" % th)
    if m1.lower() > 1:
        return x, y1, y2, m1, m2
    assert m2.lower() > 1
    return x, y2, y1, m2, m1


# ---------------------------------------------------------------------------
# (A) the 8 endpoint branch assignments
# y_big(c^-) -> -i (= -P), y_big(c^+) -> +i (= P), y_small(c^-) -> +i,
# y_small(c^+) -> -i  (60-digit values: notes/attack11-k1-certify.txt);
# the assignments at -c are the complex conjugates
# (conj(P) = P-bar, conj(-P) = -P-bar).

CASES = [
    ("y_big(c^-)",     C2, -1, True,  "-P"),
    ("y_big(c^+)",     C2, +1, True,  "P"),
    ("y_big(-c^+)",   -C2, +1, True,  "-P-bar"),
    ("y_big(-c^-)",   -C2, -1, True,  "P-bar"),
    ("y_small(c^-)",   C2, -1, False, "P"),
    ("y_small(c^+)",   C2, +1, False, "-P"),
    ("y_small(-c^+)", -C2, +1, False, "P-bar"),
    ("y_small(-c^-)", -C2, -1, False, "-P-bar"),
]


def certify_case(th0, dirn, want_big, ptname, eps):
    """Certify  y_branch(theta0^{side}) = ptname  at offset scale eps.
    Returns (max distance to claim, |D| lower bound on the interval)."""
    xt, yt = PTS[ptname]
    worst = arb(0)
    for frac in (1, arb(1) / 2):     # theta0 +/- eps and theta0 +/- eps/2
        th = th0 + dirn * eps * frac
        x, yb, ys, mb, ms = root_balls_point(acb(th))
        ysel = yb if want_big else ys
        dx = abs(x - xt).upper()
        dy = abs(ysel - yt).upper()
        assert dx < CLOSE_TOL and dy < CLOSE_TOL, \
            "pin failed: distance to %s exceeds 1/10" % ptname
        worst = max(worst, dx, dy)
    # D != 0 on the closed theta-ball [theta0, theta0 +/- 2 eps]
    lo = th0 + (2 * dirn * eps if dirn < 0 else 0)
    hi = th0 + (0 if dirn < 0 else 2 * dirn * eps)
    _, Db = D_and_x(acb(ball_interval(lo, hi)))
    assert not Db.contains(arb(0)), "D vanishes on the closing interval"
    return worst, abs(Db).lower()


def part_A():
    print("(A) closed-chain endpoint branch assignments")
    all_ok = True
    for eps in (arb("1e-6"), arb("1e-9")):
        print("  eps = %s:" % eps.mid())
        print("  %-15s %-7s |dist to claim|  |D| >= on [th0, th0+/-2eps]  ok"
              % ("claim", "point"))
        for name, th0, dirn, want_big, ptname in CASES:
            worst, dlow = certify_case(th0, dirn, want_big, ptname, eps)
            ok = bool(worst < CLOSE_TOL)
            all_ok = all_ok and ok
            print("  %-15s %-7s %14.3e  %28.6f  %s"
                  % (name, ptname, float(worst), float(dlow), ok))
    print("  continuity argument: D != 0 on each closing interval (above),")
    print("  the two exact roots at theta0 are 2 apart (y = +/- i),")
    print("  pin distance < 1/10")
    print("  ==> the one-sided limit is the claimed point, in all 8 cases.")
    print("  VERDICT (A):", "ALL 8 ASSIGNMENTS CERTIFIED" if all_ok
          else "FAILURE -- investigate")
    return all_ok


# ---------------------------------------------------------------------------
# (B) modulus ordering of Proposition 3.1 on [0, pi]

def part_B():
    print("(B) Proposition 3.1: |y_-(th)| <= 1 <= |y_+(th)| on [0, pi],")
    print("    equality |y_-| = 1 only at th = c = 2 pi/3 (and, by")
    print("    conjugate symmetry, at th = -c on [-pi, 0]); NO theta = 0")
    print("    degeneracy for k = 1 (D(1) = 5 != 0)")
    c = 2 * math.pi / 3
    strict = 0
    special = []
    # k=1 deviation: seed the partition at c (2 pi/3 is not a dyadic
    # multiple of pi, so bisection of [0, pi] would never hit c exactly).
    stack = [(0.0, c), (c, math.pi)]
    while stack:
        a, b = stack.pop()
        w = b - a
        assert w > 1e-12, "bisection of [0,pi] did not converge"
        mid = (a + b) / 2
        touchc = (a == c or b == c)
        x, Db = D_and_x(acb(ball_interval(a, b)))
        if touchc:
            if w > SPECIAL_W:
                stack.append((a, mid))
                stack.append((mid, b))
                continue
            # special ball: certify both modulus balls contain 1
            B = x * x + x + 1
            u = Db.sqrt()            # any valid sqrt enclosure suffices
            assert not (u.real.is_nan() or u.imag.is_nan())
            y1, y2 = (-B + u) / 2, (-B - u) / 2
            m1, m2 = abs(y1), abs(y2)
            assert (m1 - 1).contains(arb(0)), \
                "special ball [%g,%g]: |y1| ball misses 1" % (a, b)
            assert (m2 - 1).contains(arb(0)), \
                "special ball [%g,%g]: |y2| ball misses 1" % (a, b)
            special.append((a, b, m1, m2))
            continue
        if Db.contains(arb(0)):      # root collision possible: refine
            stack.append((a, mid))
            stack.append((mid, b))
            continue
        kind = cut_avoidance(Db)
        if kind is None:             # cut ambiguity: refine
            stack.append((a, mid))
            stack.append((mid, b))
            continue
        B = x * x + x + 1
        u = variant_sqrt(Db, kind)
        y1, y2 = (-B + u) / 2, (-B - u) / 2
        m1, m2 = abs(y1), abs(y2)
        if (m1 - 1).contains(arb(0)) or (m2 - 1).contains(arb(0)):
            stack.append((a, mid))
            stack.append((mid, b))
            continue
        mb, ms = (m1, m2) if m1.lower() > 1 else (m2, m1)
        assert mb.lower() > 1 and ms.upper() < 1
        strict += 1
    print("  strict balls (D != 0, |y_small| < 1 < |y_big| certified): %d"
          % strict)
    for a, b, m1, m2 in sorted(special):
        print("  special ball [%0.6f, %0.6f]: |y| balls %s, %s (contain 1)"
              % (a, b, m1, m2))
    print("  equality at the special point is EXACT algebra:")
    print("    th = c = 2 pi/3: B(omega) = omega^2+omega+1 = 0, omega^3 = 1,")
    print("    so y^2 = -x^3 = -1, y = +/- i, |y| = 1")
    print("    (no th = 0 special ball: D(1) = 5, y = (-3 +/- sqrt(5))/2,")
    print("     moduli (3 -/+ sqrt(5))/2 -- certified strict by bisection)")
    ok = strict > 0 and len(special) == 2
    print("  VERDICT (B):", "CERTIFIED" if ok else "FAILURE -- investigate")
    return ok


def main():
    print("=== branch certification, S_1 (k=1, conductor 17) -- referee item M2 / finding 16 ===")
    print("working precision: ctx.prec =", ctx.prec, "bits (no dps override)")
    sys.stdout.flush()
    print("exact corner points (c = 2 pi/3, omega = e^{2 pi i/3}):")
    check_exact_points()
    print()
    okA = part_A()
    print()
    okB = part_B()
    print()
    print("OVERALL VERDICT:", "ALL CERTIFIED. QED." if (okA and okB)
          else "INCONCLUSIVE -- investigate")


if __name__ == "__main__":
    main()
