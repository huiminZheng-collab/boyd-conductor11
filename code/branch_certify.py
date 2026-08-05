"""
INTERVAL (BALL-ARITHMETIC) CERTIFICATION of
  (A) the 8 endpoint branch assignments of the closed-chain lemma
      (referee major item M2), and
  (B) the modulus ordering of Proposition 3.1,
for S_0 (k=0, conductor 11).  Everything below is rigorous.
Run:  .venv/Scripts/python code/branch_certify.py
Output: notes/attack12-branch.txt

Setting.  y solves  y^2 + B y + x^3 = 0,  B = x^2 + 1,  x = e^{i theta},
discriminant D = B^2 - 4x^3 = x^4 - 4x^3 + 2x^2 + 1;  c = pi/2.

(A) At theta = c (x = i): B = 0 and y^2 = -x^3 = i, roots +/- e^{i pi/4};
    at theta = -c (x = -i): y^2 = -i, roots +/- e^{-i pi/4}.  The four
    exact endpoint points
      P = (i, e^{i pi/4}),  -P = (i, -e^{i pi/4}),
      P-bar = (-i, e^{-i pi/4}),  -P-bar = (-i, -e^{-i pi/4})
    are verified to lie on S_0 by ball arithmetic (residual ball contains
    0; the identities B(+-i) = 0 and (1+-i)^2/2 = +-i are exact).

    Each claimed assignment  y_branch(theta0^{+-}) = point  is certified
    (for eps = 1e-6 and again for eps = 1e-9) by tracking the selected
    sheet over an adaptive subdivision of the WHOLE closing interval
    (theta0, theta0 +/- 2 eps]  (referee item M2 fix: the former
    two-point pinning + continuity argument did not exclude the sheet
    drifting to the other root -- distance 2 -- and back inside the
    interval; a single whole-ball evaluation cannot close this either,
    because |y| = 1 at theta0 makes the modulus test inconclusive there):
      1. on every slab NOT touching theta0: D != 0, a certified sqrt
         variant, STRICT modulus separation |y_big| > 1 > |y_small|
         (part B's machinery), AND the selected sheet's root ball within
         certified distance < 1 of the claim (half the root separation
         at theta0);
      2. on the unique slab touching theta0 (where the modulus test is
         inconclusive): the two root balls are certified disjoint, one
         within < 1 of the claim, the other > 1 from it;
      3. D != 0 on every slab, hence on the whole closed interval, so
         the two sheets are continuous and globally defined there.  The
         sheet Y selected on the outermost slab is the modulus-selected
         root at every shared node, hence on every strict slab, with
         |Y - claim| < 1; on the touching slab Y enters within < 1 of
         the claim and disjointness + continuity keep it in the near
         ball.  So |Y - claim| < 1 on the whole open interval, while
         Y(theta0) is one of the two exact roots (distances 0 and 2
         from the claim) -- the one-sided limit IS the claim.

(B) |y_-(theta)| <= 1 <= |y_+(theta)| for theta in [0, pi], with equality
    |y_-| = 1 only at theta = 0 and theta = pi/2.  Adaptive bisection of
    [0, pi]:
      - on every ball not touching 0 or c: D != 0 certified (no root
        collision on (0, pi]) and STRICT  |y_small| < 1, |y_big| > 1
        certified (each modulus ball strictly on one side of 1);
      - on the special balls touching 0 resp. c (width <= 1e-3): both
        modulus balls are certified to CONTAIN 1, and equality AT the
        point itself is exact algebra -- double root y = -1 at theta = 0
        (p(-1) = p'(-1) = 0 exactly) and y^2 = i at theta = c.  (Recall
        y_+ y_- = x^3 has modulus exactly 1 on |x| = 1, so the moduli are
        reciprocal and equality can occur only where both equal 1.)

Verdict tables are printed for (A) and (B).
"""
from flint import acb, arb, ctx
import math
import sys

ctx.prec = 300

PI = arb.pi()
C2 = PI / 2                        # c = pi/2
I = acb(0, 1)

CLOSE_TOL = arb(1)                 # certified slab distance to the claim
                                   # (< half the root separation 2 at theta0)
SPECIAL_W = 1e-3                   # width cap for the special balls of (B)


def D_and_x(th):
    """th: acb -> (x = e^{i th}, D = x^4 - 4x^3 + 2x^2 + 1)."""
    x = acb.exp(I * th)
    return x, x**4 - 4 * x**3 + 2 * x**2 + 1


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
# exact endpoint points (x, y) on S_0

SQ2 = acb(2).sqrt()
PY = (1 + I) / SQ2                   # e^{i pi/4} = (1+i)/sqrt(2), square = i
PYB = (1 - I) / SQ2                  # e^{-i pi/4}, square = -i
PTS = {"P": (I, PY), "-P": (I, -PY),
       "P-bar": (-I, PYB), "-P-bar": (-I, -PYB)}


def check_exact_points():
    """The 4 endpoint points lie on S_0 (certified residual), and at
    theta = 0 the double root y = -1 is exact (p(-1) = p'(-1) = 0)."""
    for name, (x, y) in PTS.items():
        B = x * x + 1
        r = y * y + B * y + x**3
        assert r.contains(acb(0)) and abs(r).upper() < arb(10) ** (-80), \
            "exact point %s not on S_0" % name
        print("  %-7s on S_0: residual |y^2+By+x^3| <= %s  (contains 0)"
              % (name, abs(r).upper()))
    x, y = acb(1), acb(-1)           # theta = 0: double root y = -1
    B = x * x + 1
    r = y * y + B * y + x**3
    rp = 2 * y + B
    assert r == 0 and rp == 0
    print("  theta=0: p(-1) = 0 and p'(-1) = 0 EXACTLY (double root y = -1)")


def root_balls_point(th):
    """Root balls at an interior point th (acb), big/small by certified
    modulus test.  Returns (x, y_big, y_small, m_big, m_small)."""
    x, D = D_and_x(th)
    assert not D.contains(arb(0)), "D vanishes at %s -- roots collide" % th
    B = x * x + 1
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
    """Certify  y_branch(theta0^{side}) = ptname  at offset scale eps by
    tracking the selected sheet over an adaptive subdivision of the WHOLE
    closing interval (theta0, theta0 +/- 2 eps]  (see the module docstring
    for the argument).  Returns (max certified distance to the claim,
    |D| lower bound over the slabs, number of strict slabs)."""
    xt, yt = PTS[ptname]
    end = th0 + 2 * dirn * eps
    lo, hi = (end, th0) if dirn < 0 else (th0, end)
    worst = arb(0)
    dlow = None
    nstrict = 0
    # NOTE: arb balls never compare equal (radius > 0), so the
    # touches-theta0 property is tracked by an explicit flag: for
    # dirn = +1 theta0 is the LEFT endpoint of the touching slab,
    # for dirn = -1 the RIGHT one.
    mid0 = (lo + hi) / 2            # seed one bisection: the sheet label
    stack = [(mid0, hi, dirn < 0), (lo, mid0, dirn > 0)]
    while stack:                      # must be anchored on a strict slab
        a, b, touch = stack.pop()
        w = b - a
        assert w > arb(10) ** (-30), "slab subdivision did not converge"
        mid = (a + b) / 2

        def refine():
            stack.append((a, mid, touch and dirn > 0))
            stack.append((mid, b, touch and dirn < 0))

        x, Db = D_and_x(acb(ball_interval(a, b)))
        if Db.contains(arb(0)):          # root collision possible: refine
            refine()
            continue
        d = abs(Db).lower()
        dlow = d if dlow is None else min(dlow, d)
        B = x * x + 1
        if touch:
            # special slab touching theta0: |y| = 1 AT theta0, so the
            # modulus test is inconclusive here.  Use a cut-avoiding sqrt
            # variant (k=1 needs this: D(theta0) = -4 lies ON the
            # principal cut; for k=0 D(theta0) = 4i, kind = 'p').
            kind = cut_avoidance(Db)
            if kind is None:             # cut ambiguity: refine
                refine()
                continue
            u = variant_sqrt(Db, kind)
            y1, y2 = (-B + u) / 2, (-B - u) / 2
            d1, d2 = abs(y1 - yt), abs(y2 - yt)
            near, far = (d1, d2) if d1.mid() < d2.mid() else (d2, d1)
            if not (abs(y1 - y2).lower() > 0 and near.upper() < CLOSE_TOL
                    and far.lower() > CLOSE_TOL):
                refine()
                continue
            dx = abs(x - xt).upper()
            assert dx < CLOSE_TOL, "x-ball strays from the claimed corner"
            worst = max(worst, dx, near.upper())
            continue
        kind = cut_avoidance(Db)
        if kind is None:                 # cut ambiguity: refine
            refine()
            continue
        u = variant_sqrt(Db, kind)
        y1, y2 = (-B + u) / 2, (-B - u) / 2
        m1, m2 = abs(y1), abs(y2)
        if (m1 - 1).contains(arb(0)) or (m2 - 1).contains(arb(0)):
            refine()                     # modulus test inconclusive: refine
            continue
        yb, ys = (y1, y2) if m1.lower() > 1 else (y2, y1)
        ysel = yb if want_big else ys
        dx = abs(x - xt).upper()
        dy = abs(ysel - yt).upper()
        if not (dx < CLOSE_TOL and dy < CLOSE_TOL):
            refine()                     # ball too coarse: refine
            continue
        worst = max(worst, dx, dy)
        nstrict += 1
    assert nstrict > 0, "no strict slab -- sheet label never anchored"
    return worst, dlow, nstrict


def part_A():
    print("(A) closed-chain endpoint branch assignments")
    print("    whole-interval sheet tracking on (th0, th0+/-2eps] (referee M2 fix):")
    print("    strict slabs certify |y_big| > 1 > |y_small| AND the selected")
    print("    sheet's ball within distance < 1 of the claim; the slab touching")
    print("    th0 certifies disjoint root balls, one < 1 and one > 1 from it.")
    all_ok = True
    for eps in (arb("1e-6"), arb("1e-9")):
        print("  eps = %s:" % eps.mid())
        print("  %-15s %-7s |dist to claim|  |D| >= on [th0, th0+/-2eps]  strict slabs  ok"
              % ("claim", "point"))
        for name, th0, dirn, want_big, ptname in CASES:
            try:
                worst, dlow, nslab = certify_case(th0, dirn, want_big,
                                                  ptname, eps)
                ok = bool(worst < CLOSE_TOL)
                print("  %-15s %-7s %14.3e  %28.6f  %12d  %s"
                      % (name, ptname, float(worst), float(dlow), nslab, ok))
            except AssertionError as exc:
                ok = False
                print("  %-15s %-7s FAILED: %s" % (name, ptname, exc))
            all_ok = all_ok and ok
    print("  argument: D != 0 on every slab, so both sheets are continuous on")
    print("  the closed interval; the tracked sheet stays within distance < 1")
    print("  of the claim on the WHOLE open interval (strict slabs by modulus")
    print("  separation + ball distance, the touching slab by disjoint balls")
    print("  + continuity), and the two exact roots at th0 are 2 apart")
    print("  ==> the one-sided limit is the claimed point, in all 8 cases.")
    print("  VERDICT (A):", "ALL 8 ASSIGNMENTS CERTIFIED" if all_ok
          else "FAILURE -- investigate")
    return all_ok


# ---------------------------------------------------------------------------
# (B) modulus ordering of Proposition 3.1 on [0, pi]

def part_B():
    print("(B) Proposition 3.1: |y_-(th)| <= 1 <= |y_+(th)| on [0, pi],")
    print("    equality |y_-| = 1 only at th = 0 and th = pi/2")
    c = math.pi / 2
    strict = 0
    special = []
    stack = [(0.0, math.pi)]
    while stack:
        a, b = stack.pop()
        w = b - a
        assert w > 1e-12, "bisection of [0,pi] did not converge"
        mid = (a + b) / 2
        touch0 = (a == 0.0)
        touchc = (a == c or b == c)
        x, Db = D_and_x(acb(ball_interval(a, b)))
        if touch0 or touchc:
            if w > SPECIAL_W:
                stack.append((a, mid))
                stack.append((mid, b))
                continue
            # special ball: certify both modulus balls contain 1
            B = x * x + 1
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
        B = x * x + 1
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
    print("  equality at the special points is EXACT algebra:")
    print("    th = 0:    double root y = -1, |y| = 1  (p(-1)=p'(-1)=0)")
    print("    th = pi/2: B(i) = 0, y^2 = i, |y| = 1")
    ok = strict > 0 and len(special) == 3
    print("  VERDICT (B):", "CERTIFIED" if ok else "FAILURE -- investigate")
    return ok


def main():
    print("=== branch certification, S_0 (k=0) -- referee item M2 ===")
    print("working precision: ctx.prec =", ctx.prec, "bits (no dps override)")
    sys.stdout.flush()
    print("exact endpoint points:")
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
