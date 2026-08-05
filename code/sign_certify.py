"""
BALL-ARITHMETIC SIGN CERTIFICATION of I_split = b_11, S_0 (k=0).
Referee item: "the sign of I_split is pinned by a single floating-point
evaluation -- not a proof".  This script replaces that evaluation by a
certified enclosure.

  S_0: y^2 + (x^2+1) y + x^3 = 0,  x = e^{i th},  B = x^2+1,
  d = sqrt(B^2 - 4 x^3) (principal sqrt; the branch does not matter --
  the two ROOT BALLS are {-B+d, -B-d}/2 either way),
  y_-(th) = root of smaller modulus (continuous; |y_- y_+| = |x^3| = 1,
  so neither root ever vanishes and log|y| is well defined on both).

  I1 = (1/pi) int_0^{pi/2} log|y_-| dth,   I2 = (1/pi) int_{pi/2}^{pi} ...,
  I_split = I1 - I2.   mpmath 300 dps value: +0.1521471417259180... = b_11.

METHOD (rigorous, crude enclosures suffice -- budget: radius < 0.05):
  * The integrand g(th) = log|y_-(th)| is REAL and continuous on [0, pi]
    (at th=0 the double root y = -1 gives g = 0 with g ~ -sqrt(th/2);
    at the fold th = pi/2, |y_-| = |y_+| = 1 and g ~ -|th - pi/2|).
    Hence for any interval [a,b] and any ball G enclosing the RANGE of g
    on [a,b],  int_a^b g  lies in  (b-a) * G.  Adaptive bisection of
    [0, pi/2] and [pi/2, pi]; each piece contributes the ball h_i * G_i.
  * Range ball per piece: evaluate the whole chain th-ball -> x = e^{i th}
    -> B, d -> root balls y1, y2 with acb ball arithmetic.  If the modulus
    balls |y1|, |y2| are STRICTLY SEPARATED (upper of one < lower of the
    other, certified), the smaller root ball is a certified enclosure of
    y_- on the whole piece, and acb.log(y).real is a true enclosure of
    log|y_-|.  Near th = 0 (double root) and th = pi/2 (fold) separation
    fails; there we take the CONVEX HULL of the two roots' log balls,
    computed entirely in ball arithmetic as l1 + [0,1]*(l2-l1) with
    programmatic containment assertions (no float conversion, no empirical
    padding), so the hull is still a valid -- crude -- range enclosure
    (y_- is one of the two roots) and bisection until the hull-width
    contribution is negligible makes the refinement converge.
  * Endpoints 0, pi/2, pi are arb balls (arb.pi()); piece balls are
    arb((a+b)/2, (b-a)/2) -- NOTE python-flint: arb(a,b) = (mid, rad),
    and ctx.dps would override ctx.prec (neither pitfall present here).

ASSERTION: certified ball for I_split has mid - rad > 0  ==>  I_split > 0
(strictly).  Output archived to notes/attack14-sign-k0.txt.
"""
from flint import acb, arb, ctx
import heapq

ctx.prec = 200   # ~60 digits; do NOT set ctx.dps (it would override prec)

PI = arb.pi()
RAD_TARGET = arb("0.005")   # per-integral radius target (budget: < 0.05)

sep_pieces = 0
hull_pieces = 0
wide_pieces = 0


def logabs_ball(th_ball):
    """Ball enclosing the range of log|y_-| on the theta-ball th_ball.
    Returns (arb ball, 'sep'|'hull') or (None, 'wide') if a root ball
    contains 0 (log undefined there) and the piece must be refined."""
    x = acb.exp(acb(0, 1) * acb(th_ball))
    B = x * x + 1
    d = (B * B - 4 * x**3).sqrt()
    y1 = (-B + d) / 2
    y2 = (-B - d) / 2
    l1, l2 = y1.log().real, y2.log().real
    if l1.is_nan() or l2.is_nan():
        # a root ball contains 0 (piece too coarse): log undefined.
        # Report 'wide' so the caller refines; |y_-| >= sqrt(2)-1 > 0 on
        # [0,pi] (|y_- y_+| = 1, |y_+| <= sqrt(2)+1), so refinement
        # necessarily reaches pieces whose root balls exclude 0.
        return None, 'wide'
    m1, m2 = abs(y1), abs(y2)
    if m1.upper() < m2.lower():
        return l1, 'sep'
    if m2.upper() < m1.lower():
        return l2, 'sep'
    # not separable: convex hull of the two roots' log|.| balls, built
    # ENTIRELY within ball arithmetic: with t ranging over the ball [0,1],
    # H = l1 + t*(l2-l1) contains l1 (t=0), l2 (t=1) and every convex
    # combination, since ball subtraction contains all pairwise
    # differences.  No float conversion and no empirical padding is
    # involved; containment is asserted programmatically below.
    t = arb("0.5", "0.5")
    H = l1 + t * (l2 - l1)
    assert H.contains(l1) and H.contains(l2), "hull failed to cover inputs"
    return H, 'hull'


def piece_contrib(a, b):
    """Certified ball for int_a^b log|y_-| dth (a, b: arb endpoints).
    On too-coarse pieces (root ball contains 0) returns a huge-radius
    placeholder that forces further refinement; it is always re-bisected
    before the final summation (its radius is the largest in the heap)."""
    global sep_pieces, hull_pieces, wide_pieces
    th_ball = arb((a + b) / 2, (b - a) / 2)   # covers [a,b] incl. endpoint radii
    g, kind = logabs_ball(th_ball)
    if kind == 'wide':
        wide_pieces += 1
        return arb(0, 1e9)
    if kind == 'sep':
        sep_pieces += 1
    else:
        hull_pieces += 1
    return (b - a) * g


def integrate(lo, hi):
    """Certified ball for int_lo^hi log|y_-| dth by adaptive bisection."""
    lo, hi = arb(lo), arb(hi)
    n0 = 64
    heap = []          # (-radius, counter, a, b, contrib)
    sum_rad = 0.0      # tracked as a float: ball subtraction would ADD radii
    cnt = 0
    for k in range(n0):
        a = lo + (hi - lo) * k / n0
        b = lo + (hi - lo) * (k + 1) / n0
        c = piece_contrib(a, b)
        sum_rad += float(c.rad())
        heap.append((-float(c.rad()), cnt, a, b, c))
        cnt += 1
    heapq.heapify(heap)
    while sum_rad > float(RAD_TARGET):
        _, _, a, b, c = heapq.heappop(heap)
        sum_rad -= float(c.rad())
        m = (a + b) / 2
        for aa, bb in ((a, m), (m, b)):
            cc = piece_contrib(aa, bb)
            sum_rad += float(cc.rad())
            heapq.heappush(heap, (-float(cc.rad()), cnt, aa, bb, cc))
            cnt += 1
        assert cnt < 1 << 22, "refinement did not converge"
    assert sum_rad < 1e8, "a too-coarse placeholder piece survived refinement"
    total = arb(0)     # single final summation over all active pieces
    for _, _, _, _, c in heap:
        total += c
    return total, cnt


def main():
    print("=== ball-arithmetic sign certification of I_split = b_11 (k=0) ===")
    print("working precision: ctx.prec =", ctx.prec, "bits")
    print("per-integral radius target:", RAD_TARGET.str(5))
    I1, n1 = integrate(arb(0), PI / 2)
    print("\nI1_raw = int_0^{pi/2} log|y_-| dth  (NOT divided by pi)")
    print("  int ball: mid =", I1.mid().str(30), " rad <=", I1.rad().str(5))
    print("  pieces evaluated (cumulative):", n1)
    I2, n2 = integrate(PI / 2, PI)
    print("\nI2_raw = int_{pi/2}^{pi} log|y_-| dth  (NOT divided by pi)")
    print("  int ball: mid =", I2.mid().str(30), " rad <=", I2.rad().str(5))
    print("  pieces evaluated (cumulative):", n2)
    I_split = (I1 - I2) / PI
    print("\nI_split = (I1 - I2)/pi")
    print("  mid =", I_split.mid().str(30))
    print("  rad <=", I_split.rad().str(5))
    print("  mid - rad =", (I_split.mid() - I_split.rad()).str(25))
    print("  mpmath 300-dps reference: +0.1521471417259180... (= b_11)")
    print("  pieces: separated =", sep_pieces, " convex-hull =", hull_pieces,
          " too-coarse (auto-refined) =", wide_pieces)
    assert I_split.rad() < arb("0.05"), "radius budget exceeded"
    if I_split.lower() > 0:
        print("\nPASS: certified ball for I_split lies STRICTLY in (0, +inf).")
        print("      I_split > 0 is now a theorem (ball arithmetic, not a float).")
    else:
        print("\nFAIL: ball straddles 0 -- refine further")


if __name__ == "__main__":
    main()
