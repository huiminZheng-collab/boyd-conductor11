"""
BALL-ARITHMETIC SIGN CERTIFICATION of the structural integral, S_1 (k=1,
conductor 17).  Referee item: "the sign of int_gamma~ eta is pinned by a
single floating-point evaluation -- not a proof".  Certified enclosure here.

  S_1: y^2 + B y + x^3 = 0,  x = e^{i th},  B = x^2 + x + 1,
  d = sqrt(B^2 - 4 x^3) (principal sqrt; the two ROOT BALLS are
  {-B+d, -B-d}/2 either way),  y_small(th) = root of smaller modulus
  (continuous; |y_small * y_big| = |x^3| = 1, so neither root vanishes).

  c = 2 pi/3 (fold angle for k=1: at x = exp(+-2 pi i/3), B = 0 and
  S_1 = y^2 + 1 = 0, y = +- i, so |y_small| = |y_big| = 1 there).
  Unlike k=0, D(1) = 5 != 0: no double root at th = 0, log|y| is smooth
  there; the only non-separable point on [0, pi] is the fold th = c.

  J1 = int_0^c log|y_small| dth,  J2 = int_c^pi log|y_small| dth,
  J_split = J1 - J2.   Structural quantities (code/k1_certify.py):
      int_{gamma~} eta = 2 (J1 - J2)          mpmath: -1.8809066251248561...
      ntilde(1) = -(J1 - J2)/pi = -(1/2pi) int_{gamma~} eta
                = +0.2993555868829154... = b_17  (> 0)
  We certify J_split STRICTLY NEGATIVE, which is equivalent to
  int_{gamma~} eta < 0 and to ntilde(1) = +b_17 > 0.

METHOD: identical to code/sign_certify.py (k=0) -- real continuous
integrand, range enclosures h_i * G_i, adaptive bisection, strict
|y_small|/|y_big| separation per piece.  One k=1-specific subtlety: at
the fold th = c one has D = -4 EXACTLY on the principal-sqrt branch cut,
so acb.sqrt returns a radius-~2 ball spanning both sides and the root
balls never certify there.  On non-separable pieces we therefore use an
ALGEBRAIC hull (no sqrt): |y|^2 <= |B||y| + |x|^3 and |y1 y2| = |x|^3
imply  |x|^3/M <= |y| <= M,  M = (|B| + sqrt(|B|^2 + 4|x|^3))/2,  so
log|y_small| lies in [3 log|x| - log M, log M], a ball whose radius is
O(piece width) near the fold.  Budget: radius < 0.05.  Endpoint
c = 2*pi/3 is an arb ball.

ASSERTION: certified ball for J_split has mid + rad < 0  ==>  J_split < 0
(strictly).  Output archived to notes/attack14-sign-k1.txt.
"""
from flint import acb, arb, ctx
import heapq
import math

ctx.prec = 200   # ~60 digits; do NOT set ctx.dps (it would override prec)

PI = arb.pi()
C = 2 * PI / 3
RAD_TARGET = arb("0.005")   # per-integral radius target (budget: < 0.05)

sep_pieces = 0
hull_pieces = 0


def logabs_ball(th_ball):
    """Ball enclosing the range of log|y_small| on the theta-ball th_ball.
    Returns (arb ball, 'sep'|'hull')."""
    x = acb.exp(acb(0, 1) * acb(th_ball))
    B = x * x + x + 1
    d = (B * B - 4 * x**3).sqrt()
    y1 = (-B + d) / 2
    y2 = (-B - d) / 2
    m1, m2 = abs(y1), abs(y2)
    if m1.upper() < m2.lower():
        l = y1.log().real
        if not l.is_nan():
            return l, 'sep'
    elif m2.upper() < m1.lower():
        l = y2.log().real
        if not l.is_nan():
            return l, 'sep'
    # Not separable (or sqrt/log ball unusable).  NOTE: at the fold th = c
    # one has D = -4 exactly ON the principal-sqrt branch cut, so acb.sqrt
    # returns a ball of radius ~2 spanning both sides and the root balls
    # always contain 0 -- the sqrt chain can NEVER certify near c.  Use an
    # ALGEBRAIC hull instead (no sqrt at all): both roots of
    # y^2 + B y + x^3 = 0 satisfy  |y|^2 <= |B||y| + |x|^3, hence
    # |y| <= M := (|B| + sqrt(|B|^2 + 4|x|^3))/2,  and |y1 y2| = |x|^3
    # gives  |y| >= |x|^3 / M.  So log|y_small| in [3 log|x| - log M, log M].
    mB = float(abs(B).upper())
    mxu = float(abs(x).upper())
    mxl = float(abs(x).lower())
    M = (mB + math.sqrt(mB * mB + 4 * mxu**3)) / 2
    hi = math.log(M * 1.0000001 + 1e-30)          # padded upper bound
    lo = 3 * math.log(mxl) - hi
    lo = lo * 1.0000001 - 1e-30                   # padded lower bound (lo<0)
    mid, rad = (lo + hi) / 2, (hi - lo) / 2
    return arb(mid, rad), 'hull'


def piece_contrib(a, b):
    """Certified ball for int_a^b log|y_small| dth (a, b: arb endpoints)."""
    global sep_pieces, hull_pieces
    th_ball = arb((a + b) / 2, (b - a) / 2)   # covers [a,b] incl. endpoint radii
    g, kind = logabs_ball(th_ball)
    if kind == 'sep':
        sep_pieces += 1
    else:
        hull_pieces += 1
    return (b - a) * g


def integrate(lo, hi):
    """Certified ball for int_lo^hi log|y_small| dth by adaptive bisection."""
    lo, hi = arb(lo), arb(hi)
    n0 = 64
    heap = []
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
    print("=== ball-arithmetic sign certification, S_1 (k=1, conductor 17) ===")
    print("working precision: ctx.prec =", ctx.prec, "bits")
    print("per-integral radius target:", RAD_TARGET.str(5))
    J1, n1 = integrate(arb(0), C)
    print("\nJ1 = int_0^{2pi/3} log|y_small| dth")
    print("  int ball: mid =", J1.mid().str(30), " rad <=", J1.rad().str(5))
    print("  pieces evaluated (cumulative):", n1)
    J2, n2 = integrate(C, PI)
    print("\nJ2 = int_{2pi/3}^{pi} log|y_small| dth")
    print("  int ball: mid =", J2.mid().str(30), " rad <=", J2.rad().str(5))
    print("  pieces evaluated (cumulative):", n2)
    J_split = J1 - J2
    int_eta = 2 * J_split
    ntilde = -J_split / PI
    print("\nJ_split = J1 - J2")
    print("  mid =", J_split.mid().str(30))
    print("  rad <=", J_split.rad().str(5))
    print("  mid + rad =", (J_split.mid() + J_split.rad()).str(25))
    print("  mpmath reference: J_split = -0.9404533125624281...")
    print("\nderived structural quantities:")
    print("  int_{gamma~} eta = 2*J_split: mid =", int_eta.mid().str(30),
          " rad <=", int_eta.rad().str(5))
    print("    mpmath reference: -1.8809066251248561... = -2 pi b_17")
    print("  ntilde(1) = -J_split/pi: mid =", ntilde.mid().str(30),
          " rad <=", ntilde.rad().str(5))
    print("    mpmath reference: +0.2993555868829154... = b_17")
    print("  pieces: separated =", sep_pieces, " convex-hull =", hull_pieces)
    assert J_split.rad() < arb("0.05"), "radius budget exceeded"
    if J_split.upper() < 0:
        print("\nPASS: certified ball for J_split lies STRICTLY in (-inf, 0).")
        print("      Hence int_{gamma~} eta < 0 and ntilde(1) = +b_17 > 0 are")
        print("      theorems (ball arithmetic, not a float).")
    else:
        print("\nFAIL: ball straddles 0 -- refine further")


if __name__ == "__main__":
    main()
