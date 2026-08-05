"""
INTERVAL (BALL-ARITHMETIC) CERTIFICATION of the cycle lemma, S_1 (k=1,
conductor 17).  Adapted from n1_interval.py (k=0, conductor 11).
Everything below is rigorous:

  * All integrals are computed with Arb's acb.integral (certified enclosures).
  * SIMPLIFICATION vs k=0: D(z) = (z^2+z+1)^2 - 4z^3 = z^4 - 2z^3 + 3z^2
    + 2z + 1 satisfies D(1) = 5 != 0, and min_{|z|=1} |D(z)| = 4 > 0
    (attained at the corners z = exp(+-2 pi i/3), D = -4, certified below
    by the subdivision succeeding with strict cut avoidance).  Hence u =
    sqrt(D) never vanishes on |x| = 1: there is NO 1/sqrt(theta) endpoint
    singularity at theta=0 and the whole theta = +/- t^2 / Cauchy-tip
    machinery of k=0 is unnecessary.  The inner arc [-c, c], c = 2 pi/3,
    is integrated directly.
  * D(th) touches/crosses the negative real axis ONLY at the corner points
    theta = +- c (D = -4 there, endpoints of the integration arcs), so the
    principal sqrt is not globally analytic: same remedy as k=0 -- adaptive
    subdivision; on each piece use either sqrt(D) or i*sqrt(-D) (analytic
    square root whenever the piece's D-image avoids (-inf,0] resp.
    [0,+inf), CERTIFIED per piece by ball evaluation), and propagate the
    overall sign across nodes with a certified matching test.
  * w_anti via Carlson RF (certified) for the minimal model
    y^2 + x y + y = x^3 - x^2 - x  ([1,-1,1,-1,0], conductor 17, disc = +17
    > 0, THREE real roots of 4x^3 - 3x^2 - 2x + 1 = (x-1)(4x^2+x-1);
    roots isolated by Newton + Rouche).  For e1 > e2 > e3 real:
      w_real = 2 RF(0, e1-e2, e1-e3)   (= PARI w1, checked)
      w_imag = 2i RF(0, e2-e3, e1-e3)  (= -PARI w2 = +2.7457391... i)
    disc > 0 => conj(w2) = -w2 exactly, so the primitive anti-invariant
    period is w_anti = w_imag (k=0 had disc < 0, w_anti = 2w2 - w1).
    Same period lattice as the original cubic model (birational over Q,
    invariant differential dx/u, kappa = 1: ellfromeqn model [1,-1,-1,0,0]
    already has disc = 17 = disc_min; 16-digit numeric lock, see
    notes/attack11-k1-certify.txt; kappa = -1 would only flip the sign of
    the ratio, which stays an integer -- the verdict is unaffected).

Output: P = I_signed + A_out - A_in as a certified ball, w_anti ditto,
the ratio P/w_anti, and the verdict that the ratio ball contains an
integer n0 with |ratio - n0| < 1/2  ==>  period(C')/w_anti = n0 exactly
(a priori integer).  Expected |n0| = 2 (sign = orientation convention).
"""
from flint import acb, arb, acb_poly, ctx
import sys

ctx.prec = 300
ctx.dps = 50

PI = arb.pi()
C2 = 2 * PI / 3                    # c = 2 pi/3  (fold angle for k=1)

TOL_REL = 1e-13
TOL_ABS = 1e-25

call_count = 0


def D_and_x(th):
    """th: acb -> (x = e^{i th}, D = x^4 - 2x^3 + 3x^2 + 2x + 1)."""
    x = acb.exp(acb(0, 1) * th)
    return x, x**4 - 2 * x**3 + 3 * x**2 + 2 * x + 1


def variant_sqrt(D, kind):
    if kind == 'p':
        return D.sqrt()
    return acb(0, 1) * (-D).sqrt()   # (i*sqrt(-D))^2 = D


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


def square_ball(rad):
    r = ball_interval(-rad, rad)
    return acb(r, r)


def branch_sign(th_mid, want_big):
    """Certified branch determination at an interior point th_mid (arb):
    u = sigma*sqrt(D) with sigma in {+1,-1} chosen so y = (-B+u)/2 satisfies
    |y| > 1 (big) resp. < 1 (small).  B = x^2 + x + 1."""
    x, D = D_and_x(acb(th_mid))
    B = x * x + x + 1
    s = D.sqrt()
    out = {}
    for sigma in (1, -1):
        u = s if sigma == 1 else -s
        y = (-B + u) / 2
        out[sigma] = abs(y)
    m1, m_1 = out[1], out[-1]
    if (m1 - 1).contains(arb(0)) or (m_1 - 1).contains(arb(0)):
        raise RuntimeError("branch test inconclusive at %s" % th_mid)
    big_is_plus = m1.lower() > 1
    return 1 if (big_is_plus == want_big) else -1


def gl(f, a, b):
    global call_count
    r = acb.integral(f, a, b, rel_tol=TOL_REL, abs_tol=TOL_ABS,
                     depth_limit=300)
    call_count += 1
    if r.real.is_nan() or r.imag.is_nan():
        raise RuntimeError("acb.integral returned NaN on [%s, %s]" % (a, b))
    return r


def arc_integral(t_lo, t_hi, th_fn, jac_fn, want_big):
    """Certified integral of  i*x(th)/u(th) * jac(t) dt  over [t_lo, t_hi],
    th = th_fn(t).  The arc must contain no branch point of u (there are
    none on |x|=1 for k=1)."""
    t_lo, t_hi = arb(t_lo), arb(t_hi)
    n = 8
    while True:
        nodes = [t_lo + (t_hi - t_lo) * i / n for i in range(n + 1)]
        kinds = []
        ok = True
        for j in range(n):
            tb = acb(ball_interval(nodes[j], nodes[j + 1]))
            _, Db = D_and_x(th_fn(tb))
            k = cut_avoidance(Db)
            if k is None:
                ok = False
                break
            kinds.append(k)
        if ok:
            break
        n *= 2
        if n > 1 << 14:
            raise RuntimeError("subdivision failed on [%s,%s]" % (t_lo, t_hi))
    # branch_sign gives sigma relative to the PRINCIPAL sqrt; the first piece
    # may use the 'r' variant i*sqrt(-D), which equals +sqrt(D) or -sqrt(D)
    # depending on the side of the cut.  Convert: match sig*sqrt(D) against
    # +/- variant_sqrt(D, kinds[0]) at the same midpoint (certified test).
    th_mid = th_fn(acb((nodes[0] + nodes[1]) / 2)).mid()
    sig = branch_sign(th_mid, want_big)
    _, Dm = D_and_x(acb(th_mid))
    u_want = Dm.sqrt()
    if sig == -1:
        u_want = -u_want
    kind0 = kinds[0]
    sep = abs(variant_sqrt(Dm, kind0)).lower()
    best = None
    for cand in (1, -1):
        u_c = variant_sqrt(Dm, kind0)
        if cand == -1:
            u_c = -u_c
        if abs(u_c - u_want).upper() < sep:
            assert best is None, "ambiguous initial branch conversion"
            best = cand
    assert best is not None, "initial branch conversion failed"
    sig = best
    total = acb(0)
    for j in range(n):
        a, b, kind = nodes[j], nodes[j + 1], kinds[j]
        s = sig

        def f(t, _a, s=s, kind=kind):
            th = th_fn(t)
            x, D = D_and_x(th)
            u = variant_sqrt(D, kind)
            if s == -1:
                u = -u
            return acb(0, 1) * x * jac_fn(t) / u

        total += gl(f, a, b)
        if j < n - 1:   # propagate sign across the node
            _, Dn = D_and_x(th_fn(acb(nodes[j + 1])))
            u_arr = variant_sqrt(Dn, kind)
            if s == -1:
                u_arr = -u_arr
            nkind = kinds[j + 1]
            sep = abs(variant_sqrt(Dn, nkind)).lower()
            best = None
            for cand in (1, -1):
                u_c = variant_sqrt(Dn, nkind)
                if cand == -1:
                    u_c = -u_c
                if abs(u_c - u_arr).upper() < sep:
                    assert best is None, "ambiguous sign propagation"
                    best = cand
            assert best is not None, "sign propagation failed"
            sig = best
    return total


# ---------------------------------------------------------------------------
# root certification (generic, Newton + Rouche)

def refine_newton(p, dp, z, iters=8):
    z = acb(z)
    for _ in range(iters):
        z = z - p(z) / dp(z)
    return z.mid()


def certify_root(p, dp, z0):
    """Rouche: B(z0, r) contains a root of p if |p(z0)| < r * inf_B |p'|."""
    resid = abs(p(acb(z0))).upper()
    m0 = abs(dp(acb(z0))).lower()
    r = arb(max(1e-50, float(1000 * resid / m0)))
    B = acb(z0) + square_ball(r)
    m = abs(dp(B)).lower()
    assert resid < m * r, "Rouche certificate failed"
    return B


def certify_poly_roots(coeffs_asc, starts=None):
    p = lambda z: sum(c * z**k for k, c in enumerate(coeffs_asc))
    dp = lambda z: sum(k * c * z**(k - 1) for k, c in enumerate(coeffs_asc) if k)
    if starts is None:
        starts = [r.mid() for r in acb_poly(coeffs_asc).roots()]
    return [certify_root(p, dp, refine_newton(p, dp, s)) for s in starts]


# ---------------------------------------------------------------------------

def certify_w_anti():
    """Minimal model [1,-1,1,-1,0]: 4x^3 - 3x^2 - 2x + 1 = (x-1)(4x^2+x-1),
    three real roots (disc = +17).  Starts near 1, (sqrt(17)-1)/8,
    -(sqrt(17)+1)/8; certified by Newton + Rouche."""
    starts = [1.0, 0.3903882032022076, -0.6403882032022076]
    roots = certify_poly_roots([1, -2, -3, 4], starts)   # 4x^3-3x^2-2x+1
    roots.sort(key=lambda r: r.real.mid())               # ascending
    e3, e2, e1 = roots                                   # e1 > e2 > e3
    for r in roots:
        assert abs(r.imag).upper() < arb(10)**(-30), "root not real?!"
    w_real = 2 * acb.elliptic_rf(acb(0), e1 - e2, e1 - e3)
    w_imag = 2 * acb(0, 1) * acb.elliptic_rf(acb(0), e2 - e3, e1 - e3)
    return w_real, w_imag


def main():
    print("=== interval certification, S_1 (k=1, conductor 17) ===")
    sys.stdout.flush()

    print("[1/4] inner arc [-c,c], big branch (no theta=0 singularity):")
    inner_b = arc_integral(-C2, C2, lambda t: t, lambda t: 1, True)
    print("  inner_b  =", inner_b)
    sys.stdout.flush()

    print("[2/4] outer arcs, big branch (cut at corners handled):")
    ob_pos = arc_integral(C2, PI, lambda t: t, lambda t: 1, True)
    print("  [c,pi]   =", ob_pos)
    ob_neg = arc_integral(-PI, -C2, lambda t: t, lambda t: 1, True)
    print("  [-pi,-c] =", ob_neg)
    I_signed = inner_b - (ob_pos + ob_neg)
    print("  I_signed =", I_signed)
    sys.stdout.flush()

    print("[3/4] small branch (A_in, A_out):")
    A_in = arc_integral(-C2, C2, lambda t: t, lambda t: 1, False)
    A_out = arc_integral(C2, PI, lambda t: t, lambda t: 1, False) \
        + arc_integral(-PI, -C2, lambda t: t, lambda t: 1, False)
    print("  A_in     =", A_in)
    print("  A_out    =", A_out)
    sys.stdout.flush()

    print("[4/4] w_anti via Carlson RF (certified roots):")
    w_real, w_imag = certify_w_anti()
    print("  w_real   =", w_real, "  (expect PARI w1 = 3.0941595071022403...)")
    print("  w_imag   =", w_imag, "  (expect -PARI w2 = +2.7457391180897... i)")
    w_anti = w_imag
    sys.stdout.flush()

    P = I_signed + A_out - A_in
    print()
    print("period(C') =", P)
    ratio = P / w_anti
    print("ratio P/w_anti =", ratio)
    n0 = int(round(float(ratio.real.mid())))
    rad = abs(ratio - n0).upper()
    print("nearest integer n0 =", n0, "  |ratio - n0| <=", rad)
    ok_contains = (ratio - n0).contains(acb(0))
    ok_small = rad < 0.5
    print("contains n0:", bool(ok_contains), "  |ratio-n0| < 1/2:", bool(ok_small))
    if ok_contains and ok_small:
        print("VERDICT: period(C')/w_anti = %d  EXACTLY (a priori integer)."
              % n0)
        print("         |class(C')| = 2 * gamma^-  --  cycle lemma certified. QED.")
    else:
        print("VERDICT: INCONCLUSIVE -- investigate")
    print("GL calls:", call_count)

    print()
    print("cross-check: -2*w_anti (mid) =", (-2 * w_anti).mid())
    print("our P (mid)                =", P.mid())
    print("mpmath 60/80-digit gave period(C')/w_anti = 2.0000000000000002 "
          "(w_anti = PARI w2 = -w_imag; sign = orientation)")


if __name__ == "__main__":
    main()
