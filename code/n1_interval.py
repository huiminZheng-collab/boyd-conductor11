"""
INTERVAL (BALL-ARITHMETIC) CERTIFICATION of the cycle lemma numerics, S_0 (k=0).
Sixth wave, (b).  Everything below is rigorous:

  * All integrals are computed with Arb's acb.integral (certified enclosures).
  * The endpoint singularity u ~ sqrt(theta) at theta=0 is removed by the
    substitution theta = +/- t^2.  The substituted integrand f(t) = 2itx/u is
    ANALYTIC on a disc |t| <= rho (u(t^2) = t*v(t^2), v analytic and nonzero;
    zeros of D(t^2) other than t=0 are at t^2 = theta_j, |theta_j| certified
    > rho^2 from the roots of D(z) = z^4 - 4z^3 + 2z^2 + 1, z = e^{i theta}).
    The tip [0, delta] is summed analytically: f(t) = a0 + t^2 h(t),
      a0 = (i/sigma) exp(sign_theta * i pi/4),
      int_0^delta f = a0*delta + R,  |R| <= H delta^3 / 3,
      |h| <= H = M / (rho^2 (1 - (delta/rho)^2))   (Cauchy),
    M = max_{|t|=rho} |f| certified by covering the circle with balls and
    using |f| = 2|t|/sqrt(|D|) -- branch-independent, no cut issues.
  * D(th) crosses the negative real axis once on (c,pi) and once on (-pi,-c),
    c = pi/2, so principal sqrt is NOT analytic there.  We subdivide
    adaptively; on each piece we use either sqrt(D) or i*sqrt(-D) (an analytic
    square root of D whenever the piece's D-image avoids (-inf,0] resp.
    [0,+inf), CERTIFIED per piece by ball evaluation), and propagate the
    overall sign across nodes with a certified matching test; the initial
    sign (fixed by branch_sign relative to the PRINCIPAL sqrt) is converted
    to the first piece's variant by the same certified matching test.
  * w_anti = 2*i*RF(0, e2-e1, e3-e1) (Carlson RF, certified) for the
    Weierstrass model y^2 + y = x^3 - x^2 (11.a3), same period lattice as the
    quartic model (birational over Q, invariant differential dx/u, kappa = 1;
    see notes/proof-n1.md).  Roots of 4x^3 - 4x^2 + 1 are isolated by Newton
    iteration + a Rouche certificate.

Output: P = I_signed + A_out - A_in as a certified ball, w_anti ditto, the
ratio P/w_anti, and the verdict that the ratio ball contains an integer n0
with |ratio - n0| < 1/2  ==>  period(C')/w_anti = n0 exactly (a priori
integer).  Expected |n0| = 2 (sign = orientation convention).
"""
from flint import acb, arb, acb_poly, ctx
import sys

ctx.prec = 300
# NOTE: do NOT set ctx.dps here -- in python-flint >= 0.9 setting dps
# overrides prec (dps = 50 would silently drop the working precision to
# ~169 bits).  The certification below is meant to run at 300 bits
# (~90 digits); the tolerances TOL_REL/TOL_ABS are unchanged.

PI = arb.pi()
C2 = PI / 2                        # c = pi/2
SQRT_C = C2.sqrt()

TOL_REL = 1e-13
TOL_ABS = 1e-25

call_count = 0


def D_and_x(th):
    """th: acb -> (x = e^{i th}, D = x^4 - 4x^3 + 2x^2 + 1)."""
    x = acb.exp(acb(0, 1) * th)
    return x, x**4 - 4 * x**3 + 2 * x**2 + 1


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
    |y| > 1 (big) resp. < 1 (small)."""
    x, D = D_and_x(acb(th_mid))
    B = x * x + 1
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
    th = th_fn(t).  The arc must contain no branch point of u."""
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
    # (For k=0 the first piece always turned out to be 'p', so this latent
    # bug never fired; ported back from k1_interval.py, referee item.)
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
    roots = [certify_root(p, dp, refine_newton(p, dp, s)) for s in starts]
    # completeness certificate (referee item): the certified root balls are
    # (a) pairwise disjoint -- certified distance between the discs > 0,
    #     i.e. centre distance > sum of radii (each ball contains ONE root),
    # (b) exhaustive -- exactly deg many, so every root is accounted for.
    # Both the tip-analyticity bound (zeros of D outside |t| <= rho) and the
    # Carlson RF root identification (e1, e2, e3 ordering) rely on this.
    deg = len(coeffs_asc) - 1
    assert len(roots) == deg, \
        "root count %d != polynomial degree %d" % (len(roots), deg)
    for i in range(deg):
        for j in range(i + 1, deg):
            assert abs(roots[i] - roots[j]).lower() > 0, \
                "certified root balls %d and %d are not disjoint" % (i, j)
    return roots


# ---------------------------------------------------------------------------
# analytic tip machinery

def tip_setup():
    """Choose (rho, delta) and certify max_{|t|=rho} |f| <= M.

    f(t) = 2 i t x(t^2) / u(t^2) is analytic on |t| <= rho because
    D(t^2) = t^2 w(t^2) with w nonzero for |t^2| <= rho^2; the zeros of
    D(theta) are theta = -i Log x_j + 2 pi i k for the roots x_j of the
    quartic, all certified outside |theta| <= rho^2.
    Returns (rho, delta, M).
    """
    roots = certify_poly_roots([1, 0, 2, -4, 1])     # D(z)
    best = None
    for xj in roots:
        if abs(xj - 1).upper() < 1e-20:
            continue                                # the root x = 1 (theta = 0)
        lj = xj.log()
        for k in (-1, 0, 1):
            th = lj + 2 * acb(0, 1) * arb.pi() * k  # candidate -i log x_j...
            d = abs(acb(0, 1) * th).lower()         # |theta_j,k|
            if best is None or d < best:
                best = d
    assert best is not None and best > 0
    rho = (arb(7) / 10 * best).sqrt()
    if rho.upper() > 0.6:
        rho = arb("0.6")
    delta = rho / 5
    # certify M = max |f| on |t| = rho by covering the circle with balls;
    # |f| = 2|t| / sqrt(|D|) is branch-independent.
    N = 4096
    M = arb(0)
    twor = 2 * arb.pi() * rho / N
    for k in range(N):
        tc = rho * acb.exp(acb(0, 1) * (2 * arb.pi() * k / N))
        tb = tc + square_ball(arb(13) / 10 * twor / 2)
        _, Db = D_and_x(tb * tb)
        dlow = abs(Db).lower()
        assert dlow > 0, "D vanishes on covering ball -- reduce rho"
        val = 2 * abs(tb).upper() / dlow.sqrt()
        if val > M:
            M = val
    return rho, delta, M


def tip_ball(sign_theta, want_big, rho, delta, M):
    """Certified ball for int_0^delta f(t) dt on the given branch/side."""
    sig = branch_sign(arb("0.05") * sign_theta, want_big)  # sheet on that side
    a0 = acb(0, 1) * acb.exp(acb(0, 1) * sign_theta * arb.pi() / 4) / sig
    H = M / (rho**2 * (1 - (delta / rho)**2))
    rad = H * delta**3 / 3
    # self-check (referee item M3): a0 comes from a hand derivation and was
    # injected without machine certification; a global sign/branch error in
    # a0 would shift each tip by 2*|a0|*delta = 0.24 and could still pass
    # the final 1/2 threshold.  f(t) = a0 + t^2 h(t) with |h| <= H on
    # |t| <= delta, so a certified point evaluation of f at t = delta must
    # satisfy |f(delta) - a0| <= H*delta^2; a sign flip (distance
    # 2*|a0| = 2 >> H*delta^2) is necessarily caught.
    x, D = D_and_x(acb(sign_theta * delta * delta))
    u = D.sqrt()
    if sig == -1:
        u = -u
    f_del = 2 * acb(0, 1) * delta * x / u
    dev = abs(f_del - a0).upper()
    bound = H * delta**2
    assert dev <= bound, \
        "tip self-check failed: |f(delta)-a0| = %s > H*delta^2 = %s" \
        % (dev, bound)
    print("  tip self-check (side=%+d, %-5s): |f(delta)-a0| <= %s <= H*delta^2 = %s"
          % (sign_theta, "big" if want_big else "small", dev, bound))
    return a0 * delta + square_ball(rad), rad


def inner_half(sign_theta, want_big, tip):
    """int over the half-arc touching th=0, via th = sign_theta * t^2."""
    rho, delta, M = tip
    body = arc_integral(delta, SQRT_C,
                        lambda t: sign_theta * t * t,
                        lambda t: 2 * t,
                        want_big)
    tb, trad = tip_ball(sign_theta, want_big, rho, delta, M)
    return body + tb, trad


# ---------------------------------------------------------------------------

def certify_w_anti():
    starts = [-0.4196433776070806,
              complex(0.7098216888035403, 0.3031453646035997),
              complex(0.7098216888035403, -0.3031453646035997)]
    roots = certify_poly_roots([1, 0, -4, 4], starts)   # 4x^3 - 4x^2 + 1
    e1 = roots[0]
    # referee item (alignment with k1_interval.py): e1 must be REAL
    assert abs(e1.imag).upper() < arb(10) ** (-30), "e1 not real?!"
    neg = roots[1] if roots[1].imag.mid() < 0 else roots[2]
    pos = roots[2] if roots[1].imag.mid() < 0 else roots[1]
    e2, e3 = neg, pos
    w_real = 2 * acb.elliptic_rf(acb(0), e1 - e2, e1 - e3)
    w_imag = 2 * acb(0, 1) * acb.elliptic_rf(acb(0), e2 - e1, e3 - e1)
    return w_real, w_imag


def main():
    print("=== interval certification, S_0 (k=0) ===")
    sys.stdout.flush()

    print("[0/4] analytic-tip setup (certified quartic roots, rho, M)...")
    rho, delta, M = tip_setup()
    print("  rho =", rho.mid(), " delta =", delta.mid(), " M =", M)
    tip = (rho, delta, M)
    sys.stdout.flush()

    print("[1/4] inner halves, big branch (th = +/- t^2):")
    ib_pos, t1 = inner_half(+1, True, tip)
    ib_neg, t2 = inner_half(-1, True, tip)
    inner_b = ib_pos + ib_neg
    print("  inner_b  =", inner_b, " tip radii:", t1.mid(), t2.mid())
    sys.stdout.flush()

    print("[2/4] outer arcs, big branch (cut crossing handled):")
    ob_pos = arc_integral(C2, PI, lambda t: t, lambda t: 1, True)
    print("  [c,pi]   =", ob_pos)
    ob_neg = arc_integral(-PI, -C2, lambda t: t, lambda t: 1, True)
    print("  [-pi,-c] =", ob_neg)
    I_signed = inner_b - (ob_pos + ob_neg)
    print("  I_signed =", I_signed)
    sys.stdout.flush()

    print("[3/4] small branch (A_in, A_out):")
    ai_pos, t3 = inner_half(+1, False, tip)
    ai_neg, t4 = inner_half(-1, False, tip)
    A_in = ai_pos + ai_neg
    A_out = arc_integral(C2, PI, lambda t: t, lambda t: 1, False) \
        + arc_integral(-PI, -C2, lambda t: t, lambda t: 1, False)
    print("  A_in     =", A_in, " tip radii:", t3.mid(), t4.mid())
    print("  A_out    =", A_out)
    sys.stdout.flush()

    print("[4/4] w_anti via Carlson RF (certified roots):")
    w_real, w_imag = certify_w_anti()
    print("  w_real   =", w_real, "  (expect 6.3460465213977671...)")
    print("  w_imag   =", w_imag, "  (expect i*2.9176332338769...)")
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
    print("cross-check: 2*w_anti (mid) =", (2 * w_anti).mid())
    print("our P (mid)               =", P.mid())
    print("mpmath 80-digit gave |period(C')/w_anti| = 1.99999999999999991...")


if __name__ == "__main__":
    main()
