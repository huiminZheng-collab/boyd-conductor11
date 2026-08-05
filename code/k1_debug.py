"""Debug: per-piece kinds/signs for the outer arcs, k=1."""
from flint import acb, arb, ctx
ctx.prec = 300
ctx.dps = 50
PI = arb.pi()
C2 = 2 * PI / 3

def D_and_x(th):
    x = acb.exp(acb(0, 1) * th)
    return x, x**4 - 2 * x**3 + 3 * x**2 + 2 * x + 1

def variant_sqrt(D, kind):
    if kind == 'p':
        return D.sqrt()
    return acb(0, 1) * (-D).sqrt()

def cut_avoidance(Dball):
    re, im = Dball.real, Dball.imag
    if (not im.contains(arb(0))) or re.lower() > 0:
        return 'p'
    if re.upper() < 0:
        return 'r'
    return None

def ball_interval(lo, hi):
    lo, hi = arb(lo), arb(hi)
    return arb((lo + hi) / 2, (hi - lo) / 2)

def branch_sign(th_mid, want_big):
    x, D = D_and_x(acb(th_mid))
    B = x * x + x + 1
    s = D.sqrt()
    out = {}
    for sigma in (1, -1):
        u = s if sigma == 1 else -s
        y = (-B + u) / 2
        out[sigma] = abs(y)
    m1, m_1 = out[1], out[-1]
    print("   branch_sign at th=%s: |y(+)|=%s |y(-)|=%s" % (th_mid, m1, m_1))
    if (m1 - 1).contains(arb(0)) or (m_1 - 1).contains(arb(0)):
        raise RuntimeError("inconclusive")
    big_is_plus = m1.lower() > 1
    return 1 if (big_is_plus == want_big) else -1

def probe(t_lo, t_hi, want_big, n=8):
    t_lo, t_hi = arb(t_lo), arb(t_hi)
    nodes = [t_lo + (t_hi - t_lo) * i / n for i in range(n + 1)]
    kinds = []
    for j in range(n):
        tb = acb(ball_interval(nodes[j], nodes[j + 1]))
        _, Db = D_and_x(tb)
        kinds.append(cut_avoidance(Db))
    print("kinds:", kinds)
    sig = branch_sign((nodes[0] + nodes[1]) / 2, want_big)
    print("initial sig:", sig)
    for j in range(n - 1):
        _, Dn = D_and_x(acb(nodes[j + 1]))
        u_arr = variant_sqrt(Dn, kinds[j])
        if sig == -1:
            u_arr = -u_arr
        nkind = kinds[j + 1]
        sep = abs(variant_sqrt(Dn, nkind)).lower()
        best = None
        for cand in (1, -1):
            u_c = variant_sqrt(Dn, nkind)
            if cand == -1:
                u_c = -u_c
            if abs(u_c - u_arr).upper() < sep:
                best = cand
        print("node %s: kind %s->%s sig %s->%s (sep=%s)" %
              (nodes[j+1], kinds[j], nkind, sig, best, sep))
        sig = best

print("=== [c, pi] big ===")
probe(C2, PI, True)
print("=== [-pi, -c] big ===")
probe(-PI, -C2, True)
print("=== [-c, c] big ===")
probe(-C2, C2, True)
