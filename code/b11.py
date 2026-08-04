"""
Boyd conductor-11 attack — Step 1: high-precision b_11 = L'(E11,0) = 11/(4 pi^2) L(E11,2).

E11 = X_1(11), cusp form f(τ) = η(τ)^2 η(11τ)^2 = sum a_n q^n, weight 2, level 11, root number +1.

For weight 2, w=+1:
  Λ(f,s) = N^{s/2} (2π)^{-s} Γ(s) L(f,s),  Λ(f,s) = Λ(f,2-s).
  b_11 = L'(E,0) = Λ(f,0) = Λ(f,2)
       = Σ_n a_n [ e^{-t_n} (1/t_n + 1/t_n^2) + E_1(t_n) ],  t_n = 2πn/√11,
  since Λ(f,2) = Σ a_n ∫_1^∞ e^{-t_n y} (y + y^{-1}) dy  (∫_1^∞ e^{-ty}/y dy = E_1(t)).

Sanity check (Boyd, PNWNT 2015 slides): b_11 = 0.1521471...
"""
from mpmath import mp, mpf, exp, sqrt, pi

def coefficients_f11(nmax):
    """a_1..a_nmax of q * prod (1-q^n)^2 (1-q^{11n})^2 as exact integers."""
    # A(q) = prod_{n>=1} (1-q^n)^2  (Euler function squared), truncated to nmax
    A = [0] * (nmax + 1)
    A[0] = 1
    for n in range(1, nmax + 1):
        for j in range(nmax, n - 1, -1):
            A[j] -= 2 * A[j - n]
            if j >= 2 * n:
                A[j] += A[j - 2 * n]
    # B(q) = prod (1-q^{11n})^2
    B = [0] * (nmax + 1)
    B[0] = 1
    for n in range(1, nmax // 11 + 1):
        m = 11 * n
        for j in range(nmax, m - 1, -1):
            B[j] -= 2 * B[j - m]
            if j >= 2 * m:
                B[j] += B[j - 2 * m]
    # f/q = A*B ; a_{n} = (A*B)[n-1]
    C = [0] * (nmax + 1)
    for i in range(nmax + 1):
        if A[i]:
            for j in range(nmax + 1 - i):
                if B[j]:
                    C[i + j] += A[i] * B[j]
    return [C[n - 1] for n in range(1, nmax + 1)]

def b11(dps=80, nmax=None):
    mp.dps = dps
    if nmax is None:
        # e^{-2πn/√11}: decay ~ 10^{-0.8225 n}; take enough terms
        nmax = int(dps * 0.4343 / 0.8225) + 20
    a = coefficients_f11(nmax)
    t = 2 * pi / sqrt(11)
    s = mpf(0)
    for n in range(1, nmax + 1):
        tn = t * n
        s += a[n - 1] * (exp(-tn) * (1 / tn + 1 / tn**2) + mp.expint(1, tn))
    return s

if __name__ == "__main__":
    mp.dps = 80
    b = b11(80)
    print("b_11 =", mp.nstr(b, 70))
    print("L(E11,2) =", mp.nstr(4 * pi**2 / 11 * b, 70))
    # quick check against Boyd's value
    print("matches 0.1521471... :", mp.nstr(b, 10))
