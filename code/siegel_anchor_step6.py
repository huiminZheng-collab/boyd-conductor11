# Siegel-unit representation of x o pi, y o pi on X_1(11) via Kubert-Lang linear algebra.
# Cusps of X(11): first columns (r,t) in P^1(Z/11)/+- (60 cusps, all width 11).
# Siegel units g_{a,b}: (a,b) in ((Z/11)^2 - 0)/+- (60 units).
# div(g_{a,b}) at cusp (r,t): (11/2) B2({(a r + b t)/11})   [order in local param q^{1/11}]
# Target: div(u o proj) at cusp over rational cusp k/11 (= class (r,0), r=k mod +-):
#         11 * ord_{k/11}^{X_1(11)}(u), and 0 over non-rational cusps.
# Dividing by 11:  sum n_{a,b} (1/2) B2({(ar+bt)/11}) = ord_{k/11}(u) or 0.
from fractions import Fraction as Fr
from sympy import Matrix

N = 11

def B2h(j):   # (1/2) B_2(j/11)
    x = Fr(j, N)
    return (x*x - x + Fr(1, 6)) / 2

cusps = [(r, t) for t in range(1, 6) for r in range(11)] + [(r, 0) for r in range(1, 6)]
units = [(a, b) for a in range(1, 6) for b in range(11)] + [(0, b) for b in range(1, 6)]
assert len(cusps) == 60 and len(units) == 60

M = Matrix([[B2h((a*r + b*t) % N) for (a, b) in units] for (r, t) in cusps])

# torsion orders from Abel computation (step4): pi(k/11) = m_k A, m = (0,2,1,4,3), k=(1..5)
# div(x_E) = [A]+[2A]-[O]-[3A];  div(y_E) = 3[2A]-2[3A]-[O]
mx = {0: -1, 1: 1, 2: 1, 3: -1, 4: 0}
my = {0: -1, 1: 0, 2: 3, 3: -2, 4: 0}
m_of_k = {1: 0, 2: 2, 3: 1, 4: 4, 5: 3}

def target(ords):
    return Matrix([Fr(ords[m_of_k[r]]) if t == 0 else Fr(0) for (r, t) in cusps])

def solve(T):
    aug = M.row_join(T)
    R, piv = aug.rref()
    # check consistency
    for i in range(60):
        if all(R[i, j] == 0 for j in range(60)):
            assert R[i, 60] == 0, "system inconsistent!"
    # particular solution: free vars = 0
    sol = [Fr(0)]*60
    for i, p in enumerate(piv):
        if p < 60:
            sol[p] = R[i, 60]
    return sol, len(piv)

nx, rank = solve(target(mx))
ny, _ = solve(target(my))
print("rank(M) =", rank)
suppx = [(units[i], nx[i]) for i in range(60) if nx[i] != 0]
suppy = [(units[i], ny[i]) for i in range(60) if ny[i] != 0]
print("support x:", len(suppx))
for ab, n in suppx: print("   ", ab, n)
print("support y:", len(suppy))
for ab, n in suppy: print("   ", ab, n)

# kernel of M (should be 1-dim)
ker = M.nullspace()
print("kernel dim:", len(ker))
if ker:
    k0 = ker[0]
    print("kernel vector (nonzero entries):")
    for i in range(60):
        if k0[i] != 0:
            print("   ", units[i], k0[i])

import json
with open("code/siegel_anchor_step6_out.json", "w") as f:
    json.dump({
        "nx": [[list(ab), str(n)] for ab, n in suppx],
        "ny": [[list(ab), str(n)] for ab, n in suppy],
        "kernel": [[[list(units[i])], str(k0[i])] for i in range(60) if ker and k0[i] != 0] if ker else [],
    }, f, indent=1)
print("written code/siegel_anchor_step6_out.json")
