"""
Exact diamond-product computation for S_1 (k=1), cross-check of hand work.

Cubic model S_1: y^2 + (x^2+x+1)y + x^3 = 0, identity O = [0:1:0],
A = (0,0) has exact order 4 (chord-tangent, notes/proof-k1.md):
  2A = (0,-1) = Q,  3A = T = [1:-1:0],  4A = O.
Divisors (orders verified by local expansion, notes/proof-k1.md):
  div(x) = [A] + [2A] - [O] - [3A]
  div(y) = 3[A] - 2[O] - [3A]
Diamond: (x)<>(y) = sum_ij m_i n_j [P_i - Q_j] in Z[E], points in Z/4.
Convention: point k means k*A, k in Z/4; O = 0.
Then evaluate D_E:  D_E(O) = D_E(2A) = 0 (2-torsion), D_E(3A) = -D_E(A).

Also do k=0 (Z/5) as a control: must reproduce 6(O)+5(A)-5(2A).
"""
from collections import defaultdict

def diamond(divf, divg, n):
    """divf, divg: dict point(coeff in Z/nA, O=0) -> multiplicity.
    Returns dict of raw coefficients for [point] in Z[E]."""
    raw = defaultdict(int)
    for P, m in divf.items():
        for Q, k in divg.items():
            raw[(P - Q) % n] += m * k
    return dict(raw)

# ---- k=1, Z/4 ----
divx = {1: 1, 2: 1, 0: -1, 3: -1}
divy = {1: 3, 0: -2, 3: -1}
raw = diamond(divx, divy, 4)
print("k=1 raw diamond:", raw)
# anti-symmetrize in Z[E]^- : [k] -> -[(-k) mod 4]
asym = defaultdict(int)
for P, c in raw.items():
    asym[P] += c
    asym[(-P) % 4] -= c
print("k=1 anti-symmetrized (divided by 2 in Q):",
      {k: v // 2 for k, v in asym.items()})
# D_E evaluation: D_E(0)=D_E(2)=0, D_E(3) = -D_E(1)
DEA = 0.4702266562812140308266032191876266124613
val = raw.get(1, 0) * DEA + raw.get(3, 0) * (-DEA)
print("k=1 D_E(diamond) =", val, " = ", val / (3.1415926535897932385 * 0.29935558688291539), "pi b_17")

# ---- k=0 control, Z/5 ----
raw0 = diamond(divx, divy, 5)
print("k=0 raw diamond:", raw0)
asym0 = defaultdict(int)
for P, c in raw0.items():
    asym0[P] += c
    asym0[(-P) % 5] -= c
print("k=0 anti-symmetrized (expect 6(O)+5(A)-5(2A)):",
      {k: v // 2 for k, v in asym0.items()})
DEP0 = 0.1911937370843316957549544343121738161012
DE2P0 = 1.5 * DEP0
val0 = (raw0.get(1, 0) - raw0.get(4, 0)) * DEP0 + (raw0.get(2, 0) - raw0.get(3, 0)) * DE2P0
print("k=0 D_E(diamond) =", val0, " vs -pi b_11 =", -3.1415926535897932385 * 0.15214714172591805)
