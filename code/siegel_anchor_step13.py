# attack17 / step13:  PRIMITIVITY of the anti-invariant cycle
#   delta = {0,3/11} - {0,8/11}  on  X_1(11) = 11.a3.
#
# Pure integer linear algebra (Manin symbols for +-Gamma_1(11), index 60 in
# PSL_2(Z)), following Manin / Stein "Modular Forms: A Computational
# Approach" sec.3:
#   * cosets  +-Gamma_1(11)\PSL_2(Z)  <->  nonzero (c,d) in (Z/11)^2 / {+-1}
#     (bottom row of a matrix, mod 11, up to sign):  60 cosets.
#   * Manin symbols x_i = {g_i oo, g_i 0}, relations
#       x + xS = 0,        S = (0 -1; 1 0)
#       x + xR + xR^2 = 0, R = (0 -1; 1 1)   (order 3 in PSL_2(Z)).
#   * boundary d(x_i) = [g_i 0] - [g_i oo] to the cusp divisor group
#     (cusps = double cosets +-Gamma_1(11)\PSL_2(Z)/<T>, T = (1 1; 0 1)).
#   * H_1 = ker d, integral basis via Smith normal form (with transforms).
#   * complex conjugation c: tau |-> -conj(tau) acts on bottom rows by
#     (c,d) |-> (-c,d)  (matrix J g J, J = diag(1,-1)).
# We then express the 7-symbol chain of delta (from
# notes/attack16-siegel-anchor.txt step 8-9) in the integral H_1 basis and
# check that its coordinate vector is primitive (gcd = 1) and
# anti-invariant.  Since H_1(X_1(11),Z)^- has rank 1, primitivity is
# exactly "delta generates H_1^-".
#
# Cross-check (optional, needs gp in PATH): PARI mseval/mfsymboleval gives
#   int_delta omega = w1 - 2 w2 = 2.9176332338769904... * I  (= w_anti)
# and the normalized anti-invariant PARI symbol evaluates to exactly 1 on
# delta.
#
# Run:  .venv/Scripts/python.exe code/siegel_anchor_step13.py
# Output: console PASS/FAIL summary + code/siegel_anchor_step13_out.json

import json
import math
import os
import subprocess

from sympy import Matrix
from sympy.matrices.normalforms import smith_normal_decomp

N = 11
HERE = os.path.dirname(os.path.abspath(__file__))
OUT_JSON = os.path.join(HERE, "siegel_anchor_step13_out.json")

# ---------------------------------------------------------------- cosets
def canon(c, d):
    """Canonical representative of +-(c,d) mod N, (c,d) != (0,0)."""
    c %= N
    d %= N
    assert (c, d) != (0, 0)
    if c > 5 or (c == 0 and d > 5):
        c = (-c) % N
        d = (-d) % N
    return (c, d)

COSETS = sorted({canon(c, d) for c in range(N) for d in range(N) if (c, d) != (0, 0)})
assert len(COSETS) == 60
IDX = {cd: i for i, cd in enumerate(COSETS)}

def mul_cd(cd, M):
    """Right action of 2x2 matrix M on bottom row (c,d), then canon."""
    (c, d) = cd
    ((p, q), (r, s)) = M
    return canon(c * p + d * r, c * q + d * s)

S = ((0, -1), (1, 0))
R = ((0, -1), (1, 1))
T = ((1, 1), (0, 1))
J = ((-1, 0), (0, 1))   # (c,d) -> (-c,d): complex conjugation

PERM_S = [IDX[mul_cd(cd, S)] for cd in COSETS]
PERM_R = [IDX[mul_cd(cd, R)] for cd in COSETS]
PERM_T = [IDX[mul_cd(cd, T)] for cd in COSETS]
PERM_J = [IDX[mul_cd(cd, J)] for cd in COSETS]

# ---------------------------------------------------------------- cusps
# cusps = double cosets +-Gamma_1(11)\PSL_2(Z)/<T>: orbits under right T.
cusp_of = list(range(60))
for i in range(60):
    j = PERM_T[i]
    while True:
        if cusp_of[j] < cusp_of[i]:
            cusp_of[i] = cusp_of[j]
        j = PERM_T[j]
        if j == i:
            break
# flatten + relabel
for i in range(60):
    while cusp_of[i] != cusp_of[cusp_of[i]]:
        cusp_of[i] = cusp_of[cusp_of[i]]
CUSP_IDS = sorted(set(cusp_of))
NCUSP = len(CUSP_IDS)
cusp_label = {c: k for k, c in enumerate(CUSP_IDS)}
cusp_of = [cusp_label[c] for c in cusp_of]

# a point a/c of P^1(Q) -> cusp label: bottom-row coset of any matrix with
# second column... rather: cusp of g*oo has bottom row of g; cusp a/c =
# cusp of matrix (a b; c d) -> coset (c,d).
def cusp_of_point(a, c):
    return cusp_of[IDX[canon(c, a)]]
    # matrix (a b; c d) has bottom row (c,d); gcd(a,c)=1 => (c,d)!=(0,0) mod 11

# ---------------------------------------------------------------- boundary
# d(x_i) = [g_i 0] - [g_i oo];  g 0 = gS oo  =>  d(e_i) = cusp(S i) - cusp(i)
B = Matrix.zeros(NCUSP, 60)
for i in range(60):
    B[cusp_of[i], i] -= 1
    B[cusp_of[PERM_S[i]], i] += 1

# ---------------------------------------------------------------- relations
def e(i):
    v = [0] * 60
    v[i] = 1
    return v

rel_rows = set()
for i in range(60):
    rel_rows.add(tuple([a + b for a, b in zip(e(i), e(PERM_S[i]))]))
    j, k = PERM_R[i], PERM_R[PERM_R[i]]
    rel_rows.add(tuple([a + b + c for a, b, c in zip(e(i), e(j), e(k))]))
REL = [list(r) for r in sorted(rel_rows)]

# relations must be closed (lie in ker d)
rel_closed = all((B * Matrix(r)).is_zero_matrix for r in REL)

# ---------------------------------------------------------------- H_1
# K = ker(d) : SNF of B, U B V = D; kernel = columns of V beyond rank.
Db, Ub, Vb = smith_normal_decomp(B)
rankB = sum(1 for i in range(min(Db.rows, Db.cols)) if Db[i, i] != 0)
Kmat = Vb[:, rankB:]                      # 60 x (60-rankB), columns = Z-basis of ker d
rankK = Kmat.cols

def coords_in_K(v):
    """v: 60-dim integer vector in ker d -> coords in K-basis (exact)."""
    x, _params = Kmat.gauss_jordan_solve(Matrix(v))
    assert Kmat * x == Matrix(v)
    return x

# relation sublattice Lambda inside K (coords), as columns of Acol
Acol = Matrix.hstack(*[coords_in_K(r) for r in REL])
rankA = Acol.rank()

# H_1 = K / Lambda : SNF of Acol, U Acol V = D
Da, Ua, Va = smith_normal_decomp(Acol)
diag_a = [Da[i, i] for i in range(min(Da.rows, Da.cols))]
nonzero_diag = [int(d) for d in diag_a if d != 0]
rankH1 = rankK - rankA
torsion_free = all(d == 1 for d in nonzero_diag)

# H_1 basis (in K-coords): columns of Ua^{-1} beyond rankA
Ua_inv = Ua.inv()
H1_basis_K = [Ua_inv[:, rankA + j] for j in range(rankH1)]       # in K-coords
H1_basis_60 = [Kmat * h for h in H1_basis_K]                     # as symbol vectors

def H1_coords(v):
    """v in ker d (60-dim) -> coordinate vector in H_1 basis (last rankH1
    entries of Ua * coords_in_K(v)); requires v's class well-defined."""
    y = Ua * coords_in_K(v)
    return [int(y[rankA + j]) for j in range(rankH1)], y

# ---------------------------------------------------------------- delta
# 7-symbol continued-fraction chain (attack16 step 8-9):
#   {0,3/11} = +[[1,0],[3,1]] - [[1,1],[3,4]] + [[3,1],[11,4]]
#   {0,8/11} = +[[1,0],[1,1]] - [[1,2],[1,3]] + [[3,2],[4,3]] - [[3,8],[4,11]]
# symbol [[a,b],[c,d]] = matrix (a b; c d) -> coset (c,d).
CHAIN_311 = [(((1, 0), (3, 1)), +1), (((1, 1), (3, 4)), -1), (((3, 1), (11, 4)), +1)]
CHAIN_811 = [(((1, 0), (1, 1)), +1), (((1, 2), (1, 3)), -1),
             (((3, 2), (4, 3)), +1), (((3, 8), (4, 11)), -1)]

def sym_index(mat):
    (a, b), (c, d) = mat
    return IDX[canon(c, d)]

delta = [0] * 60
chain_descr = []
for chain, chsgn in ((CHAIN_311, +1), (CHAIN_811, -1)):
    for mat, sgn in chain:
        delta[sym_index(mat)] += chsgn * sgn
        chain_descr.append({"matrix": [list(mat[0]), list(mat[1])],
                            "coeff": chsgn * sgn, "coset": list(COSETS[sym_index(mat)])})

delta_closed = (B * Matrix(delta)).is_zero_matrix

# telescoping sanity: d(chain(3/11)) = [3/11]-[0], d(chain(8/11)) = [8/11]-[0],
# and [3/11] = [8/11] as +-Gamma_1(11) cusps:
cusp_311 = cusp_of_point(3, 11)
cusp_811 = cusp_of_point(8, 11)

delta_H1, y_delta = H1_coords(delta)
# integrality of delta's class in H_1: y_delta[i] must be divisible by the
# i-th SNF invariant factor d_i of the relation lattice (all d_i = 1 here,
# so the class is automatically integral; the check is kept for honesty).
delta_class_integral = all(int(y_delta[i]) % int(diag_a[i]) == 0 for i in range(rankA))

# ------------------------------------------------- complex conjugation
# permutation matrix action on symbol vectors
def apply_J(v):
    w = [0] * 60
    for i in range(60):
        w[PERM_J[i]] = v[i]
    return w

# J preserves ker d -> 51x51 matrix Mc on K-coords
Mc = Matrix.hstack(*[coords_in_K(list(apply_J(list(Kmat[:, j])))) for j in range(rankK)])
# J preserves the relation lattice:  Mc*Acol = Acol*W
Wcols = []
J_preserves_rel = True
for j in range(Acol.cols):
    try:
        w, _p = Acol.gauss_jordan_solve(Mc * Acol[:, j])
        if Acol * w != Mc * Acol[:, j]:
            J_preserves_rel = False
            break
        Wcols.append(w)
    except Exception:
        J_preserves_rel = False
        break

# induced 2x2 matrix on H_1 = coker(Acol): bottom-right block of Ua Mc Ua^{-1}
Mc_bar = Ua * Mc * Ua_inv
Cmat = Matrix([[int(Mc_bar[rankA + i, rankA + j]) for j in range(rankH1)]
               for i in range(rankH1)])
lower_left_zero = all(Mc_bar[rankA + i, j] == 0
                      for i in range(rankH1) for j in range(rankA))

z_delta = Matrix(delta_H1)
anti_inv_delta = (Cmat * z_delta == -z_delta)
C_sq = (Cmat * Cmat == Matrix.eye(rankH1))

# anti-invariant sublattice ker(C+I): primitive generator
CpI = Cmat + Matrix.eye(rankH1)
ns = CpI.nullspace()
anti_rank = len(ns)
if anti_rank == 1:
    v_anti = ns[0]
    g = math.gcd(*[abs(int(x)) for x in v_anti])
    v_anti = Matrix([int(x) // g for x in v_anti])
    k = None
    for j in range(rankH1):
        if v_anti[j] != 0:
            k = int(z_delta[j]) // int(v_anti[j])
            break
    delta_is_generator = (z_delta == k * v_anti) and abs(k) == 1
else:
    v_anti = None
    k = None
    delta_is_generator = False

gcd_coords = math.gcd(*[abs(x) for x in delta_H1]) if any(delta_H1) else 0

# ------------------------------------------------- optional gp cross-check
gp_result = {"status": "skipped"}
try:
    script = r"""
default(realbitprecision, 128);
E = ellinit([0,-1,1,0,0]);
mfE = mffromell(E);
fs = mfsymbol(mfE[1], mfE[2]);
v = mfsymboleval(fs, [0, 3/11]) - mfsymboleval(fs, [0, 8/11]);
w = E.omega; wanti = w[1]-2*w[2];
printf("VINT %.30f\n", real(v));
printf("WANTI %.30f\n", imag(wanti));
printf("DEV %.5e\n", real(abs(2*Pi*I*v/wanti - 1)));
Mm = msfromell(E, -1);
vm = mseval(Mm[1], Mm[2], [0, 3/11]) - mseval(Mm[1], Mm[2], [0, 8/11]);
Mp = msfromell(E, 1);
vp = mseval(Mp[1], Mp[2], [0, 3/11]) - mseval(Mp[1], Mp[2], [0, 8/11]);
printf("VM %d\n", vm);
printf("VP %d\n", vp);
"""
    r = subprocess.run(["gp", "-q", "-f"], input=script, capture_output=True,
                       text=True, timeout=300)
    vals = {}
    for line in r.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 2:
            vals[parts[0]] = "".join(parts[1:])   # gp prints "1.2 E-30"
    dev = float(vals["DEV"])
    gp_result = {"status": "ok", "int_f_dtau_real": vals["VINT"],
                 "wanti_imag": vals["WANTI"], "vm_exact": int(vals["VM"]),
                 "vp_exact": int(vals["VP"]),
                 "abs_dev_ratio_minus_1": dev, "dev_lt_1e-25": dev < 1e-25}
except Exception as exc:  # gp missing or failed: cross-check only
    gp_result = {"status": "failed", "error": str(exc)}

# ---------------------------------------------------------------- report
checks = {
    "60 cosets": len(COSETS) == 60,
    "10 cusps (=(1/2)sum phi(d)phi(11/d))": NCUSP == 10,
    "rank(d) = #cusps - 1 = 9": rankB == NCUSP - 1,
    "relations closed (subset ker d)": rel_closed,
    "rank K = 60 - 9 = 51": rankK == 51,
    "rank relations = 49": rankA == 49,
    "H_1 rank = 2": rankH1 == 2,
    "H_1 torsion-free (all SNF diag = 1)": torsion_free,
    "cusp(3/11) = cusp(8/11)": cusp_311 == cusp_811,
    "delta closed in M": delta_closed,
    "delta class integral in H_1": delta_class_integral,
    "J preserves ker d and relations": J_preserves_rel and lower_left_zero,
    "C^2 = I on H_1": C_sq,
    "rank H_1^- = 1": anti_rank == 1,
    "delta anti-invariant in H_1": anti_inv_delta,
    "delta H_1 coords gcd = 1 (PRIMITIVE)": gcd_coords == 1,
    "delta = +- generator of H_1^-": bool(delta_is_generator),
}

print("=" * 72)
print("step13: primitivity of delta = {0,3/11} - {0,8/11} in H_1(X_1(11),Z)^-")
print("=" * 72)
print(f"#cosets (+-Gamma_1(11)\\PSL_2(Z)) : {len(COSETS)}")
print(f"#cusps                            : {NCUSP}")
print(f"#unique Manin relation rows       : {len(REL)}  (rank {rankA})")
print(f"rank M = 60 - rank(rel)           : {60 - rankA}  (expect 2g+#cusps-1 = 11)")
print(f"rank ker d                        : {rankK}")
print(f"rank H_1                          : {rankH1}   SNF diag of relations: {nonzero_diag[:5]}... (all 1: {torsion_free})")
print(f"cusp(3/11) == cusp(8/11)          : {cusp_311 == cusp_811} (cusp label {cusp_311})")
print(f"delta H_1 coords (in basis below) : {delta_H1}   gcd = {gcd_coords}")
print(f"H_1 basis vector 1 (60-dim)       : {list(H1_basis_60[0])}")
print(f"H_1 basis vector 2 (60-dim)       : {list(H1_basis_60[1])}")
print(f"conjugation C on H_1              : {[list(r) for r in Cmat.tolist()]}")
if v_anti is not None:
    print(f"H_1^- generator (basis coords)    : {list(v_anti)}   delta = {k} * generator")
print(f"gp cross-check                    : {gp_result}")
print("-" * 72)
allpass = True
for name, ok in checks.items():
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    allpass &= bool(ok)
print("-" * 72)
print("OVERALL:", "PASS - delta is a generator of H_1(X_1(11),Z)^-"
      if allpass else "FAIL - see above")

out = {
    "cosets_60": [list(cd) for cd in COSETS],
    "n_cusps": NCUSP,
    "cusp_of_coset": cusp_of,
    "rank_boundary": rankB,
    "rank_kernel": rankK,
    "rank_relations": rankA,
    "rank_H1": rankH1,
    "snf_diag_relations": nonzero_diag,
    "H1_basis_symbol_vectors": [list(map(int, h)) for h in H1_basis_60],
    "delta_chain": chain_descr,
    "delta_symbol_vector": delta,
    "delta_H1_coords": delta_H1,
    "delta_H1_coords_gcd": gcd_coords,
    "conjugation_matrix_on_H1": [list(map(int, r)) for r in Cmat.tolist()],
    "H1_minus_generator_coords": [int(x) for x in v_anti] if v_anti is not None else None,
    "delta_multiple_of_generator": k,
    "checks": {k2: bool(v2) for k2, v2 in checks.items()},
    "gp_cross_check": gp_result,
}
with open(OUT_JSON, "w") as f:
    json.dump(out, f, indent=1)
print("json written:", OUT_JSON)
