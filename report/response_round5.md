# Response to the technical note after `rev4` (round 5)

We thank the referee for the careful reading of the new Siegel-unit anchor.
All three main points were valid; each is now repaired by an exact argument,
in every case along one of the routes the referee suggested. Two findings made
during the repair are worth flagging up front:

- On point 1, the premise of the "preferred" route fails: Brunault's Lemma 11
  states that `e_{a,b}` has weight 1 and level **N²** (i.e. `e_{a,b} ∈
  M₁(Γ₁(121))`, not `M₁(Γ₁(11))`), so no clean finite relabeling check under
  `Γ₀(11)` exists. We therefore carried out the referee's **fallback** route,
  which gives an unconditional result.
- On point 3, the exact computation revealed that our archive's claim
  "D_U = D_V = 0" was **wrong for U**: in fact `D_U = 2π` (U has winding
  number 1 along the cycle). The correction nevertheless vanishes, because the
  constants `C_x = −1`, `C_y = +1` are now determined **exactly**, so
  `log|C_x| = log|C_y| = 0` kills `log|C_x|·D_V − log|C_y|·D_U = 0·0 − 0·2π`.
  We have added a prominent erratum to `notes/attack16-siegel-anchor.txt`
  (line 72 and item 6 of the lessons list). The conclusion of the anchor is
  unchanged, but the reason is now correct — and stronger than before.

## Point 1: membership of F_total before Sturm

Agreed — the 251-coefficient comparison did not by itself prove membership in
`M₂(Γ₀(11))`, and the previous inference was circular. Repair (new
`code/siegel_anchor_step12.py`, archive `notes/attack17-membership.txt`):

1. By Brunault's Definition 10 and Lemma 11, every `e_{a,b}` occurring in
   `F_total` is an Eisenstein series of weight 1 on `Γ₁(121)`, holomorphic on
   ℋ and at every cusp. Hence the **unconditional** membership
   `F_total ∈ M₂(Γ₁(121))`, and `D := F_total + 2f₁₁ ∈ M₂(Γ₁(121))`.
2. The Sturm bound for this space is `(2/12)·[PSL₂(ℤ):Γ̄₁(121)] = (2/12)·7260
   = 1210` (−I acts trivially in even weight). We computed the q-expansion of
   `D` in exact rational arithmetic (higher e-coefficients are integral,
   `α₀(a,b) ∈ (1/22)ℤ`; Kronecker-packed exact convolution, cross-validated
   against the earlier 251 coefficients) and found `a_n(D) = 0` for all
   `0 ≤ n ≤ 2420` — twice the sharp bound and also beyond the conservative
   SL₂-index convention 2420.
3. By Sturm's theorem `D = 0`, i.e. `F_total = −2f₁₁` as modular forms; a
   fortiori `F_total ∈ M₂(Γ₀(11))`, so the old `Sturm bound 2` remark is
   retroactively legitimized (and removed from the text).

The stronger per-symbol claims are handled the same way: each of the seven
per-symbol forms lies in `M₂(Γ₁(121))`, and the identical computation
certifies that the two named symbols contribute exactly `−f₁₁` each and the
other five vanish identically (checked through q^2420). The main proof now
uses only the total; the per-symbol decomposition is stated as supplementary
but is likewise proved, not merely observed.

## Point 2: primitivity of the cycle

Repaired exactly as the referee's first suggestion (new
`code/siegel_anchor_step13.py`, archive `notes/attack17-primitivity.txt`).
We compute the integral homology of X₁(11) by Manin symbols for ±Γ₁(11)
(index 60 in PSL₂(ℤ), parametrized by bottom rows `(c,d) ∈
((ℤ/11)²∖0)/±1`), relations `x + xS = 0`, `x + xR + xR² = 0`, boundary map to
the 10 cusps, and Smith normal forms over ℤ:

- `H₁(X₁(11),ℤ) = ker ∂ ≅ ℤ²`, torsion-free, with explicit integral basis;
- conjugation `(c,d) ↦ (−c,d)` induces `C = [[0,1],[1,0]]` on H₁, so
  `H₁⁻ = ker(C+I) = ℤ·(−1,1)` has rank exactly one;
- the displayed seven-symbol chain telescopes to `[3/11] − [8/11] = 0` (the
  two cusps coincide for ±Γ₁(11)), so δ is closed; its coordinates in the
  integral basis are **(1,−1)**, i.e. δ is exactly ±1 times a generator of
  `H₁⁻` (equivalently, intersection number ±1 with an integral class).

Sanity checks in the same computation: 60 cosets, 10 cusps
(= ½∑φ(d)φ(11/d)), rank M = 11 = 2g + #cusps − 1. The 60-digit period
agreement and PARI's exact `v⁻(δ) = 1` (msfromell/mseval) are retained as
checks only, per the referee's remark.

## Point 3: exact constants and the cusp table

Repaired via both of the referee's routes combined (new
`code/siegel_anchor_step14.py`, archive `notes/attack17-constants.txt`).

(a) Argument periods. By Brunault's Lemma 5 (whose exact statement we verified
against the article's LaTeX source), each piece of each of the seven Manin
symbols contributes an explicit rational multiple of 2π (the root-of-unity
factors in `g_{a,b}∘γ = w·g_{(a,b)γ}` are constant and invisible to d arg).
Exact summation gives `D_U = 2π`, `D_V = 0` — see the erratum note above.

(b) Constants determined exactly. At cusps where both sides of a presentation
have order 0 one may evaluate directly: `π(4/11) = 4A`, `π(3/11) = A`, and
`x(4A) = 1`, `y(A) = −1`. The leading q-coefficients `κ_c(U)`, `κ_c(V)` are
explicit roots of unity, computed exactly by writing `γi∞ = c` as a word in
S, T and applying Brunault's Lemma 4 and eq. (3); the exact cyclotomic
computation gives `κ_{4/11}(U) = κ_{3/11}(V) = −1`, hence `C_x = 1/(−1) = −1`
and `C_y = (−1)/(−1) = +1` — exact equalities, not 70-digit observations.
Then `log|C_x| = log|C_y| = 0` exactly, and the correction
`log|C_x|·D_V − log|C_y|·D_U` vanishes despite `D_U = 2π`.

(c) Cusp–torsion table. The tuple `(m₁,…,m₅) = (0,2,1,4,3)` is now derived in
the text from Brunault's thesis table (3.152)–(3.153) (after Lecacheux):
`P₁ = ∞, P₂ = (1,0), P₃ = (0,−1), P₄ = (0,0), P₅ = (1,−1)`, `P_{4^a} = a·P₄`,
read against the group law `2A = (1,−1)`, `3A = (1,0)`, `4A = (0,−1)`. The
label conversion is stated explicitly: `P_v` is the cusp `k/11` iff
`v ≡ k^{−1} (mod 11)` up to ±1 — the anchor `P₁ = i∞ = 1/11` (since
`[[1,0],[11,1]] ∈ Γ₁(11)`) and diamond-operator transitivity restrict the
dictionary to `v ≡ k^{±1}`, and the exact divisor-order match selects the
inverse convention. Independently, the Kubert–Lang orders of U, V matched
against `div(x)`, `div(y)` force the same tuple uniquely. The 60-digit
Abel–Jacobi evaluation is described only as a check.

## Small items

- `step10` is now listed in the code table as a discarded, defective
  experiment and excluded from the proof/reproduction sequence; the table
  entry reads `step1–9, 11` plus the new `step12/13/14`.
- The remaining "three independent tempered symbols" occurrence has been
  changed to "three different tempered symbols" (this time verified by grep).
- "251 coefficients prove" no longer occurs: the membership argument of
  point 1 is supplied, and the text cites the 2420-coefficient exact Sturm
  computation instead.

## Verification summary

- `step12.py`: F_total = −2f₁₁ exact for n = 0..2420, all 7 per-symbol
  identities PASS, cross-check vs the old 251 coefficients PASS (11 s).
- `step13.py`: 17 checks PASS; δ has coordinates (1,−1), generates H₁⁻;
  PARI cross-check v⁻(δ) = 1 exact.
- `step14.py`: all exact items PASS (D_V = 0, D_U = 2π ∈ 2πℤ, endpoint-argument
  consistency, group law, cusp table via two independent exact routes,
  C_x = −1, C_y = +1, correction term exactly 0); numerical cross-checks to
  ~1e-60.
- `paper.tex` compiles cleanly (34 pages, no undefined references).

The three lattice/space/constant identifications are now proofs.
