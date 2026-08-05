# Response to the second referee report

We thank the referee for the careful re-reading. All four essential items
and all four reproducibility items are addressed in revision `rev3`.
Below we answer point by point.

## Item 1 — exact torus-intersection criterion

Accepted. The displayed size condition in `prop:torus` was derived by an
incorrect modulus step. The proposition is restated in the two-case form
that the referee suggests: a torus point at parameter θ exists iff

- (i) θ = 0 and |k+2| ≤ 2, or
- (ii) k + 2cos θ = 0 (solvable iff |k| ≤ 2),

and the proof now compares imaginary and real parts of
`2cosψ·e^{iθ/2} = −(k+2cosθ)` directly. The final classification
(intersection iff k ∈ [−4,2], folds at ±arccos(−k/2), the θ=0
branch-exchange for −4 ≤ k ≤ 0, tangencies at k = −4, 0, and the merger at
k = 2) is unchanged; every downstream description (fold angles in
Table `tab:family`, the branch-exchange terminology, the tangency rows) was
rechecked against the corrected equations, and the classification was
additionally verified numerically on a k-grid (intersection scan, fold
moduli, and the gap outside [−4,2]).

## Item 2 — all-Arb sign certificate

Accepted. The `hull` branch of `code/sign_certify.py` no longer converts
Arb endpoints to binary floats. The convex hull is now formed entirely
inside ball arithmetic as

```
H = l1 + [0,1]·(l2 − l1)        (t = arb("0.5","0.5"))
```

which contains both input balls because ball subtraction contains all
pairwise differences; containment is asserted programmatically
(`assert H.contains(l1) and H.contains(l2)`). No empirical inflation
factor remains. The same audit was applied to `code/k1_sign_certify.py`,
whose algebraic hull (|y| ≤ M := (|B|+√(|B|²+4|x|³))/2 etc.) is now
computed entirely in `arb` as well, with the same containment assertion.
The duplicate `x = exp(iθ)` evaluation was removed and the `I1`/`I2`
output labels now state explicitly that the printed balls are the
unnormalised integrals. Both archived certificates
(`notes/attack14-sign-k0.txt`, `notes/attack14-sign-k1.txt`) were
regenerated: the certified conclusions are unchanged
(I_split ∈ [0.1489, 0.1553] ⊂ (0,∞); J_split strictly negative,
ñ(1) ∈ [0.2962, 0.3025]). Theorem `thm:cert`(6) now states that all range
enclosures, including the non-separated convex hulls, are formed and
verified within ball arithmetic.

## Item 3 — conductor-17 modular-unit assertion

Accepted; the assertion was indeed wrong on the natural modular model, and
it is removed. Proposition `thm:family`(1) now states: for k = 0 the
divisors of x, y are supported exactly on the rational cusps of X_1(11)
(modular units, proved); for k = 1 the divisors are supported on rational
torsion and the symbol is tempered (Newton face polynomials cyclotomic),
so {x,y} ∈ K₂(E)⊗ℚ — which is all the conductor-17 argument uses. The
proposition now says explicitly that we do *not* assert x, y to be
modular units for k = 1, and why (X_0(17) has two cusps; divisors
supported on all four rational torsion points cannot be cuspidal there).
The appendix status paragraph is rewritten to list precisely the two
hypotheses used — exact temperedness, and Bloch's theorem in the factor-1
normalization of [Lalín–Ramamonjisoa, Thm. 6] (conditional, via the
consistency argument of Remark `rem:normalisation`) — and states that no
modular-unit property is used or claimed. The abstract, the family
dichotomy remark, introduction item 5, and Table `tab:family` (k = 1 row
now "proved (cond.)") are updated accordingly.

## Item 4 — Brunault citation and exact anchor

Accepted.

- The citations are separated: the curve model and the D_E-coefficients
  are cited as [Brunault, Thm. 8 and Cor. 101, eq. (3.151)]; the
  regulator anchor is cited as [Brunault, Thm. 118 and its proof,
  eqs. (3.210)–(3.211)]. We verified in the thesis text that Théorème 118
  (§3.9) proves (3.208) and that its proof contains (3.210)–(3.211),
  reducing via (3.151) to Bertin's identity |r({x,y})| = (5/2π)D_E(P)
  [Bertin, Crelle 569, Th. 6 et Cor. 6.1 — Brunault's reference [10]].
- The ratio lemma is no longer conditional: the existence of a single
  symbol-independent constant c is now supplied by Bloch's theorem
  [Bloch] in the formulation of [Lalín–Ramamonjisoa, Thm. 6]
  (r({f,g})[γ] = D_E((f)⋄(g)), i.e. c = 1 in their normalization), with
  a note that any consistent change of conventions rescales both sides by
  fixed scalars; the ratio identity, the only thing used, is an immediate
  corollary. ℚ-linearity extends it from K₂(E) to K₂(E)⊗ℚ, hence to all
  tempered symbols.
- The model/function dictionary is now stated at the level of the two
  theorems, with verbatim quotes of the thesis sentences fixing the
  meaning of x, y (Thm. 8's equation y²+y = x³−x²; §3.9's opening
  sentence referring back to §3.7); the divisor-wise check (exact) and
  the 9-digit numerical exclusion are retained as independent
  corroborations, with the numerical item explicitly labelled a check,
  not part of the proof. The certified-computation theorem's closing
  statement is updated to match (Bloch + LalinRam Thm. 6 + Brunault
  Thm. 118/Cor. 101).

## Reproducibility and presentation

1. `LICENSE` (MIT) added at the repository root; the Permanence paragraph
   now points to it.
2. `requirements.txt` added (mpmath, sympy, python-flint; PARI/GP
   separately). The reproduction commands no longer use the checked-in
   `.venv`; the appendix now gives a clean-checkout recipe
   (`pip install -r requirements.txt`) and notes that `.venv` is not part
   of the archive (it is gitignored and untracked).
3. All overfull boxes above 10 pt are fixed: the cusp table is typeset in
   two lines (was 106.6 pt), the regulator-anchor display is rewrapped
   (was 78.3 pt), the certification paragraph is set with `\sloppy`
   (was 40.3 pt), and the family/torsion/code tables are compressed
   (was 22.7 pt). The largest remaining overfull is 9.6 pt.
4. Kept as suggested: the 366-digit value is presented as a consistency
   check; the coarse certified intervals are the proof objects.

## Updated cross-reference

The version of record for this revision is the tagged commit `rev3`.
