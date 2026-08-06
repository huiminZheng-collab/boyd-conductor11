# Response to the fifth referee report (on `rev5`, commit `4978565`)

We thank the referee for the recommendation "accept subject to minor
revision", and for independently rerunning the exact Sturm computation and the
integral Manin-symbol computation in a separate copy — we are glad both
reproduced the archived results. All requested changes are implemented in
`rev6`.

## Required minor correction: cusp-label conversion made direct

Agreed, and thank you for the one-paragraph argument — it is strictly better
than our diamond-equivariance-plus-divisor-match derivation, which did have
the circular appearance the referee describes (presentations selecting the
cusp map, cusp map proving the presentations). The manuscript now derives the
dictionary directly from Brunault's definition `P_v = ⟨v⟩i∞ = [0,v]`: a matrix
`γ = [[k,b],[11,d]] ∈ SL₂(ℤ)` representing the cusp `k/11` has determinant
`kd − 11b = 1`, hence `d ≡ k^{−1} (mod 11)`, so the bottom-row label of `k/11`
is `[0,d] = [0,k^{−1}]`, i.e. `k/11 = P_{k^{−1}}` up to the `±1` convention
(new eq. (labelconv) in step (1) of the proof of the anchor theorem). No
divisor comparison enters the identification; the Kubert–Lang order match is
kept in the text explicitly as "an independent exact check that does not enter
the identification above."

The helper `cusp_of_point(a,c)` in `code/siegel_anchor_step13.py` is fixed the
same way: it now computes a genuine Bézout coefficient `d = a^{−1} (mod c)`
(via `gcdex`, with `det = ad − bc = 1` asserted) instead of substituting the
bottom row `(c,a)`, and the comment records the fix. As the referee
anticipated, the rerun reproduces the archived result verbatim: the seven
explicit matrices, the exact `B·δ = 0`, the Smith computation, the coordinates
`(1,−1)`, and the PARI cross-check (`v⁻(δ) = 1` exact) are all unchanged — for
the two tested cusps the old representatives indeed landed in the same cusp
orbit. The full suite again reports `OVERALL: PASS - delta is a generator of
H_1(X_1(11),Z)^-`.

## Other minor corrections

1. The sentence "the regulator integral does not even depend on these values"
   is replaced by "Because the constants have modulus one, their correction to
   the regulator integral vanishes" — accurate now that `D_U = 2π ≠ 0`.
2. `code/siegel_anchor_step14.py` now asserts `D_U == 1` in units of `2π`
   (check name: "D_U = 2*pi exactly (winding number 1)"), replacing the weaker
   denominator test. Rerun: PASS.
3. "Beyond the conservative bound 2420" is now "in particular through the
   conservative SL₂-index convention 2420", matching the actual range
   `0 ≤ n ≤ 2420`.
4. The Bloch subgroup-formulation discussion is neutralized: it now reads as a
   "normalization/lattice discrepancy — subgroup versus coinvariant-quotient
   generator, related by the index 2 — which the present data resolve in
   favour of the quotient reading; whether this reflects a misnormalization in
   the transmitted subgroup statement of [Lalín–Ramamonjisoa, Thm. 6] (a
   statement their paper never uses) or a convention we have not
   reconstructed, we leave as an open question for the literature." The claim
   that the transmitted theorem is "imprecise" is withdrawn.
5. Typesetting: the remaining overfull boxes are all < 10 pt; we will do a
   final typesetting pass with the copy-editor's requirements if and when the
   paper is accepted.

## Verification summary for `rev6`

- `paper.tex` compiles cleanly: 34 pages, no errors, no undefined references,
  no overfull ≥ 10 pt.
- `siegel_anchor_step13.py` (with the Bézout fix): all 17 checks PASS,
  results identical to the archived ones.
- `siegel_anchor_step14.py` (with the stronger `D_U == 1` assertion): PASS.

The conductor-17 appendix remains labelled conditional, as before.
