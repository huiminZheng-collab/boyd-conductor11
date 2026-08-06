# Fifth referee report on “Boyd’s conductor-11 Mahler measure conjecture”

## Recommendation

**Accept subject to minor revision.**

I reviewed the clean revision tagged `rev5` (commit `4978565`), the response to the technical note, the three new exact scripts and archives, the revised proof, and the rendered manuscript. I also independently reran the exact Sturm computation and the integral Manin-symbol computation in a separate temporary copy. Both reproduced the archived results.

The central objection is now resolved. The manuscript no longer relies on the disputed normalization of the Bloch--diamond formula. Instead, Brunault's proved Siegel-unit regulator theorem is applied directly, and the three exactness gaps in `rev4` have been repaired:

1. `F_total` is first placed unconditionally in `M_2(Gamma_1(121))`; exact vanishing through `q^2420` then exceeds the sharp Sturm bound 1210 and proves `F_total=-2f_11`.
2. The modular symbol `delta={0,3/11}-{0,8/11}` is proved primitive by an integral Manin-symbol presentation and Smith normal form: its coordinates are `(1,-1)` in an integral basis and it generates the rank-one anti-invariant lattice.
3. The constants in the Siegel presentations are now determined exactly as `C_x=-1`, `C_y=+1`, while Brunault's Lemma 5 gives `D_U=2pi`, `D_V=0`; hence the constant correction vanishes exactly.

These changes convert the former high-precision evidence into a coherent exact regulator anchor. The discovery and correction of the earlier false assertion `D_U=0` also increases confidence in the audit trail: the revised proof uses the correct nonzero winding and explains why it causes no change in the regulator.

One small mathematical presentation issue should be corrected before publication; the other comments below are editorial or defensive clarifications.

## Required minor correction: make the cusp-label conversion direct and non-circular

The manuscript currently derives the conversion between Brunault's cusp `P_v` and `k/11` by saying that diamond equivariance leaves the alternatives `v=k^{±1}` and that the Kubert--Lang divisor-order match selects the inverse convention. The same divisor match is then used to establish the Siegel presentations. This makes the exposition appear circular: the candidate presentations help select the cusp map, and the selected cusp map helps prove the presentations.

There is a shorter direct proof from Brunault's definition, requiring no divisor comparison. In the thesis he defines

\[
P_v=\langle v\rangle\infty=[0,v].
\]

Represent the cusp `k/11` by a matrix

\[
\gamma=\begin{pmatrix}k&b\\11&d\end{pmatrix}\in SL_2(\mathbf Z).
\]

The determinant equation gives `kd-11b=1`, hence `d=k^{-1} (mod 11)`. The bottom-row label of the cusp is therefore `[0,d]=[0,k^{-1}]`, so

\[
k/11=P_{k^{-1}}
\]

up to the existing `±1` convention. This proves the inverse-label dictionary directly. The Kubert--Lang order comparison may then remain as an independent check rather than participating in the proof of the dictionary.

The helper `cusp_of_point(a,c)` in `siegel_anchor_step13.py` should be corrected for the same reason. It currently substitutes the bottom row `(c,a)`. For a matrix sending infinity to `a/c`, the second bottom-row entry is a Bézout coefficient `d` satisfying `ad-bc=1`, not generally `a`. This does not affect the reported primitivity result: the seven explicit matrices, the exact check `B*delta=0`, the Smith computation, and the PARI cross-check are independent of that helper, while for the two tested cusps the erroneous representatives happen to land in the same cusp orbit. Nevertheless, the helper and its comment should be fixed.

## Other minor corrections

1. After obtaining `D_U=2pi`, the sentence “the regulator integral does not even depend on these values” is inaccurate. It does depend on the moduli of the constants, particularly on `|C_y|` because `D_U` is nonzero. Replace it by: “Because the constants have modulus one, their correction to the regulator integral vanishes.”
2. In `siegel_anchor_step14.py`, assert `D_U == 1` in units of `2pi`, not merely that its denominator is one. The exact output is already `1`; the stronger assertion makes the certificate match the theorem stated in the paper.
3. “Beyond the conservative bound 2420” should read “through (or at) the conservative bound 2420.” The coefficients are checked for `0 <= n <= 2420`, which is sufficient.
4. The discussion asserting that the transmitted subgroup formulation of Bloch's theorem is “imprecise” is not needed for the proof and remains stronger than what this paper establishes. I recommend describing it neutrally as a normalization/lattice discrepancy, or moving it to a short question for the literature. The direct Siegel-unit proof is now strong enough that the manuscript gains nothing by claiming to correct a general theorem whose original normalization has not been fully reconstructed here.
5. The relevant PDF pages render cleanly. The LaTeX log still contains several small overfull boxes (the largest about 9.6 pt) and PDF-string warnings; these are not mathematical issues but can be cleaned in the final typesetting pass.

## Final assessment

The coefficient-bearing heart of the conductor-11 argument is now supported by an appropriate published theorem and exact integral, modular-symbol, cyclotomic, and Sturm computations. The remaining cusp-label correction has a one-paragraph exact solution and does not require changing any calculated value. Subject to that correction and the small wording changes above, I would recommend acceptance of the main conductor-11 result.

The conductor-17 appendix should continue to be labelled conditional, as it already is; this does not affect the main theorem.
