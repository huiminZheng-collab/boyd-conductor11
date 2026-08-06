# Fourth referee report on “Boyd’s conductor-11 Mahler measure conjecture”

## Recommendation

**Major revision.  The proposed index-two mechanism is mathematically plausible, but the coefficient-bearing step of the proof is still not established by the cited theorem.**

I reviewed the clean revision tagged `rev3` (commit `00905f0`), including `response_round3.md`, the new coinvariant lemma, `code/verify_coinvariant.gp`, its archived output, and the original statement of Lalín--Ramamonjisoa, Theorem 6.  The revision is a useful advance: the elementary integral-lattice calculation is correct.  For a real elliptic curve with

\[
 c_*a=a,\qquad c_*b=a-b,
\]

one indeed has

\[
 H_1^+=\mathbf Za,\qquad H_1^-=\mathbf Z(a-2b),\qquad
 H_1/H_1^+=\mathbf Z[b],
\]

and the image of the anti-invariant subgroup in the quotient is
`2 Z[b]`.  Thus an anti-invariant differential pairing satisfies
\(\int_{a-2b}\eta=-2\int_b\eta\).  This supplies a credible *possible location* for the observed factor two.

It does not, however, prove that Bloch’s formula has factor one on the quotient generator `[b]`.  That assertion is precisely the point still requiring proof.

## 1. The cited source states the theorem for the anti-invariant subgroup, not the quotient

Equation `(eq:bloch)` is attributed to Lalín--Ramamonjisoa, Theorem 6, but changes the theorem’s cycle.  The cited paper first observes that the regulator pairing vanishes on the invariant part and then says that it suffices to consider the regulator on \(H_1(E,\mathbf Z)^-\).  Its Theorem 6 takes \(\gamma\) to be a generator of the anti-invariant subgroup \(H_1(E,\mathbf Z)^-\).  It does not state the formula for a generator of \(H_1/H_1^+\).

The fact that the pairing descends to \(H_1/H_1^+\) does not permit one to replace the subgroup generator by the quotient generator while retaining the same right-hand side.  In the present lattice, if \(\bar b=[b]\) and \(\gamma^-=a-2b\), then \(\gamma^-=-2\bar b\) in the quotient.  Consequently a formula

\[
 \int_{\gamma^-}\eta=D_E((f)\diamond(g))
\]

would imply

\[
 \int_{\bar b}\eta=-\tfrac12D_E((f)\diamond(g)),
\]

up to orientation, not the factor-one formula asserted in `(eq:bloch)`.  Conversely, the manuscript’s quotient formula would imply a factor two on the subgroup.  Choosing which statement is correct is not an elementary consequence of the index lemma; it is the missing regulator-normalization theorem.

The conductor-17 example cannot by itself transport this integral normalization to the conductor-11 topology.  At conductor 17 the subgroup-to-quotient map has index one, so that example cannot distinguish the two formulations.  Agreement there therefore does not prove which integral lattice is selected when the index becomes two.

To close the proof, the authors should do one of the following:

1. cite a primary-source regulator theorem explicitly formulated on the free sign-coinvariant quotient, with the integral normalization used here; or
2. derive that formulation carefully from the relevant Deligne-cohomology/regulator construction, including the integral lattice and the comparison with the elliptic dilogarithm.

Until this is supplied, the statements at lines corresponding to `(eq:bloch)`, `(eq:anchor)`, and the conclusion of (C3) remain unsupported.  The numerical observations make the quotient formulation plausible, but they cannot substitute for the missing theorem.

There is also a terminology issue worth correcting.  `H_1/H_1^+` should be defined explicitly as the quotient by the invariant lattice (or as the free sign-coinvariant lattice in the precise convention intended).  “Coinvariants” has several standard involution conventions, and two-torsion/index-two distinctions are material in exactly this argument.

## 2. The new numerical verification is overstated and its extrapolation does not match the observed error

The archived output does not support the accuracy claims in `response_round3.md` or the manuscript.

From `notes/attack15-coinvariant.txt`, the reported extrapolated values are:

* `int_a = -0.001128977...`, not zero at the displayed numerical level;
* the error of `int_b` relative to `-pi*b11` is `6.787797...e-5`, not `1e-6`;
* `-2*int_b - int_{a-2b} = -0.0001357559...`, not agreement to `1e-22`;
* only `int_{a-2b}` agrees with its target to about `1e-22`.

Moreover, the errors of `int_a` and `int_b` approximately halve when `N` doubles, indicating a leading \(O(1/N)\) error.  But the function `rich2` uses the factors `4` and `16`, appropriate to cancelling even powers such as \(N^{-2}\) and \(N^{-4}\), rather than performing the claimed linear-in-`1/N` Richardson extrapolation.  The script, archive, prose, and response should be made mutually consistent.  A correct first-order extrapolation should be implemented and its actual residuals reported without promoting a numerical check to a proof.

This issue is secondary to the theorem mismatch, but it matters because the manuscript currently presents the computation as direct confirmation of the disputed quotient normalization.

## 3. Minor corrections

1. Two passages cite `Lemma (lem:coinvariant)(3)`, although the lemma has only items (1) and (2).  The intended reference appears to be item (1), the \(\Delta>0\) case.
2. The phrase “three independent tempered symbols” should be justified if “independent” means linear independence in a specified \(K_2\)-group; otherwise “three different tempered symbols” is safer.
3. Statements that a 60-digit numerical comparison shows two normalizations “coincide verbatim” should be separated into a textual comparison of definitions and an explicitly non-rigorous numerical check.

## Items resolved since round three

The revision does resolve the previous *internal* presentation in which a single constant was simultaneously called 1 and 2.  It now identifies an exact index-two lattice map and uses one clearly described proposed mechanism.  The revision is cleanly committed and tagged, and the earlier torus-intersection, interval-arithmetic, temperedness, citation, licensing, and typesetting corrections remain in place.

## Updated assessment

The manuscript is closer, and the index calculation is not a cosmetic patch: it is the right kind of structure to investigate.  Nevertheless, the main theorem still depends on replacing the generator appearing in the cited Bloch formula by a generator of a larger quotient lattice.  The cited source does not make that replacement, and the conductor-17 control cannot detect it.  I therefore cannot yet recommend acceptance.

This looks potentially correctable in one focused revision if the authors can supply the missing integral regulator theorem (or a complete derivation) and repair the numerical verification and its claims.  Without that theorem, the central coefficient remains conjectural.
