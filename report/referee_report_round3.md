# Third referee report on “Boyd's conductor-11 Mahler measure conjecture”

## Recommendation

**One central blocker remains. The other round-two objections have been satisfactorily addressed.**

I reviewed the uncommitted working-tree revision prepared in response to `referee_report_round2.md`, including `response_round2.md`, the new all-Arb sign scripts and certificates, the revised manuscript, the original pages of Brunault's thesis, and the current LaTeX log.

The correction of the torus-intersection proposition is sound. The conductor-17 modular-unit claim has been withdrawn and replaced by the precise temperedness statement actually used. The sign scripts no longer leave ball arithmetic when forming their hulls, and the new containment assertions remove the specific rigor defect identified in round two. The license, dependency recipe and typesetting fixes are also present. I therefore regard round-two items 1–3 and the reproducibility comments as resolved in substance.

The regulator normalization, however, is still internally contradictory. This is not a cosmetic issue: it controls the coefficient in the main theorem.

## Blocking issue: the manuscript simultaneously proves `c=1` and asserts `c=2`

The revised Ratio Lemma states that, with the definitions of Lalín–Ramamonjisoa,

\[
\int_\gamma\eta(f,g)=D_E\big((f)\diamond(g)\big),
\]

so the common constant is `c=1`. This is quoted correctly from their Theorem 6: their regulator is defined by the unnormalised integral \(\int_\gamma\eta\), and the theorem identifies it with the elliptic dilogarithm of the diamond product for a generator \(\gamma\in H_1(E,\mathbb Z)^-\).

Immediately afterward, the manuscript says that on the conductor-11 curve, “the same conventions” give `c=2` exactly. The manuscript's own displayed values make the contradiction explicit:

\[
D_E((x)\diamond(y))=\pi b_{11},
\qquad
\int_{\gamma^-}\eta(x,y)=2\pi b_{11},
\]

and, for the purported Brunault anchor,

\[
D_E((x_W)\diamond(y_W))=-\pi b_{11},
\qquad
\left|\int_{\gamma^-}\eta(x_W,y_W)\right|=2\pi b_{11}.
\]

These cannot coexist with the cited Bloch theorem under one fixed set of definitions. A change of convention cannot repair the argument as presently written: a consistent convention change may replace `c=1` by another fixed value, but the manuscript first uses the cited theorem to establish `c=1` in its declared definitions and then asserts `c=2` in those same definitions.

The issue directly affects the proof. If the cited factor-one formula is applied to the manuscript's diamond value and its certified class \([C']=2\gamma^-\), then

\[
\int_{\tilde\gamma}\eta
=\frac12\int_{C'}\eta
=\frac12\cdot2\cdot\pi b_{11}
=\pi b_{11},
\]

which gives \(I_{\mathrm{split}}=b_{11}/2\), not \(b_{11}\). The observed identity and high-precision computation indicate that some normalization or identification in this chain is off by exactly two, but the current Ratio Lemma does not locate or resolve it.

### Brunault's page does not hide a doubled cycle

I checked the original page 134 of Brunault's thesis. Equation (3.211) defines

\[
r(\{x,y\})=\frac1{2\pi}\int_\gamma\eta(x,y)
\]

with \(\gamma\) generating \(H_1(E(\mathbb C),\mathbb Z)^-\). The superscript `2` after “engendrant” is a footnote marker. The footnote says that the rank-one group has two generators—namely the two choices \(\pm\gamma\)—and hence the regulator is defined up to sign. It does not say that Brunault integrates over twice a primitive generator. Thus the factor two cannot be dismissed as a nonprimitive-cycle convention on this page.

### What must be done

The authors need to locate the factor at source level, not postulate a curve-dependent constant after citing an exact factor-one theorem. At least one of the following data must differ from what the manuscript currently asserts:

1. the rational functions denoted by \(x,y\) in Brunault/Bertin's regulator calculation;
2. the divisor or diamond product of the anchored symbol;
3. the definition of \(D_E\) used in the cited Bloch formula versus the series used in the manuscript;
4. the integral regulator normalization;
5. the homology generator; or
6. the applicability of the quoted theorem to the root-of-unity tame-symbol classes after passage to \(K_2(E)\otimes\mathbb Q\).

The weakest point remains the symbol dictionary. Brunault's Theorem 118 says that Bertin calculated the Mahler measure in terms of a regulator and then writes \(r(\{x,y\})\); the page itself does not explicitly declare that these are the bare Weierstrass coordinate functions. The fact that Section 3.9 refers back to the curve of Section 3.7 fixes the elliptic curve, but not automatically the two rational functions imported from Bertin's calculation. A textual “only functions named \(x,y\)” argument, an exact divisor computation for a candidate pair, and a numerical integral agreeing with the desired anchor do not prove that the candidate pair is Bertin's pair—especially when that identification produces a direct contradiction with Bloch's theorem.

The revision should obtain the exact functions and birational map from Bertin's Theorem 6/Corollary 6.1 (or another explicit primary-source statement), compute their divisors and diamond product, and then apply one precisely stated regulator theorem without changing normalizations mid-argument. Once the correct anchored symbol is known, the ratio argument may well work, but its common constant must have one value throughout.

## Items now resolved

Subject to rerunning the final clean revision, I consider the following satisfactory:

1. The torus-intersection criterion is now correctly split into the cases \(\theta=0\), \(|k+2|\le2\), and \(k+2\cos\theta=0\).
2. `sign_certify.py` and `k1_sign_certify.py` now form non-separated hulls entirely in Arb and assert containment; the former float-padding objection is removed.
3. The manuscript no longer claims that the conductor-17 torsion points are all cusps or that its coordinate functions are modular units.
4. Numerical family observations and the conductor-53 discussion retain appropriate evidence labels.
5. Brunault's Theorem 8 and Theorem 118 are now distinguished bibliographically.
6. `LICENSE` and `requirements.txt` have been added, and the reproduction instructions no longer rely on a nonportable virtual environment.
7. The large overfull boxes have been reduced; the log now has no undefined citations/references and no overflow above 10 pt.

## Repository state

The response calls this revision the tagged commit `rev3`, but at the time of review the changes are still uncommitted and `git tag --list` contains only `rev2`. This is merely a release-process issue: commit the exact reviewed sources and regenerated certificates, then create/update the `rev3` tag only after the mathematical blocker has been resolved.

## Updated assessment

The manuscript has become much cleaner, and nearly all peripheral claims are now calibrated correctly. I cannot recommend acceptance while the main regulator step simultaneously requires `c=1` and `c=2`; that is the coefficient-bearing heart of the claimed proof. If the authors identify the precise Bertin symbol and reconcile it with Bloch's theorem under one consistent normalization, the rest of the conductor-11 argument appears close to publishable form.
