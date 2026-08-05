# Second referee report on “Boyd's conductor-11 Mahler measure conjecture”

## Recommendation after revision

**Substantially improved, but a further major revision is still required.**

I reviewed commit `0946f73` (`rev2: major revision after GPT referee report`) against the first report. The authors have responded seriously and have corrected most of the problems concerning evidence levels. In particular, the revised paper now distinguishes exact arguments, interval-certified arguments and numerical observations; replaces the floating-point sign choice by a proposed interval certificate; weakens the conductor-53 conclusion to failure of the present mechanism; separates exact torus geometry from conjectural family evaluations; defines the branch used in the split integral; proves primitivity of the anti-invariant generator; gives the cusp dictionary for the conductor-11 curve; and moves the conductor-17 extension to an explicitly conditional appendix.

These are major improvements. The conductor-11 core now looks close to a referee-verifiable proof. Nevertheless, the revision contains two new exact statements that are incorrect as written, and the new sign certificate is not yet fully rigorous. The Brunault/Bertin anchor also needs one final bibliographic and logical correction. Items 1–4 below should be addressed before acceptance.

## Essential remaining corrections

### 1. The exact torus-intersection criterion contains an algebraic error

Proposition `prop:torus` states that a torus point exists if and only if

\[
(k+2\cos\theta)\sin(\theta/2)=0,
\qquad
|k+2\cos\theta|\,|\cos(\theta/2)|\le 2.
\]

The second condition does not follow from the displayed equation in the proof. Writing

\[
2\cos\psi\,e^{i\theta/2}=-(k+2\cos\theta),
\]

and comparing imaginary and real parts gives

\[
\cos\psi\sin(\theta/2)=0,
\qquad
2\cos\psi\cos(\theta/2)=-(k+2\cos\theta).
\]

After eliminating \(\psi\), the relevant size condition is

\[
|k+2\cos\theta|\le 2|\cos(\theta/2)|,
\]

together with the phase condition (equivalently, in the two resulting cases, \(k+2\cos\theta=0\) or \(\theta=0\)). The manuscript instead multiplies by \(|\cos(\theta/2)|\). The eventual union of parameter ranges may remain \([-4,2]\), because the phase equation restricts the possibilities strongly, but the advertised “if and only if” formula and its proof are false as written. Correct the proposition and recheck every downstream description of folds, boundary points and tangencies from the corrected equations.

### 2. The interval sign certificate leaves Arb through ordinary floating-point arithmetic

The new script `code/sign_certify.py` is intended to repair the former uncertified sign choice. Its overall strategy is appropriate, but its `hull` branch is not a rigorous ball-arithmetic operation:

```python
lo = min(float(l1.lower()), float(l2.lower()))
hi = max(float(l1.upper()), float(l2.upper()))
mid, rad = (lo + hi) / 2, (hi - lo) / 2
rad = rad * 1.000001 + 1e-25
return arb(mid, rad), 'hull'
```

Converting outward-rounded Arb endpoints to binary `float`, then applying an empirical inflation factor, does not prove that the returned ball contains the original endpoints. A conversion may round inward, and no theorem in the paper bounds the required padding. Consequently the statements “strictly covering,” “every certification step is carried out within interval arithmetic,” and Theorem `thm:cert`(6) are not justified by the current code.

Construct the convex hull entirely with Arb operations (or convert exact endpoint strings with directed rounding and prove containment), then assert programmatically that the final hull contains both original balls. The same audit should be applied to `k1_sign_certify.py`. Regenerate the archived outputs after the correction. There is also a duplicate evaluation of `x = exp(i*theta)` in `logabs_ball`, and the output labels `I1` and `I2` as already divided by \(\pi\) although the printed balls are the unnormalised integrals; these are minor but should be cleaned up.

### 3. The conductor-17 modular-unit assertion is not established and appears false on the stated modular model

Proposition `thm:family` now claims that for both \(k=0\) and \(k=1\), all rational torsion points supporting the divisors of \(x,y\) are “exactly the rational cusps under the modular identification.” The revision supplies the required explicit cusp table for \(X_1(11)\), so the \(k=0\) claim is now adequately supported. No corresponding table or modular identification is supplied for \(k=1\).

For the conductor-17 strong Weil curve, the natural genus-one modular curve is \(X_0(17)\), which has only two cusps. The divisors displayed in the appendix involve \(O,A,2A,3A\), all four rational torsion points. Torsion does not imply cuspidal support; Manin–Drinfeld has no converse. Thus the statement that \(x,y\) are modular units does not follow and, on the natural \(X_0(17)\) identification, cannot be true for divisors supported at all four distinct points.

This does not directly damage the conductor-11 theorem, and the conductor-17 regulator argument may require only temperedness/K-theory rather than modular units. The remedy is therefore to remove \(k=1\) from the exact “modular units in the family” proposition unless an explicit modular curve, map and cusp-preimage computation are provided. Rewrite the abstract and the claimed family “locus” accordingly. The conditional conductor-17 appendix should state precisely which hypothesis it uses; it should not inherit an unproved modular-unit assertion.

### 4. Correct the Brunault theorem citation and make the anchor statement exact

The revised symbol dictionary is useful, and the fixed-curve ratio principle is a sensible way to eliminate a convention constant. However, the citations currently conflate two different theorems in Brunault's thesis:

- **Théorème 8** in Section 3.7 gives the conductor-11 elliptic-dilogarithm identities, including equation (3.151).
- **Théorème 118** in Section 3.9 proves the Mahler-measure identity (3.208); its proof contains the regulator formula (3.210), defined by (3.211), and invokes Bertin's Theorem 6 and Corollary 6.1.

Lines 1012–1019 and the concluding certified-computation theorem cite “Thm. 8” as if it were the source of (3.210)–(3.211). That is incorrect. Cite Théorème 118 and its proof for the regulator anchor, and Théorème 8/(3.151) only for the \(D_E(P)\)-to-\(L(E,2)\) coefficient.

The ratio lemma is presently conditional on the existence of one symbol-independent constant \(c\). Its proof largely restates that assumption. To make the application complete, cite the precise Bloch regulator theorem that supplies this common \(c\) for the relevant subgroup of tempered symbols with the chosen \(D_E\), diamond quotient and primitive cycle. Then give the exact model/function dictionary used in Brunault–Bertin. The current textual argument (“the only functions named \(x,y\)”) and a nine-digit numerical exclusion are useful checks but should not be the primary identification proof.

The unexplained conductor-11 factor 2 is no longer automatically fatal: if the common-constant theorem and the anchored symbol dictionary are supplied correctly, it cancels in the ratio. But those two premises must be documented at theorem level.

## Points successfully resolved from the first report

The following responses are satisfactory in principle, subject to the corrections above:

1. The main theorem now clearly declares its computer-assisted components.
2. PSLQ failures are no longer presented as proofs of irrationality or nonexistence.
3. The conductor-53 section now proves only failure of the specified closed-cycle/modular-unit mechanism and leaves other identities open.
4. The family evaluations at \(k=2,3,-4,-5,-6\) are labelled numerical conjectures.
5. The branch \(y_-\), branch continuation and chain orientations are described much more clearly.
6. The primitive anti-invariant generator for \(\Delta<0\) is derived explicitly.
7. For conductor 11, the divisor support is matched to the rational cusps of \(X_1(11)\), and tame-symbol/residue issues are addressed.
8. The homology-class computation is isolated as a certified integer-recognition argument with a wide safety margin.
9. The conductor-17 material is separated from the main proof and its normalization status is disclosed.

## Reproducibility and presentation

1. The repository has a `rev2` tag, but the paper says the code is released under the MIT license while no `LICENSE` or `COPYING` file is present at repository root. Add the actual license file.
2. Do not rely on a checked-in Windows `.venv` as a reproduction environment. In this review environment `.venv/Scripts/python.exe` could not launch because it points to a user-specific base interpreter, while the bundled Python lacked `python-flint`. Supply `requirements.txt`, `pyproject.toml`/lock data, or an equivalent clean-environment recipe and test the commands from a fresh clone.
3. `paper.log` has no undefined references or citation errors, but it reports several large overfull boxes, including one over 100 pt near the cusp table and one over 40 pt in the certification discussion. The first pages render legibly, but these large overflows should be corrected before publication.
4. Keep the 366-digit computation as a consistency check; the much coarser rigorous intervals are the actual proof objects.

## Updated assessment

The revision materially increases my confidence in the central conductor-11 strategy. I no longer view the evidence-label problems or the conductor-53 overclaim as obstacles. If items 1–4 are corrected—especially the all-Arb sign certificate and the exact Brunault/Bertin anchor—I would expect the conductor-11 main theorem to merit a much more favorable recommendation. The family modular-unit proposition and the conductor-17 appendix should not delay the main result: they can be weakened or removed if their extra modular claims cannot be proved cleanly.
