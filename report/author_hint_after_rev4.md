# Technical note to the authors after `rev4`

The new Siegel-unit route is a substantial improvement. It genuinely bypasses the unresolved subgroup-versus-quotient normalization in the Bloch--diamond formula, and Brunault's Theorem 1 is the right kind of primary theorem for an absolute regulator calculation. I believe this route is potentially sufficient to prove the conductor-11 identity.

At present, however, three places advertised as exact still rely on finite-precision or finite-truncation evidence. These are focused, repairable issues.

## 1. Prove the modular-form space membership before applying the Sturm bound

The most important remaining gap is the statement

\[
F_{\mathrm{total}}=-2f_{11}\in M_2(\Gamma_0(11)).
\]

`code/siegel_anchor_step9.py` computes 251 exact Fourier coefficients and tests that they agree with a linear combination of `f11` and `E2,11`. It does **not** prove that `F_total` belongs to `M_2(Gamma_0(11))`: the line printed as an `M_2(Gamma_0(11)) test` is only a finite coefficient comparison. Sturm's bound 2 becomes available only after membership in that space has been proved independently.

The ambient fact supplied by Brunault is that each `e_{a,b}` has level `N^2`; for `N=11` this gives an a-priori much larger group. Agreement through `q^250` does not by itself lower the level. Thus the current order of inference is circular:

1. assume/announce that `F_total` lies in `M_2(Gamma_0(11))`;
2. use the dimension and Sturm bound of that space;
3. infer the identity from the first 251 coefficients.

Two clean repairs are available.

* **Preferred:** use the exact transformation law for the Eisenstein series `e_{a,b}` to prove invariance of the finite symbolic expression defining `F_total` under a set of generators/coset relations for `Gamma_0(11)`. Holomorphy then follows from the Eisenstein-series construction. Once this descent is written down, Sturm bound 2 proves the identity from the first few coefficients.
* **Fallback:** keep `F_total` in a rigorously established ambient space such as the relevant weight-2 level-121 space and apply that space's actual Sturm bound. A safe unreduced `Gamma_1(121)` bound is on the order of 2420 coefficients (the precise value depends on the `±I` convention), not 2. This is computationally less elegant but logically straightforward.

The same repair must justify the stronger per-symbol claims that two summands equal `-f11` and the other five vanish identically; presently these are also finite coefficient observations.

## 2. Certify the modular symbol as a primitive integral generator

The cycle

\[
\delta=\{0,3/11\}-\{0,8/11\}
\]

is plausibly the desired anti-invariant generator, but the proof currently says that its period is an a-priori nonzero integer multiple of `w_anti` and agrees with `w_anti` to 60 digits. Ordinary 60-digit numerical agreement does not identify that integer rigorously.

The shortest exact repair is an integral modular-symbol computation: impose the Manin relations for `Gamma_1(11)`, compute the cuspidal homology lattice by Smith normal form, express the displayed seven-symbol chain in an integral basis, and show that its anti-invariant coordinate has gcd 1. Equivalently, exhibit an integral class whose intersection with `delta` is `±1`.

An interval alternative would also work. Since `delta=n gamma^-` with `n in Z`, compute its period with Arb and enclose `period(delta)/w_anti` in a ball of radius less than `1/2` around 1, exactly as the manuscript already does for `C'`. The present mpmath calculation can remain as a check but should not be the proof of primitivity.

This point matters because the Siegel-unit calculation evaluates the regulator on the specific modular symbol `delta`; the main proof later needs that symbol to be the same primitive generator used in `[C']=2 gamma^-`.

## 3. Replace the numerical constants in the Siegel presentations by exact arguments

Matching cusp divisors proves only

\[
x\circ\pi=C_x\frac{G_4G_5}{G_2^2},\qquad
y\circ\pi=C_y\frac{G_1G_5^3}{G_2^3G_3}
\]

for nonzero constants `C_x,C_y`. The statement that these constants are roots of unity because numerical evaluation gives `-1,+1` to 70 digits is not an exact proof.

There are again two easy ways to close this.

* Compare the leading local `q`-terms at the infinite cusp. The Siegel products have explicit leading powers and cyclotomic leading coefficients, while the modular parametrization and Weierstrass equation determine the leading terms of `x∘pi` and `y∘pi`. This should prove `C_x=-1` and `C_y=1` algebraically.
* Or avoid determining the constants: write explicitly
  \[
  \eta(C_xU,C_yV)=\eta(U,V)+\log|C_x|\,d\arg V-\log|C_y|\,d\arg U
  \]
  and prove by the seven Manin-symbol decomposition and Brunault's Lemma 5 that both argument periods vanish exactly. The archive says `Du=Dv=0`; that exact calculation should appear in the proof rather than the unsupported root-of-unity assertion.

The cusp-to-torsion map should likewise be made exact in the text. Brunault's equations (3.152)--(3.153) give an exact cusp table. Spell out the conversion between his label `P_v` and the cusp `k/11` (the numerical tuple used here corresponds to the inverse-label convention `v=k^{-1} mod 11`, up to `±`). Then derive `(m_1,...,m_5)=(0,2,1,4,3)` from that table. The 60-digit Abel integral should be described only as a check.

## Small reproducibility and wording fixes

* `code/siegel_anchor_step10.py` is explicitly marked in the archive as defective and numerically unreliable, but the manuscript advertises `step1--11` collectively as the proof archive. Remove step 10 from the proof/reproduction sequence, rename it as a discarded experiment, or correct it.
* One occurrence of “three independent tempered symbols” remains in `report/paper.tex`; `response_round4.md` says all such occurrences were changed to “three different tempered symbols.”
* Replace “251 coefficients prove” by “251 coefficients verify” until the space-membership argument above is supplied.

## Bottom line

The former Bloch-normalization blocker has been bypassed in the correct way, and the numerical/exact evidence strongly suggests that the Siegel-unit evaluation is right. The remaining work is not to discover a new formula; it is to upgrade three lattice/space/constant identifications from high-confidence computation to proof. If the authors add those exact arguments, the new regulator anchor would be substantially more convincing and could plausibly close the central objection.
