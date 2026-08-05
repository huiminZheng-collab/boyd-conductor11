# Referee report on “Boyd's conductor-11 Mahler measure conjecture”

## Recommendation

**Major revision; in its present form I cannot recommend acceptance.**

The manuscript contains an interesting and potentially publishable idea: Boyd's split path is completed by a compensating small-root chain, producing a closed anti-invariant cycle, and a certified period computation is then used to identify its integral homology class. If every normalization and certification step is supplied in a self-contained, independently checkable form, this may indeed yield a proof of the conductor-11 identity. The present manuscript, however, repeatedly blurs four logically different levels of evidence: exact proof, computer-assisted proof with rigorous bounds, high-precision numerical confirmation, and failure of a PSLQ search. Several statements advertised as propositions or theorems are supported only by one of the last two. There is also a central identification/normalization issue on the regulator side that must be resolved explicitly rather than discussed as an unexplained factor-of-two phenomenon.

I regard items 1–5 below as essential.

## Major comments

### 1. The sign in the main theorem is not proved as written

The regulator argument in lines 791–850 yields only

\[
\int_{\gamma^-}\eta(x,y)=\pm 2\pi b_{11}.
\]

The final step (lines 902–917) says that the sign is “pinned by a single numerical evaluation.” An ordinary floating-point evaluation, even at high precision, is not a proof of a signed equality. As written, the argument proves at most

\[
|I_{\mathrm{split}}|=b_{11}.
\]

The same defect occurs in the conductor-17 theorem (lines 1047–1060). This is likely repairable: either derive the sign analytically from the chosen orientations and the signs of the logarithms on the two arcs, or give an interval-arithmetic evaluation whose enclosure is strictly contained in the positive (respectively negative) half-line. The theorem must not depend on an uncertified decimal.

### 2. The Brunault/Bertin anchor and the factor-of-two discrepancy require a proof-level reconciliation

The manuscript's main regulator anchor is the claim in lines 791–850 that Brunault's equations (3.210)–(3.211), together with Bertin, evaluate the regulator of the Weierstrass coordinate symbol \(\{x_{\mathrm W},y_{\mathrm W}\}\), after which a ratio of diamond products transfers the value to the manuscript's symbol.

This identification is central and is not sufficiently documented. Brunault's thesis proves Theorem 118 by referring to Bertin's regulator computation, but the manuscript must spell out exactly which rational functions Brunault/Bertin denote by \(x,y\), on which model, how they pull back to the chosen minimal model, which closed cycle is used, and why that cycle is the same primitive generator (up to sign) used here. A citation to equation numbers plus a computed divisor table is not enough for such a delicate normalization.

The need for this clarification is made acute by Remark 8: the manuscript reports that a published “factor-1” diamond formula misses the conductor-11 integrations by a factor exactly 2, while it works at conductor 17, and admits that no mechanism-level explanation has been obtained. It may be true that a symbol-independent constant cancels in the ratio, but this requires a precise theorem stated with all hypotheses and conventions. In particular, the authors must establish:

1. both symbols define classes in the same regulator domain (including tame-symbol issues);
2. the same regulator map, elliptic dilogarithm, orientation and primitive homology generator are used for both;
3. the diamond-product convention is identical after the explicit birational maps;
4. Brunault's/Bertin's anchored symbol is exactly the claimed symbol, not the coordinate symbol on Bertin's plane cubic or another pullback; and
5. the proportionality constant is genuinely independent of the symbol in the precise setting being used.

Until these points are proved in the text, the assertion that the unexplained factor 2 “does not affect the main theorem” is not yet referee-verifiable. The cleanest revision would formulate a normalization-independent ratio lemma and prove it, then give a complete dictionary from the published anchored symbol to the manuscript's functions.

### 3. Numerical non-detection is repeatedly promoted to a mathematical impossibility

Several statements must be weakened substantially.

- Lines 116–121 call it a “fact” that \(m(S_0)\) is not a rational multiple of \(b_{11}\), but the only evidence given is a bounded PSLQ search. This establishes no such fact. Write “no relation was found up to the stated bound” or provide a proof.
- Proposition 4 (“Negative adjudication of the conductor-53 case”) claims that no Boyd-type identity can exist. The absence of the particular torus cycle and the failure of the modular-unit mechanism rule out that mechanism, not every possible identity between the numerical quantities. A PSLQ failure likewise proves nothing. The conclusion in lines 518–519 that Samart's remark is “most likely a low-precision numerical false positive” is inappropriate without first identifying Samart's precise proposed formula and testing exactly that formula.
- The assertion that a non-torsion point means \(x,y\) are not modular units requires the displayed divisors and a precise identification of the cusps; trivial rational torsion alone is not a substitute for that argument.

I recommend replacing the conductor-53 “theorem” by a section titled, for example, “Failure of the present closed-cycle/modular-unit mechanism at \(k=-1\),” with the exact topological computation retained as a proposition and all universal nonexistence claims removed.

### 4. Proposition 3 is not a theorem and its stated dichotomy is internally inconsistent

The result titled “Structural dichotomy for \(S_k\), numerically established” (lines 443–467) is presented as a proposition about the integer family, but its evidence consists of a few computed values and PSLQ searches. A proposition cannot be “numerically established” in the sense used here.

There are also conflicts within the statement. It says that for \(-4<k<2\) the Boyd-type identity is proved “precisely” for \(k=0,1\), while the range contains infinitely many nonintegral \(k\) unless the restriction to integer \(k\) is carried through explicitly; even for integer \(k\), failure at \(-3,-2,-1\) is only numerical. Conversely, at \(k=2,3,-4,-5,-6\) the displayed identities are conjectural, not established consequences of the “classical Deninger mechanism.”

Separate the content into:

- an exact proposition describing torus intersections as a function of real \(k\);
- exact torsion/modular-unit statements for individually proved cases;
- a table of numerical observations; and
- clearly labelled conjectures.

The abstract and introduction must use the same evidence labels as the body.

### 5. The computer-assisted portion is not yet self-contained enough for a proof paper

The homology identification \([C']=2\gamma^-\) is the novel bridge in the argument. The manuscript relies on Python/Arb scripts and archived output, but a journal proof must state a finite, auditable certificate or a precise algorithm with proved error bounds. At present, important facts are dispersed between prose, scripts and logs, and some descriptions are heuristic in tone (“at two scales,” “a sign or branch error would move the value by 2,” “re-audited”).

The revision should include a dedicated computer-assisted theorem containing:

- exact algebraic input and the exact model/birational maps;
- definitions and orientations of every arc;
- a proof that the listed arcs form the claimed integral singular cycle;
- complete branch-continuation criteria on every subdivision, not merely endpoint proximity tests;
- the exact interval enclosures used for each period integral and the primitive period;
- a rigorous derivation of the integer a priori statement;
- the final interval showing distance less than \(1/2\) from exactly one integer; and
- software/library versions plus one deterministic command that regenerates the certificate.

The current radius \(4.33\times10^{-3}\) is more than sufficient for integer recognition, so extreme numerical precision is irrelevant here. What matters is a short, transparent certificate. Ideally the certificate output should be included as an appendix or deposited in a permanent archive, rather than referenced only by a local repository path.

## Further mathematical comments

### 6. Temperedness, tame symbols and descent to homology

Lines 550–552 say that cyclotomic face polynomials imply \(\{x,y\}\in K_2(E)\otimes\mathbb Q\), while lines 748–755 say that the residues are “torsion, so the pairing descends to homology.” Please give the exact localization argument. State whether the symbol lies in \(K_2(E)\), in \(K_2(E)\otimes\mathbb Q\), or only after multiplying by an integer, and explain explicitly why root-of-unity tame symbols give zero real residue for \(\eta\). “Tempered” and “taming trivially” are not literally equivalent unless conventions and tensoring are specified.

### 7. Modular-unit assertion on \(X_1(11)\)

Lines 540–550 compute divisors supported on rational 5-torsion points and then identify these with the rational cusps. Supply the explicit isomorphism with \(X_1(11)\) and a table mapping all points in the two divisors to cusps. Manin–Drinfeld says cusp differences are torsion; its converse is false in general. Thus “torsion support” alone does not prove “modular unit.” The special identification for \(X_1(11)\) may make the conclusion true, but it must be demonstrated.

### 8. Branch notation and definition of \(y_-\)

The definition of \(I_{\mathrm{split}}\) uses \(y_-\), whereas the proof switches among “small,” “big” and signed continuations through fold points. Give an explicit initial value and continuation rule for \(y_-\). At the fold the labels exchange, so a reader must be able to reconstruct exactly which algebraic branch is integrated on each interval. Expressions such as \(y_{\mathrm{big}}(c^-)= -P\) also conflate a \(y\)-coordinate with a point of the elliptic curve; use point-valued lifts or coordinate-valued notation consistently.

### 9. Exact versus certified structural identity

Proposition 1 is called an exact identity, but the crucial modulus ordering is proved only by a computer-assisted subdivision. That is acceptable if advertised as a computer-assisted proposition and if the certificate satisfies comment 5. Otherwise an elementary analytic proof of the modulus ordering would be preferable. The phrase “certified numerically” should be replaced by “proved by interval arithmetic,” accompanied by the precise certificate.

### 10. The homology generator and real topology

The claim that \(H_1(E,\mathbb Z)^-=\mathbb Z\gamma^-\) and that the displayed purely imaginary period is primitive should be justified, especially because the manuscript itself attributes the factor-of-two anomaly to the one-component real topology. Give the conjugation matrix on a chosen symplectic basis and derive the primitive anti-invariant generator. This will also clarify whether the unexplained factor 2 is actually a nonprimitive-cycle issue.

### 11. The conductor-17 extension should be separated unless fully checked to the same standard

The conductor-17 theorem repeats the uncertified sign step and introduces a second normalization discussion whose reconstruction is not present in the cited paper. It substantially enlarges the burden of verification without being necessary for the conductor-11 result. I recommend either moving it to a separate paper/appendix labelled conditional on the normalization lemma, or supplying the same complete chain, regulator and interval certificates demanded above.

## Exposition and presentation

1. The manuscript currently reads partly as a research log. References to “waves,” latent implementation bugs, local filenames, PSLQ experiments and very long decimal records should be moved to a reproducibility appendix or repository documentation.
2. State one main theorem near the beginning with all hypotheses, conventions and whether it is computer-assisted.
3. Distinguish rigorously among “proved,” “proved with interval arithmetic,” “numerically verified,” “conjectured” and “no relation found within a search bound.”
4. Remove citation syntax such as `\cite{Brunault}*{(3.151)}` and use standard LaTeX citation/page or equation references.
5. Check notation in the divisor discussion: line 773 switches from the torsion point \(A\) to \(P\).
6. The 366-digit computation is impressive as a check but not mathematically stronger than a moderate-precision certified enclosure. It should not dominate the abstract.
7. The paper should cite a permanent source for code and exact output, with a commit hash and license.

## Suggested minimum revision path

A focused revision could become publishable if it does the following:

1. restricts the principal paper to the conductor-11 identity;
2. proves a normalization-independent regulator-ratio lemma and supplies the exact Brunault/Bertin symbol dictionary;
3. replaces both numerical sign choices by analytic or interval proofs;
4. presents the closed-chain and period computation as a self-contained computer-assisted theorem with a compact certificate;
5. demotes all PSLQ-based nonexistence statements and family “dichotomies” to observations/conjectures; and
6. removes or resolves the factor-of-two discussion before using the relevant normalization elsewhere.

The central construction deserves serious consideration, but the paper must first be rewritten so that every advertised theorem is supported by proof-level evidence of the same stated strength.
