# Response to the fourth referee report

We thank the referee for forcing the normalization question to its source.
The revision tagged `rev4` addresses all three items. In summary: (1) the
regulator anchor of the main theorem is no longer an appeal to Bloch's
diamond theorem in any normalization — the integral is computed **directly**
by Brunault's proved regulator formula for Siegel units, so the "which
lattice" question is bypassed rather than stipulated; (2) the numerical
verification has been redone with a centred scheme, empirically measured
convergence order, and honestly reported residuals; (3) the minor
corrections are implemented.

## 1. The regulator anchor no longer relies on Bloch's diamond theorem

The referee is right on every point of fact: [Lalín–Ramamonjisoa, Thm. 6]
states Bloch's formula for a generator of the anti-invariant *subgroup*
$H_1(E,\mathbb Z)^-$; our index lemma alone cannot decide which lattice
carries the factor-one formula; and the conductor-17 control (index 1)
cannot distinguish the two formulations. Rather than adjudicating the
correct general statement of Bloch's theorem for $\Delta<0$ curves, the
revised proof **does not use Bloch's theorem at all** — nor Bertin's
Theorem 6, nor the Brunault symbol dictionary.

The point is that the symbol $\{x,y\}$ of the main theorem is a pair of
**modular units**: its divisors
$\operatorname{div}(x)=[A]+[2A]-[O]-[3A]$,
$\operatorname{div}(y)=3[2A]-2[3A]-[O]$ are supported on the rational
5-torsion subgroup, which is exactly the set of rational cusps of
$X_1(11)$ (the paper's "Modular units" subsection, with the exact
birational identification and the cusp table). Regulator integrals of
Siegel units are computed by a proved theorem involving no diamond product
and no choice of homology lattice:

> **Brunault [J. Number Theory 163 (2016) 542–569, Thm. 1].**
> $\int_0^{i\infty}\eta(g_u,g_v)=\pi\,\Lambda^*(e_{a,d}e_{b,-c}+e_{a,-d}e_{b,c},0)$
> for Siegel units $g_u,g_v$, with $e_{a,b}$ explicit weight-1 level-$N^2$
> Eisenstein series; arbitrary modular symbols by linearity (his Rem. 2).

The new Theorem (thm:anchor) in §(subsec:anchor) proves
$\int_{\gamma^-}\eta(x,y)=\pm2\pi b_{11}$ in four steps, all exact or
archived (`code/siegel_anchor_step1–11`,
`notes/attack16-siegel-anchor.txt`):

1. **Cusps and divisors.** Under the modular parametrization
   $X_1(11)\to E$, the rational cusps map as $k/11\mapsto m_kA$,
   $(m_k)=(0,2,1,4,3)$, consistent with Brunault's cusp table
   (3.152)–(3.153); the cusp orders of $x\circ\pi$, $y\circ\pi$ follow.
2. **Siegel presentations.** Exact rational linear algebra on the 60 cusps
   of $X(11)$ (Kubert–Lang orders) gives
   $x\circ\pi=-G_4G_5/G_2^2$, $y\circ\pi=G_1G_5^3/(G_2^3G_3)$,
   $G_a=\prod_b g_{a,b}$; each ratio is a priori constant, and the
   constants are roots of unity ($-1$, $+1$ to 70 digits), hence invisible
   to $\eta$.
3. **The cycle.** $\gamma^-=\{0,\tfrac{3}{11}\}-\{0,\tfrac{8}{11}\}$ is
   closed ($3\equiv-8\bmod11$), anti-invariant, and primitive (period
   $=w_{\mathrm{anti}}$, a priori a nonzero integral multiple thereof),
   hence the class $\pm(a-2b)$ of the anti-invariant lemma; continued
   fractions decompose it into seven Manin symbols.
4. **Evaluation.** Term-wise application of Brunault's formula and exact
   $q$-expansion arithmetic (251 coefficients, far beyond the Sturm bound
   2 for $M_2(\Gamma_0(11))$) prove the total weight-2 form is exactly
   $F_{\mathrm{total}}=-2f_{11}$; hence the $\Lambda$-sum is
   $-2\Lambda(f_{11},0)=-2b_{11}$ (functional equation, root number $+1$),
   and the factor $\pi$ of the theorem gives
   $\int_{\gamma^-}\eta(x,y)=-2\pi b_{11}$ with the chosen orientation.

Two verification remarks. First, the factor $\pi$ in Brunault's Theorem 1
was checked against the article's LaTeX source (text extraction had
silently dropped it), and our implementation of the formula was validated
by reproducing Brunault's own conductor-14 application to 60 digits by
three independent methods (Jensen's formula, direct $\eta$-integration
along the Deninger path, and the $\Lambda$-chain). Second, the result is
independently corroborated by the paper's existing 366-digit value of
$I_{\mathrm{split}}$ and by direct numerical integration along $\gamma^-$
(45 digits, opposite orientation).

The factor-of-two question is now purely a remark clarifying the
literature (Lemma (lem:coinvariant), the paragraph "The status of Bloch's
diamond formula (not used)", and Remark (rem:diamondk0)): all verified
data — the new anchor, Brunault–Bertin's value for $\{x_W,y_W\}$, direct
numerical integration on both generators, and the proved conductor-17
computation (where subgroup and quotient coincide) — are consistent with
the factor-one formula holding for the **coinvariant quotient**
$H_1/H_1^+$, and with the transmitted subgroup formulation of
[Lalín–Ramamonjisoa, Thm. 6] being imprecise for $\Delta<0$ curves, a case
their paper never uses. We state this explicitly as the situation of the
literature; it is not an ingredient of the proof.

Consequently, Theorem (thm:cert) now records that the only external input
to the proof of (C3) is Brunault's Theorem 1 (proved in loc. cit.) together
with the functional-equation evaluation $\Lambda(f_{11},0)=b_{11}$; the
Siegel presentations, the Manin decomposition, and the identity
$F_{\mathrm{total}}=-2f_{11}$ are exact computations archived with the
paper. The conductor-17 appendix keeps its conditional status (its symbol
is not cuspidal on $X_0(17)$, so the Siegel-unit method does not apply
there), now explicitly contrasted with the conductor-11 situation.

## 2. The numerical verification, redone honestly

The referee's reading of the archive was correct: the old script used a
left-endpoint rule (hence the $O(1/N)$ error) while applying Richardson
factors appropriate to even-power error models, and the reported accuracies
were overstated. `code/verify_coinvariant.gp` has been rewritten with a
centred discretization ($\log|\cdot|$ sampled at the subinterval midpoint,
multiplied by the exact argument increment), and
`notes/attack15-coinvariant.txt` now reports, verbatim from the run:

- $\int_a\eta$: converges to $0$ spectrally, reaching the 70-digit working
  floor ($<10^{-75}$) from $N=2000$ on — the pairing vanishes on $H_1^+$;
- $\int_b\eta$ and $\int_{a-2b}\eta$: raw errors decay from
  $-1.4\times10^{-5}$, $-6.5\times10^{-6}$ at $N=500$ to
  $-5.4\times10^{-8}$, $-2.5\times10^{-8}$ at $N=8000$, with empirically
  measured order $p=2.000000$ at every doubling (no error model assumed);
- a single $p=2$ Richardson step on the finest pair $N=4000/8000$ leaves
  residuals $-1.8\times10^{-15}$ against $-\pi b_{11}$ and
  $-9.2\times10^{-16}$ against $2\pi b_{11}$; the identity
  $-2\int_b\eta-\int_{a-2b}\eta$ holds to $4.6\times10^{-15}$ on the
  extrapolants.

The manuscript now presents this computation explicitly as a numerical
consistency check, not as evidence for any normalization (and per item 1
the proof no longer involves the quotient normalization at all).

## 3. Minor corrections

1. The references to Lemma (lem:coinvariant)(3) now point to item (1), the
   $\Delta>0$ case (three occurrences).
2. "Three independent tempered symbols" replaced by "three different
   tempered symbols".
3. The normalization comparison with [Lalín–Ramamonjisoa, Def. 5] is now
   stated as a textual, term-by-term identity of the defining series, with
   the 60-digit agreement explicitly labelled a non-rigorous numerical
   check.
4. Terminology (the referee's closing remark of item 1): the lemma now
   defines the coinvariant quotient explicitly as the quotient by the
   invariant sublattice $H_1(E,\mathbb Z)^+$, and states that $\ker(c_*+1)$
   and the quotient are not used interchangeably, as they differ by the
   material index.

## Status of the revision

`rev4` (paper.pdf: 32 pages, clean compile, no undefined references).
Logical dependencies of (C3), end to end: exact divisor/tame/torsion
algebra; the closed-chain lemma and the exact doubling identity
$\int_{C'}\eta=2\int_{\tilde\gamma}\eta$; the Arb-certified class
$[C']=2\gamma^-$; the Siegel-unit anchor (new Theorem (thm:anchor));
the Arb-certified sign. External theorems used: Rodriguez Villegas'
temperedness criterion, Brunault's Theorem 1 [BrunaultSiegel], the
functional equation of $L(f_{11},s)$, and (only for the cross-check
remarks) Brunault's thesis Cor. 3.5.101 and Bertin's Theorem 6.
