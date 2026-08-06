# Cover Letter — Experimental Mathematics

Dear Editors,

Please find enclosed our submission entitled

**"Boyd's conductor-11 Mahler measure conjecture: proof of the split-integral identity (C3), with an exact structural analysis of the family S_k"**

by Huimin Zheng (College of Information and Network Engineering, Anhui Science and Technology University, Fengyang, Anhui 233100, P. R. China; zhhm@ahstu.edu.cn), which we submit for consideration in *Experimental Mathematics*.

This manuscript is original, has not been published previously, and is not under consideration by any other journal. All authors (a single author) are aware of and approve this submission.

## Summary and significance

This is a story of experimental mathematics carried through to a rigorous proof. Boyd's 1998 numerical experiments produced conjectural identities between Mahler measures of two-variable polynomials and L-values of elliptic curves, beginning with the smallest possible conductor, N = 11. Of the three conductor-11 identities, the first two were proved by Brunault; the third, (C3), concerns a polynomial vanishing on the unit torus and asserts that a signed split integral of log|y| around the branch cut equals b_11 = L'(E_11, 0). This paper proves (C3), completing the conductor-11 block of Boyd's table. Highlights:

- **Computation-driven proof**: the proof strategy — closing Samart's signed open chain by a small-branch compensating arc — was found through numerical experimentation; the crucial homology statement (period ratio equal to exactly 2) is then pinned down by certified interval arithmetic (Arb ball arithmetic with rigorous error bounds), not merely observed numerically.
- **Full reproducibility**: every certification step is carried out within interval arithmetic, and all scripts, frozen certificates, and run logs are openly available at https://github.com/huiminZheng-collab/boyd-conductor11 and archived at https://doi.org/10.5281/zenodo.21820650 .
- **New experimental results**: a complete structural analysis of the family S_k, new Boyd-type numerical identities (recorded as conjectures, to up to 70 digits, PARI/GP cross-checked), and a mechanistic analysis of the conductor-53 case showing precisely why the closed-cycle/modular-unit mechanism provably fails there — a negative result that delimits the method.
- A direct regulator computation via Brunault's proved Siegel-unit formula bypasses Bloch's diamond theorem and its normalization ambiguities.

We believe this combination of conjecture discovery, certified computation, structural classification, and mechanism-level counterexample analysis is exactly the kind of work *Experimental Mathematics* aims to publish.

## Transparency

The research was carried out by the author with the assistance of the AI system Kimi (Moonshot AI); a formal Declaration on the use of AI tools appears on the title page of the article, and the author has verified all mathematical content and takes full responsibility for it. A preprint has been submitted to arXiv (identifier to be announced; available on request).

## Suggested referees

- Matilde N. Lalín (Université de Montréal)
- François Brunault (ENS de Lyon)
- Mathew Rogers
- Wadim Zudilin (Radboud University)
- Detchat Samart
- Marie-José Bertin

The author has no conflicts of interest with the suggested referees.

Thank you for your consideration.

Sincerely,

Huimin Zheng
Anhui Science and Technology University
zhhm@ahstu.edu.cn
