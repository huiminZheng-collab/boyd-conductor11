# arXiv submission metadata（Boyd conductor-11 篇）

按 `arxiv-submit` skill 的规则准备：Title 逐字去换行；abstract 删 \cite/\label/\eqref（本文摘要无 \cite），展开自定义宏（\Z→\mathbb Z，\period→\mathrm{period}），行内数学保留。

## Title

Boyd's conductor-11 Mahler measure conjecture: proof of the split-integral identity (C3), with an exact structural analysis of the family S_k

## Authors

Huimin Zheng (zhhm@ahstu.edu.cn)

## Abstract（复制下面整段）

Boyd's 1998 tables of conjectural identities $m(P)=r\,|L'(E,0)|$ between Mahler measures of two-variable polynomials and $L$-values of elliptic curves begin with the smallest possible conductor, $N=11$. Two of the three conductor-$11$ identities were proved by Brunault via an explicit version of Beilinson's theorem on modular units. The third one, (C3), concerns the polynomial $S_0=y^2+(x^2+1)y+x^3$, which vanishes on the unit torus, and asserts that a signed split integral $I_{\mathrm{split}}$ of $\log|y|$ around the branch cut equals $b_{11}=L'(E_{11},0)$. We prove (C3). The proof identifies the split integral with a regulator integral along Samart's signed open chain $\tilde\gamma$; we close $\tilde\gamma$ by a small-branch compensating arc $\beta_0$ into a closed anti-invariant integral cycle $C'=\tilde\gamma+\beta_0$, and prove that its homology class is $2\gamma^-$, where $\gamma^-$ generates $H_1(E,\mathbb{Z})^-$: the period ratio $\mathrm{period}(C')/w_{\mathrm{anti}}$ is a-priori an integer, and a ball-arithmetic (Arb) computation with certified error bounds pins it to $2$ (in particular it is non-zero). Combined with the exact integral identity $\int_{\beta_0}\eta=\int_{\tilde\gamma}\eta$ and a direct regulator computation via Brunault's proved Siegel-unit formula---the symbol $\{x,y\}$ being a pair of modular units on $X_1(11)$, so no use of Bloch's diamond theorem is needed---this yields $I_{\mathrm{split}}=b_{11}$; the sign is certified in interval arithmetic, and the identity agrees with the numerical value to $366$ digits. For the family $S_k=y^2+(x^2+kx+1)y+x^3$ we determine exactly the torus intersections (as a function of real $k$), the single case where $x,y$ are modular units ($k=0$), and the cases where the symbol is tempered with torsion-supported divisors; all further Boyd-type evaluations we find are numerical observations, recorded as conjectures ($k=2,3$ to $70$ digits, $k=-4,-5,-6$ to $25$ digits, PARI/GP cross-checked). At $k=-1$ (conductor $53$) the closed-cycle/modular-unit mechanism provably fails; Samart's suggested analogue, for which no formula is on record, is left open. An appendix applies the same method to Samart's conductor-$17$ analogue $\tilde n(1)=b_{17}$, conditional on a normalization lemma in the conventions of Lal\'in--Ramamonjisoa. All computations are reproducible from the accompanying code; every certification step is carried out within interval arithmetic.

## Categories

- Primary: **math.NT**（账号已获背书，直接可投）
- Cross-list：可不选（前两篇都只选 math.NT）；如想加，math.AG

## Comments

34 pages. Certification scripts available at https://github.com/huiminZheng-collab/boyd-conductor11 and archived at https://doi.org/10.5281/zenodo.21820650 . Research carried out with the assistance of the AI system Kimi (Moonshot AI); see the declaration in the article.

## Journal-ref

留空（暂未投期刊）。

## License

arXiv 标准 perpetual non-exclusive（与前两篇一致）。

## 上传文件

- `report/paper.tex` 单文件（内嵌 thebibliography，无 .bbl 依赖；无 \input/\include/\includegraphics 外部依赖——提交前最后一查）
- 不传 PDF（arXiv 自动编译）；amsart/tikz/booktabs/hyperref(hidelinks) 均为 arXiv 标配
