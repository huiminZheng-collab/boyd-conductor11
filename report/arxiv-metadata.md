# arXiv submission metadata（Boyd conductor-11 篇）

按 `arxiv-submit` skill 的规则准备：Title 逐字去换行；abstract 删 \cite/\label/\eqref（本文摘要无 \cite），展开自定义宏（\Z→\mathbb Z，\period→\mathrm{period}），行内数学保留。

## Title

Boyd's conductor-11 Mahler measure conjecture: proof of the split-integral identity (C3), with an exact structural analysis of the family S_k

## Authors

Huimin Zheng (zhhm@ahstu.edu.cn)

## Abstract（复制下面整段；已压缩到 1610 字符，符合 arXiv 1920 上限）

Boyd's 1998 tables of conjectural identities $m(P)=r\,|L'(E,0)|$ between Mahler measures of two-variable polynomials and $L$-values of elliptic curves begin with the smallest possible conductor, $N=11$. The third conductor-$11$ identity, (C3), concerns the polynomial $S_0=y^2+(x^2+1)y+x^3$, which vanishes on the unit torus, and asserts that a signed split integral $I_{\mathrm{split}}$ of $\log|y|$ around the branch cut equals $b_{11}=L'(E_{11},0)$. We prove (C3). The proof identifies the split integral with a regulator integral along Samart's signed open chain $\tilde\gamma$, closed by a small-branch arc $\beta_0$ into a closed anti-invariant cycle $C'$, and proves its homology class is $2\gamma^-$ with $\gamma^-$ generating $H_1(E,\mathbb{Z})^-$: the period ratio is a-priori integral, and ball arithmetic (Arb) pins it to $2$. A direct regulator computation via Brunault's proved Siegel-unit formula---the symbol $\{x,y\}$ being a pair of modular units on $X_1(11)$, so Bloch's diamond theorem is not needed---then yields $I_{\mathrm{split}}=b_{11}$; the sign is certified in interval arithmetic, and the identity agrees with the numerical value to $366$ digits. For the family $S_k=y^2+(x^2+kx+1)y+x^3$ we determine exactly the torus intersections and the modular-unit/tempered cases; further Boyd-type evaluations are recorded as conjectures. At $k=-1$ (conductor $53$) the mechanism provably fails; an appendix treats Samart's conductor-$17$ analogue conditionally. All computations are reproducible from the accompanying code; every certification step is carried out within interval arithmetic.

（注：论文内部 `\begin{abstract}` 保持长版不变，仅网站表单用这段。）

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
