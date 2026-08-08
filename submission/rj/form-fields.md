# RJ 投稿表单字段（Editorial Manager 复制粘贴用）

投稿入口：https://www.editorialmanager.com/rama （用你 formal.tex 的账号登录）

## Article Type

Original Paper

## Title

Boyd's conductor-11 Mahler measure conjecture: proof of the split-integral identity (C3), with an exact structural analysis of the family S_k

## Abstract（1610 字符，EM 摘要框足够放）

Boyd's 1998 tables of conjectural identities $m(P)=r\,|L'(E,0)|$ between Mahler measures of two-variable polynomials and $L$-values of elliptic curves begin with the smallest possible conductor, $N=11$. The third conductor-$11$ identity, (C3), concerns the polynomial $S_0=y^2+(x^2+1)y+x^3$, which vanishes on the unit torus, and asserts that a signed split integral $I_{\mathrm{split}}$ of $\log|y|$ around the branch cut equals $b_{11}=L'(E_{11},0)$. We prove (C3). The proof identifies the split integral with a regulator integral along Samart's signed open chain $\tilde\gamma$, closed by a small-branch arc $\beta_0$ into a closed anti-invariant cycle $C'$, and proves its homology class is $2\gamma^-$ with $\gamma^-$ generating $H_1(E,\mathbb{Z})^-$: the period ratio is a-priori integral, and ball arithmetic (Arb) pins it to $2$. A direct regulator computation via Brunault's proved Siegel-unit formula---the symbol $\{x,y\}$ being a pair of modular units on $X_1(11)$, so Bloch's diamond theorem is not needed---then yields $I_{\mathrm{split}}=b_{11}$; the sign is certified in interval arithmetic, and the identity agrees with the numerical value to $366$ digits. For the family $S_k=y^2+(x^2+kx+1)y+x^3$ we determine exactly the torus intersections and the modular-unit/tempered cases; further Boyd-type evaluations are recorded as conjectures. At $k=-1$ (conductor $53$) the mechanism provably fails; an appendix treats Samart's conductor-$17$ analogue conditionally. All computations are reproducible from the accompanying code; every certification step is carried out within interval arithmetic.

## Keywords

Mahler measure; elliptic curve; L-function; Beilinson regulator; modular units; interval arithmetic; certified computation

## MSC 2020

- Primary: 11R06, 11G40
- Secondary: 11G05, 11G55, 65G30

## Cover Letter

用 `cover_letter_rj.md` 的内容，粘贴时把 Markdown 标记（**、##、-）去掉变纯文本即可；日期和称呼自行确认。arXiv 占位符维持 "submitted to arXiv, identifier to be announced"（目前 on hold，这是准确表述）。

## 系统单独询问的声明（与文内 Statements and Declarations 口径一致）

- Funding: The author did not receive support from any organization for the submitted work.
- Competing Interests: The author has no relevant financial or non-financial interests to disclose.
- Data Availability: All certification scripts, frozen certificates, and run logs are available at https://github.com/huiminZheng-collab/boyd-conductor11 , archived at https://doi.org/10.5281/zenodo.21820650 .
- Use of AI tools: documented in the Declaration on the use of AI tools on the title page of the manuscript.

## Suggested Reviewers（可选字段；机构/邮箱请投稿前从公开主页核对后填写）

- Matilde N. Lalín（Université de Montréal）
- François Brunault（ENS de Lyon）
- Mathew Rogers
- Wadim Zudilin（Radboud University）
- Detchat Samart
- Marie-José Bertin

（EM 一般允许跳过或只填 2–4 位；建议至少填 Lalín 和 Brunault。）

## 上传文件

- Manuscript（主文件）：`paper_rj.pdf`（34 页）
- LaTeX Source（如要求/可选）：`paper_rj.tex`（自包含单文件，无外部依赖）
- Cover Letter 文件（如系统要求传文件而非粘贴）：把 cover_letter_rj.md 另存为 .txt 上传

## 提交前最后一查

- [ ] PDF 首页：作者 Huimin Zheng + 完整单位 + Declaration on the use of AI tools 段在首页
- [ ] 参考文献前有 Statements and Declarations 小节
- [ ] 无一稿多投：formal.tex（arXiv:2608.02255）是不同文章，编辑如问起如实说明即可
