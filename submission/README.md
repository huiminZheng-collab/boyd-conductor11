# submission/ — 期刊投稿包（Boyd conductor-11 论文，冻结于 rev9）

源文件 `report/paper.tex`（amsart, 34 页, 内嵌 thebibliography）为冻结版本，本目录四个投稿包均由它派生，源文件未做任何改动。

## 四个目标刊

| 目录 | 期刊 | 出版方 | tex 处理 | 编译结果 |
|------|------|--------|----------|----------|
| `jnt/` | Journal of Number Theory | Elsevier | 转 elsarticle `[preprint,12pt]`，frontmatter 重写 | pdflatex ×2，42 页，0 error，无 undefined |
| `rj/` | The Ramanujan Journal | Springer | 复制 + 参考文献前加 Statements and Declarations 小节 | pdflatex ×2，34 页，0 error，无 undefined |
| `ijnt/` | International Journal of Number Theory | World Scientific | 逐字节复制（接收后再换 ws-ijnt 样式） | pdflatex ×2，34 页，0 error，无 undefined |
| `expmath/` | Experimental Mathematics | Taylor & Francis | 逐字节复制 | pdflatex ×2，34 页，0 error，无 undefined |

## 各包内容

- `jnt/`：`paper_jnt.tex` + 已编译 PDF、`cover_letter_jnt.md`、`highlights.txt`（Elsevier Highlights，5 条 ≤85 字符）、`checklist.md`
- `rj/`：`paper_rj.tex` + PDF、`cover_letter_rj.md`、`checklist.md`
- `ijnt/`：`paper_ijnt.tex` + PDF、`cover_letter_ijnt.md`、`checklist.md`
- `expmath/`：`paper_expmath.tex` + PDF、`cover_letter_expmath.md`、`checklist.md`

目录内的 `.aux/.log/.out/.pdf` 为编译产物，可留可删；重新编译只需在对应目录跑两遍 `pdflatex <name>`。

## 投稿优先级建议

1. **JNT 首选**：主题（Mahler measure / 椭圆曲线 L 值恒等式）在 JNT 有最强发表史（Mellit、Brunault 相关工作均在此），prestige 与契合度最佳。
2. **RJ 并列首选**：作者已有在投稿件（formal.tex = arXiv:2608.02255），编辑部对作者工作线熟悉；Springer 初投格式自由，成本最低。
3. IJNT 备选（Lalín–Ramamonjisoa 文的发表地，契合但接收后需换 WS 样式）。
4. Exp Math 备选（实验数学叙事最契合，但数论核心读者群相对小）。

**一次只投一家**；被拒后再按顺序转投，转投时更新 cover letter 日期与 arXiv 号。

## 状态表

| 期刊 | 状态 | 日期 | 备注 |
|------|------|------|------|
| JNT | 未投 | — | 包已就绪 |
| RJ | 未投 | — | 包已就绪 |
| IJNT | 未投 | — | 包已就绪 |
| Exp Math | 未投 | — | 包已就绪 |

## 统一事项

- 四封 cover letter 均含：原创/未一稿多投声明、(C3) 证明三大卖点、AI 辅助透明度声明（Kimi, Moonshot AI；正文首页有正式 Declaration）、代码开源（GitHub + Zenodo DOI 10.5281/zenodo.21820650）、6 位 suggested referees（Lalín, Brunault, Rogers, Zudilin, Samart, Bertin）及无利益冲突声明。
- arXiv 预印本号为占位符（"to be announced"），公布后统一回填。
- 移植细节与各刊注意事项见各包 `checklist.md`。
