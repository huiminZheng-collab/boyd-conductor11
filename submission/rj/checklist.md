# Checklist — The Ramanujan Journal (Springer)

## 投稿入口

- Editorial Manager for The Ramanujan Journal: https://www.editorialmanager.com/rama
- Springer 初投稿格式自由，普通 article 类可接受（作者前一篇 formal.tex 即如此投稿）。

## 所需文件

- `paper_rj.tex`（= 冻结源 report/paper.tex + 参考文献前新增无编号 "Statements and Declarations" 小节）
- `paper_rj.pdf`（本目录已编译好，34 页，0 error，无 undefined reference/citation）
- `cover_letter_rj.md`（提交前微调）

## 移植记录

- 唯一改动：`\begin{thebibliography}` 之前插入 `\section*{Statements and Declarations}`，含四个 `\paragraph`：Funding / Competing Interests / Data Availability / Use of AI tools（文案按 Springer 标准模板）。
- 其余与冻结源逐字节一致；amsart 11pt 原样编译通过。

## 注意事项

- 作者已有一篇在投 RJ 稿件（formal.tex = arXiv:2608.02255，Boyd conductor-11 形式化相关前篇）。两篇主题相邻但内容不同、无一稿多投冲突；cover letter 无需特别说明，如编辑问起如实答复即可。
- Springer 投稿系统会单独询问 Funding/Competing Interests/Data Availability——与文内 Statements and Declarations 小节口径保持一致。
- arXiv 预印本号出来后回填 cover letter 占位符。
