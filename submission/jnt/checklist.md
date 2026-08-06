# Checklist — Journal of Number Theory (Elsevier)

## 投稿入口

- Editorial Manager for JNT: https://www.editorialmanager.com/jnt
- 需要注册/登录 ORCID 关联账号。

## 所需文件

- `paper_jnt.tex`（elsarticle, preprint 12pt；连同编译所需：单文件，内嵌 thebibliography，无外部图片/样式依赖）
- `paper_jnt.pdf`（本目录已编译好，42 页，0 error，无 undefined reference/citation）
- `cover_letter_jnt.md`（提交前转 PDF 或粘贴到 EM 文本框，按需要微调）
- `highlights.txt`（Elsevier Highlights，5 条，每条 ≤85 字符，EM 中单独上传）

## 移植记录（相对冻结源 report/paper.tex）

- `amsart` → `elsarticle [preprint,12pt]`。
- frontmatter 重写：`\author[ahstu]{Huimin Zheng}` + `\ead` + `\affiliation`（elsarticle v3.x key-value 语法）；`\journal{Journal of Number Theory}` 置于 preamble。
- `\subjclass[2020]{...}` → keyword 环境内 `\MSC[2020] 11R06 \sep 11G40 \sep 11G05 \sep 11G55 \sep 65G30`。
- `\keywords{...}` → `\begin{keyword} ... \sep ... \end{keyword}`。
- `\begin{abstract}` 整段移入 `\begin{frontmatter}`（内容一字未改）。
- amsart 的 `\date{August 5, 2026}` 无对应物，删除（elsarticle 自管收稿日期）。
- 正文、附录、thebibliography、AI 声明段、所有自定义宏与 theorem 定义原样保留；amsthm 与 elsarticle 无冲突（preamble 内照常 `\newtheorem`，编译验证通过）。
- 已知坑：elsarticle 的 abstract 必须在 frontmatter 内，否则编译报错——已按要求处理；未发现其他宏冲突。

## 注意事项

- JNT 初投稿接受单文件 PDF；EM 系统会让你分别上传 manuscript、cover letter、highlights。
- 页数从 34（amsart 11pt）涨到 42（elsarticle preprint 12pt），属正常。
- arXiv 预印本号出来后回填 cover letter 占位符。
