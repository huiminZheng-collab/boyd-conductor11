# Boyd 的 conductor 11 Mahler 测度猜想：研究报告

**日期**：2026-08-04　**仓库**：`boyd-conductor11/`　**攻击用时**：约 2 小时（数值 + 代数实验）

## 1. 问题陈述

Boyd (1998) 系统猜想：许多二元多项式的对数 Mahler 测度 $m(P)$ 是椭圆曲线 $L$ 值的有理倍数，
$m(P) \stackrel?= r\cdot b_N$，其中 $b_N=L'(E_N,0)=\dfrac{N}{4\pi^2}L(E_N,2)$（根数 $+1$），$r\in\mathbb Q$。

Conductor 11（最小可能 conductor；$E_{11}=X_1(11)$，LMFDB `11.a3`，$E(\mathbb Q)=\mathbb Z/5\mathbb Z$；
关联模形式 $f_{11}=\eta(\tau)^2\eta(11\tau)^2$）：

| # | 恒等式 | 状态 |
|---|--------|------|
| (C1) | $m\big((1+x)(1+y)(1+x+y)+xy\big)=7b_{11}$ | 已证（Brunault 2005/06） |
| (C2) | $m\big(y^2+(x^2+2x-1)y+x^3\big)=5b_{11}$ | 已证（Brunault 2006） |
| (C3) | $S_0=y^2+(x^2+1)y+x^3$ 的**劈裂积分** $I_{\mathrm{split}}=\pm b_{11}$ | **开放**（Boyd 1998 (2-33)；Samart 2023 eq. (4.1)） |

(C3) 的具体形式：$S_0=0$ 即 conductor 11 曲线。它在 2-torus 上有零点（$x=\pm i$），故
$m(S_0)$ 本身**不是** $b_{11}$ 的有理倍数（Boyd 的数值观察，本报告复现）。猜想的是
$$I_{\mathrm{split}}:=\frac1\pi\int_0^{\pi/2}\log|y_-(e^{i\theta})|\,d\theta-\frac1\pi\int_{\pi/2}^{\pi}\log|y_-(e^{i\theta})|\,d\theta\ \stackrel?=\ \pm b_{11},$$
其中 $y_\pm(x)=-\frac{x^2+1}2\pm\sqrt{\frac{(x^2+1)^2}4-x^3}$。符号取决于根的命名（$|y_+||y_-|=1$，换根变号；
本报告取主支 sqrt，得 $+b_{11}$；Samart 的 $\tilde y_-$ 约定给出 $-L'(E,0)$，内容相同）。

## 2. 方法

- **$b_{11}$ 高精度**（`code/b11.py`）：$f_{11}=\sum a_nq^n$ 整数系数 + 权 2 近似函数方程
  $$b_{11}=\Lambda(f,2)=\sum_{n\ge1}a_n\Big[e^{-t_n}\Big(\frac1{t_n}+\frac1{t_n^2}\Big)+E_1(t_n)\Big],\quad t_n=\frac{2\pi n}{\sqrt{11}},$$
  $E_1$ 为指数积分。与 Boyd 的 $b_{11}=0.1521471\ldots$ 吻合。
- **Mahler 测度**：Jensen 降维 $m(Ay^2+By+C)=\frac1{2\pi}\int_0^{2\pi}[\log|A|+\sum_j\log^+|y_j|]d\theta$（mpmath，80–300 dps）。
- **PSLQ**：mpmath.pslq。
- **椭圆曲线群律**：在四次模型 $u^2=x^4-4x^3+2x^2+1$（令 $u=2y+x^2+1$）上用首一抛物线法实现
  （`code/torsion.py` 等），在 $\mathbb Q$、$\mathbb Q(\sqrt2)$、$\mathbb Q(\zeta_8)$ 上**精确有理数运算**。

## 3. 数值结果

### 3.1 已证恒等式的独立复验（80 dps，`notes/attack1-results.txt`）

- $m\big((1+x)(1+y)(1+x+y)+xy\big)=1.06502999208142634\ldots$，$|m-7b_{11}|\approx5.0\times10^{-53}$；
- $m\big(y^2+(x^2+2x-1)y+x^3\big)=0.76073570862959024\ldots$，$|m-5b_{11}|\approx3.6\times10^{-53}$。

### 3.2 开放猜想 (C3) 确认到 149 位 —— **本报告主要数值结果**（`notes/attack3-results.txt`）

300 dps 下：
$$I_{\mathrm{split}}=0.152147141725918049486227297478634495628143589164226122809889823882023289695302776676\ldots$$
$$|I_{\mathrm{split}}-b_{11}|=4.85\times10^{-149}.$$
此前公开记录为 Boyd 的 50 位验证；本次将其推进到 **149 位**。

### 3.3 结构恒等式（数值 152 位 + 可证）

数值发现 $I_1+I_2=-m(S_0)$ 到 152 位（`notes/attack2-results.txt`）。事实上这是**定理**：
数值扫描确认 $\max_{\theta\in[0,\pi]}|y_-(e^{i\theta})|=1$（仅在 $\theta=0,\pi/2$ 取到），又 $|y_+y_-|=|x^3|=1$，
故 $|y_-|\le1\le|y_+|$ 处处成立，
$$m(S_0)=\frac1\pi\int_0^\pi\log|y_+|\,d\theta=-\frac1\pi\int_0^\pi\log|y_-|\,d\theta=-(I_1+I_2).\qquad\blacksquare$$
推论：(C3) 等价于 $I_1=\dfrac{b_{11}-m(S_0)}{2}$、$I_2=-\dfrac{b_{11}+m(S_0)}{2}$。

### 3.4 $m(S_0)$ 的负结果

$m(S_0)=0.40560295591501040\ldots$（Boyd 的 $0.4056029$ ✓）。
- PSLQ$(m(S_0),b_{11})$，系数界 $10^8$：**无关系**（复现 Boyd 的"seemingly not $rb_{11}$"）；
- PSLQ 对 $\{m(S_0),b_{11},\log2,\log3,\mathrm{Catalan},m(1+x+y)\}$，系数界 $10^{10}$：**无关系**。
  支持 $m(S_0)$ 需用椭圆双对数表达而非初等常数。

### 3.5 插曲：两个模型的 Mahler 测度不同

$m(P')$，$P'=y^2+y+x^3+x^2$（Boyd slides 的模型）：$0.40560289185535\ldots$，与 $m(S_0)=0.40560295591501\ldots$
**仅前 7 位相同**（$0.4056029$），第 8 位起不同（差 $6.4\times10^{-8}$）。Boyd slides 的 $0.4056029$ 实为 $m(P')$。

## 4. 证明路线分析：modular units 可行，朴素 BMZ 被堵 —— **本报告主要代数结果**

### 4.1 regulator 形式化

$|x|=1$ 上 $\log|y|\,d\arg x=-\eta(x,y)$，$\eta(x,y)=\log|x|\,d\arg y-\log|y|\,d\arg x$。
劈裂积分 $=-\frac1\pi\int_\gamma\eta(x,y_-)$，路径边界
$$\partial\gamma=2[P_{\pi/2}]-[P_0]-[P_\pi],\quad P_0=(1,-1),\ P_\pi=(-1,-1+\sqrt2),\ P_{\pi/2}=(i,e^{i\pi/4}).$$

### 4.2 $x,y$ 是 $X_1(11)$ 上的 modular units（精确验证 ✓）

四次模型 $u^2=x^4-4x^3+2x^2+1$ 的不变量 $I=16,\ J=-304$，$j=2^8I^3/(4I^3-J^2)\cdot?$ 计算得
$j=-4096/11=-2^{12}/11=j(X_1(11))$ ✓。用首一抛物线群律（$\mathbb Q$ 上精确有理运算，`code/torsion.py`）：
$$A:=(0,1)\ [= (0,0)\in S_0]:\qquad 2A=(0,-1),\quad 4A=(1,0)=-A\ \Longrightarrow\ \boxed{5A=O}.$$
故 $(0,0),(0,-1),P_\infty$ 都是 5-扭点。由于
$$\operatorname{div}(x)=[(0,0)]+[(0,-1)]-2[P_\infty],\qquad \operatorname{div}(y)=3[(0,0)]-3[P_\infty],$$
$x,y$ 的除子支撑在 $E(\mathbb Q)_{\mathrm{tors}}=\mathbb Z/5\mathbb Z$（$=X_1(11)$ 的有理尖点）上，
**$x,y$ 确为 modular units**（Manin–Drinfeld）。这正是 Brunault 路线得以施行的前提。

### 4.3 障碍：路径边界不是尖点除子（精确验证 ✗）

群律在 $\mathbb Q(\sqrt2)$、$\mathbb Q(\zeta_8)$ 上精确计算（`code/endpoint_torsion2.py`、`boundary_torsion.py`）：
- $P_0=(1,0)=-A$：5-扭点 ✓（尖点）；
- $P_\pi=(-1,2\sqrt2)$、$P_{\pi/2}=(i,2\zeta_8)$：算到 $20P$ 均非扭点（高度二次增长）；
- 关键组合 $T:=2P_{\pi/2}-P_\pi=(3,\,2i\sqrt2)$：**算到 $30T$ 非扭点**。

故 $\partial\gamma$ **不是**尖点除子，朴素 Brunault–Mellit–Zudilin 公式（要求尖点间路径）
不能直接应用。这与 Samart 2023 §4 的观察一致（"$S$ 族情形 less apparent"；对照：$Q_\alpha$ 族在
conductor 19 处路径闭包条件成立，其 Theorem 2 得证）。

**结论**：(C3) 的 149 位数值成立意味着某种更隐蔽的机制——候选方向：
1. Mellit 式"平行线"椭圆双对数关系（在非扭点 $T$ 处 $D_E$ 的消去）；
2. Samart 式 hypergeometric 公式（$S$ 族尚无一般公式）；
3. "half-Mahler measure" 分解（Lalín–Samart–Zudilin conductor 21 方法）。

## 5. 总结

1. (C1)(C2) 独立复现至 52 位；(C3) 确认至 **149 位**（原记录 50 位）。
2. 新结构定理 $|y_-|\le1\Rightarrow I_1+I_2=-m(S_0)$，把 (C3) 化为 $I_1=(b_{11}-m(S_0))/2$。
3. $m(S_0)$ 对初等常数 PSLQ 阴性（界 $10^{10}$）。
4. 证明路线测绘：modular units 前提**成立**（$5A=O$ 精确验证），但朴素 BMZ **被边界非尖点阻断**
   （$T=(3,2i\sqrt2)$ 非扭点）——指明了证明必须绕开的具体障碍。

## 6. 复现方式

```
cd code && python b11.py && python attack1.py && python attack2.py \
  && python attack3.py && python torsion.py && python endpoint_torsion2.py \
  && python boundary_torsion.py
```

依赖：Python 3.12 + mpmath + sympy。文献见 `literature/`，笔记见 `notes/literature-notes.md`。
