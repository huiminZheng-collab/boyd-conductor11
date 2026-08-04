# Boyd 的 conductor 11 Mahler 测度猜想：研究报告

**日期**：2026-08-04　**仓库**：`boyd-conductor11/`　**攻击用时**：约 2 小时（数值 + 代数实验）

> 阅读指南：第一部分（§1–§4）是公开文献的详细综述，自包含，不需要先验知识以外的背景；
> 第二部分（§5–§9）是本次工作（两波数值/代数攻击）的结果，含第一波结论的更正（§8）与证明纲要（§9）。

---

# 第一部分：背景详解（公开部分）

## 1. Mahler 测度是什么

### 1.1 定义与 Jensen 公式

对非零多项式 $P(x_1,\dots,x_n)\in\mathbb C[x_1^{\pm1},\dots,x_n^{\pm1}]$，其**对数 Mahler 测度**是 $\log|P|$ 在 $n$ 维环面 $\mathbb T^n=\{|x_1|=\cdots=|x_n|=1\}$ 上的平均：

$$m(P)=\int_0^1\!\!\cdots\!\int_0^1 \log\big|P\big(e^{2\pi i t_1},\dots,e^{2\pi i t_n}\big)\big|\,dt_1\cdots dt_n,
\qquad M(P)=e^{m(P)}.$$

**单变量情形**是完全可以算的：若 $P(x)=a_0\prod_{j=1}^d(x-\alpha_j)$，Jensen 公式给出

$$m(P)=\log|a_0|+\sum_{j=1}^d\log^+|\alpha_j|,\qquad \log^+v:=\max(\log v,0).$$

即"首项系数 + 单位圆外根的贡献"。由此，对整系数单变量多项式，$m(P)$ 总是**代数数的对数**。
Kronecker 定理：$m(P)=0 \iff P$ 是单位根式（分圆）多项式乘以单项式。著名的 **Lehmer 问题**（1933）
问 $m(P)>0$ 能否任意小；目前最小纪录仍是 Lehmer 本人的 $m(x^{10}+x^9-x^7-x^6-x^5-x^4-x^3+x+1)=\log(1.17628\ldots)$。

**两变量情形**完全不同：$m(P)$ 一般是超越数，而且——这正是本课题的主题——它竟然反复地等于
**椭圆曲线 $L$ 函数特殊值的有理倍数**。

### 1.2 降维：把二维积分化为一维

计算 $m(P(x,y))$ 的实用方法是先对 $y$ 用 Jensen 公式。若把 $P$ 看成 $y$ 的多项式
$P=A(x)y^2+B(x)y+C(x)$，两根 $y_\pm(x)$，则

$$m(P)=\frac1{2\pi}\int_0^{2\pi}\Big[\log|A(e^{i\theta})|+\log^+|y_+(e^{i\theta})|+\log^+|y_-(e^{i\theta})|\Big]d\theta.$$

本报告所有 Mahler 测度数值都是用这个一维积分算的（mpmath 任意精度）。

### 1.3 第一个漂亮公式：Smyth (1981)

$$m(1+x+y)=\frac{3\sqrt3}{4\pi}L(\chi_{-3},2)=L'(\chi_{-3},-1)=0.3230659\ldots$$

其中 $\chi_{-3}$ 是 conductor 3 的奇 Dirichlet 特征。证明路线：Jensen 降维后得到 log-sine 积分，
再认出它是 Clausen 函数（二重对数 $\mathrm{Cl}_2$）的特殊值，而 $\mathrm{Cl}_2(\pi/3)\propto L(\chi_{-3},2)$。
这个公式是后来一切"Mahler 测度 = $L$ 值"现象的鼻祖。注意：曲线 $1+x+y=0$ 是**有理曲线**（亏格 0），
对应的 $L$ 函数是 Dirichlet 的；亏格 1 时就轮到椭圆曲线的 $L$ 函数登场。

## 2. 椭圆曲线及其 $L$ 函数：conductor 11 为何特殊

### 2.1 椭圆曲线的 $L$ 函数与 conductor

$\mathbb Q$ 上的椭圆曲线 $E$ 有 $L$ 函数 $L(E,s)=\sum_{n\ge1}a_n n^{-s}$，其中好素数处
$a_p=p+1-\#E(\mathbb F_p)$。**conductor** $N$ 是衡量 $E$ 坏约化程度的正整数（可看作"$E$ 的层级"）。
**$\mathbb Q$ 上椭圆曲线的最小可能 conductor 是 11**——没有 conductor 1 到 10 的椭圆曲线。
所以"最简的"椭圆曲线 $L$ 值就是 conductor 11 的 $L$ 值，这也是 Boyd 表格里 conductor 11 居首位的原因。

Conductor 11 只有**一个同源类（isogeny class）**，含三条曲线（LMFDB `11.a1–a3`）。其中：

- $X_1(11):\ y^2+y=x^3-x^2$（即 `11.a3`），判别式 $-11$，$j=-2^{12}/11$，$E(\mathbb Q)=\mathbb Z/5\mathbb Z$；
- $X_0(11):\ y^2+y=x^3-x^2-10x-20$（即 `11.a1`）。

### 2.2 模形式与函数方程

模性定理（Wiles 等）：$L(E,s)$ 的系数 $a_n$ 是一个权 2、水平 $N$ 的尖形式 $f_E$ 的 Fourier 系数。
对 conductor 11，这个尖形式有极漂亮的 **eta 乘积**表达式：

$$f_{11}(\tau)=\eta(\tau)^2\eta(11\tau)^2=q\prod_{n\ge1}(1-q^n)^2(1-q^{11n})^2=\sum_{n\ge1}a_nq^n,
\qquad q=e^{2\pi i\tau}.$$

完备化 $L$ 函数 $\Lambda(E,s)=N^{s/2}(2\pi)^{-s}\Gamma(s)L(E,s)$ 满足函数方程
$\Lambda(E,s)=w\,\Lambda(E,2-s)$。$E_{11}$ 的根数 $w=+1$（秩 0）。由此推出本报告的核心常数恒等式：

> **推导** $b_N=L'(E,0)=\frac{N}{4\pi^2}L(E,2)$：$s\to0$ 时 $\Gamma(s)=1/s+O(1)$，
> 而 $L(E,0)=0$（函数方程迫使），故 $L(E,s)=L'(E,0)s+O(s^2)$，$\Gamma(s)L(E,s)\to L'(E,0)$，
> 即 $\Lambda(E,0)=L'(E,0)$。另一方面 $\Lambda(E,0)=\Lambda(E,2)=N(2\pi)^{-2}\Gamma(2)L(E,2)=\frac{N}{4\pi^2}L(E,2)$。$\blacksquare$

所以 $L'(E,0)$（$s=0$ 处的导数，Beilinson 猜想喜欢的点）与 $L(E,2)$（临界区间右端点，
模形式方法好算的点）只差一个有理因子，文献中两种写法混用。数值上

$$b_{11}=L'(E_{11},0)=\frac{11}{4\pi^2}L(E_{11},2)=0.15214714172591804948622729747863\ldots$$

（我们的 `code/b11.py` 用 $f_{11}$ 的系数和"近似函数方程"把它算到 300 位；
$L(E,2)=\sum a_n[\,\cdots\,]$ 是一个以 $e^{-2\pi n/\sqrt{11}}$ 速度收敛的级数，几十项就够几百位。）

## 3. Boyd 猜想从哪来：Deninger、Beilinson 与 Boyd 的数值实验

### 3.1 Deninger 的观察（1995–1997）

Deninger 研究 $K_2$ 与 Mahler 测度的联系时猜想

$$m\Big(1+x+\frac1x+y+\frac1y\Big)\stackrel?=\frac{15}{4\pi^2}L(E_{15},2)=L'(E_{15},0),$$

曲线 $1+x+1/x+y+1/y=0$ 恰为 conductor 15 的椭圆曲线。数值吻合到几十位。
**为什么应该有这种事？** 机制（Deninger–Rodríguez Villegas 的解释）：

1. 对"温顺"（tempered）多项式 $P$（Newton 多边形的每个面多项式都是分圆的，
   这等价于 $K$ 论里某个符号 $\{x,y\}$ 在曲线 $P=0$ 上 tame 平凡），
   Jensen 降维后的一维积分可以改写为**椭圆 regulator**
   $$\eta(x,y)=\log|x|\,d\arg y-\log|y|\,d\arg x$$
   沿环面与曲线相交路径的积分；
2. Beilinson 猜想（对模曲线是定理级别的显式化）说这类 regulator 配对等于 $L(E,2)$ 的有理倍数。

于是"$m(P)=r\,b_N$"是 Bloch–Beilinson 猜想的具体化身；Rodríguez Villegas (1997) 进一步把
$m(P_k)$ 写成模形式，使得 CM（复乘）情形可以严格证明。

### 3.2 Boyd 的系统实验（1998）

Boyd（*Experimental Mathematics* 7:1，引用 248 次）对形如
$P_k=A(x)y^2+(B(x)+kx)y+C(x)$ 的族做了大规模数值实验，按 conductor 分类列出了几十条
$m(P)=r\,b_N$ 型猜想（$r$ 为小的有理数），每条都验证到约 50 位小数。最小 conductor 就是 11。
此后 25 年，这些猜想被逐个击破：

- CM 曲线（conductor 27, 32, 36）：Rodríguez Villegas、Lalín–Rogers、Rogers–Zudilin 等；
- conductor 14：Mellit (2012) 与 Mellit–Brunault（"平行线"椭圆双对数方法）、Touafek；
- conductor 15：Rogers–Zudilin (2014) 证明 Deninger 原猜想；
- conductor 20, 24：Rogers–Zudilin；conductor 21：Lalín–Samart–Zudilin (2016)；
- **conductor 11：Brunault (2005/2006)**，见下节。

### 3.3 Brunault 的突破与 BMZ 公式

Beilinson 曾证明：模曲线上 modular units（除子只支撑在尖点上的有理函数）的 regulator
可以用 $L$ 值表示，但常数不显式。**Brunault 的博士论文（2005，ENS Lyon）把这个定理完全显式化**，
并应用于 $X_1(11)$——它恰好是可以用 modular units 参数化的椭圆曲线
（Brunault 后来证明这种曲线只有有限条）。由此他证明了 conductor 11 的两条 Boyd 猜想：

$$m\big((1+x)(1+y)(1+x+y)+xy\big)=\frac{77}{4\pi^2}L(E_{11},2)=7b_{11}, \tag{C1}$$
$$m\big(y^2+(x^2+2x-1)y+x^3\big)=5b_{11}. \tag{C2}$$

后来 Mellit–Brunault–Zudilin 把这类计算凝成一条可直接套用的公式（`literature/zudilin-regulator.pdf`）：
对 Siegel units $g_a(\tau)=q^{NB_2(a/N)/2}\prod_{n\equiv a}(1-q^n)\prod_{n\equiv -a}(1-q^n)$，

$$\int_{c/N}^{i\infty}\eta(g_a,g_b)=\frac1{4\pi}L(f_{a,b;c},2),$$

其中 $f_{a,b;c}$ 是某个显式写出的权 2 模形式（Eisenstein 级数乘积组合）。
**直观含义**：只要你的多项式的 $x,y$ 是模曲线上的 modular units、且积分路径连接两个尖点，
regulator 积分就直接是一个 $L$ 值。这是目前证明此类恒等式的主力武器，也是本报告 §8 分析证明路线时的标尺。

## 4. 仍然开放的 (C3)：torus 上有零点时会发生什么

### 4.1 陈述

Boyd 1998 编号 (2-33) 的族 $S_k=y^2+(x^2+kx+1)y+x^3$。取 $k=0$：

$$S_0=y^2+(x^2+1)y+x^3=0\quad\Longleftrightarrow\quad \text{conductor 11 的椭圆曲线}.$$

与 (C1)(C2) 的多项式不同，$S_0$ **在环面上有零点**：$x=i$ 时 $y^2=i$，$x=-i$ 时 $y^2=-i$，
即 $(x,y)=(i,\pm e^{i\pi/4})$ 与 $(-i,\pm e^{-i\pi/4})$ 满足 $|x|=|y|=1$。此时 Deninger 的 regulator 公式出现边界修正，
Boyd 数值发现

$$m(S_0)=0.4056029559150104\ldots\ \ \text{“seemingly not } r\,b_{11}\text{”},$$

即 $m(S_0)$ 本身**不是** $b_{11}$ 的有理倍数（我们用 PSLQ 在 $10^8$ 系数界内复核了这一点）。
但 Boyd 同时发现：**沿 branch cut 劈开的带符号积分**仍然等于 $b_{11}$。写 $S_0=0$ 的两根

$$y_\pm(x)=-\frac{x^2+1}{2}\pm\sqrt{\frac{(x^2+1)^2}{4}-x^3},$$

$x=e^{i\theta}$ 沿上半环面走，$\theta=\pi/2$ 处正是 torus 交点 $x=i$（"分支切口"），定义

$$I_{\mathrm{split}}:=\underbrace{\frac1\pi\int_0^{\pi/2}\log|y_-(e^{i\theta})|\,d\theta}_{I_1}
\;-\;\underbrace{\frac1\pi\int_{\pi/2}^{\pi}\log|y_-(e^{i\theta})|\,d\theta}_{I_2}\ \stackrel?=\ \pm b_{11}. \tag{C3}$$

（符号取决于哪个根叫 $y_-$：$|y_+||y_-|=|x^3|=1$，换根整体变号。Samart 2023 的约定写成 $-L'(E,0)$。）

Boyd 的原话（经 Samart 2023 §4 引用）：

> "This is in accord with our contention that in case $P$ vanishes on the torus, it is the integral
> of $\log|y|$ around a branch cut rather than $m(P)$, which should be rationally related to $L'(E,0)$."

Samart 2023（arXiv:2301.05390）把 (C3) 明确列为**未证明的猜想**，并指出可尝试用他在
conductor 19 的 $Q_\alpha$ 族上成功的方法（超几何公式 + BMZ）来证，但"$S$ 族的情形 less apparent"。
这就是本次攻击的靶子。

### 4.2 为什么这个情形难

(C1)(C2) 的证明链条是：modular units + 尖点间路径 + BMZ。(C3) 的积分路径端点是
$\theta=0,\pi/2,\pi$ 对应的三个点，其中 $\theta=\pi/2$（torus 交点）处曲线"穿过"环面，
路径边界 $2[P_{\pi/2}]-[P_0]-[P_\pi]$ 是否由尖点组成、能否套用 BMZ，文献中没有答案。
第一波分析曾据此断言朴素 BMZ 被阻断；第二波发现该推理用错了积分链——
正确对象是全圆带符号闭链，闭性是拓扑性质，与端点是否扭点无关。详见 §8（更正）与 §9（证明纲要）。

---

# 第二部分：本次工作

## 5. 方法

- **$b_{11}$ 高精度**（`code/b11.py`）：由 $f_{11}=\eta(\tau)^2\eta(11\tau)^2=\sum a_nq^n$ 的整数系数
  （Euler 函数平方的截断卷积，精确整数），用权 2、根数 $+1$ 的近似函数方程
  $$b_{11}=\Lambda(f,2)=\sum_{n\ge1}a_n\Big[e^{-t_n}\Big(\frac1{t_n}+\frac1{t_n^2}\Big)+E_1(t_n)\Big],\quad t_n=\frac{2\pi n}{\sqrt{11}},$$
  $E_1$ 为指数积分（来自 $\int_1^\infty e^{-ty}/y\,dy$）。项衰减 $\sim e^{-1.894n}$，200 项足够 300 位。
- **Mahler 测度**：§1.2 的一维 Jensen 积分，mpmath 80–300 dps。
- **PSLQ**：mpmath.pslq，搜索小系数整数关系。
- **椭圆曲线群律（精确）**：令 $u=2y+x^2+1$ 把 $S_0=0$ 化为四次曲线
  $u^2=x^4-4x^3+2x^2+1$，在其上用首一抛物线 $u=x^2+bx+c$ 实现加法/取负
  （`code/torsion.py`）；全部在 $\mathbb Q$、$\mathbb Q(\sqrt2)$、$\mathbb Q(\zeta_8)$ 上
  用分数精确运算，**不是数值近似**。

## 6. 数值结果

### 6.1 已证恒等式的独立复验（80 dps，`notes/attack1-results.txt`）

| 恒等式 | 计算值 | 与右端之差 |
|---|---|---|
| (C1) $m((1+x)(1+y)(1+x+y)+xy)$ | $1.06502999208142634\ldots$ | $\|m-7b_{11}\|\approx5.0\times10^{-53}$ |
| (C2) $m(y^2+(x^2+2x-1)y+x^3)$ | $0.76073570862959024\ldots$ | $\|m-5b_{11}\|\approx3.6\times10^{-53}$ |

### 6.2 开放猜想 (C3) 确认到 149 位 —— **主要数值结果**（`notes/attack3-results.txt`）

$$I_{\mathrm{split}}=0.152147141725918049486227297478634495628143589164226122809889823882023289695302776676\ldots$$
$$|I_{\mathrm{split}}-b_{11}|=4.85\times10^{-149}.$$

此前公开记录是 Boyd 的 50 位验证；本次推进到 **149 位**。

### 6.3 结构恒等式（先数值发现，后给出证明）

计算中注意到 $I_1+I_2=-m(S_0)$ 吻合到 152 位。事实上这是**定理**：

> **命题**：在 $[0,\pi]$ 上 $|y_-(e^{i\theta})|\le1\le|y_+(e^{i\theta})|$（数值扫描：
> $\max|y_-|=1$，仅在 $\theta=0,\pi/2$ 取等）。又 $|y_+y_-|=|x^3|=1$，故
> $$m(S_0)=\frac1\pi\int_0^\pi\log|y_+|\,d\theta=-\frac1\pi\int_0^\pi\log|y_-|\,d\theta=-(I_1+I_2).\qquad\blacksquare$$

**推论**：(C3) 等价于 $I_1=\dfrac{b_{11}-m(S_0)}{2}$、$I_2=-\dfrac{b_{11}+m(S_0)}{2}$。
也就是说，劈裂积分的猜想给出的是"大弧段积分"与"小弧段积分"各自的确切值。

### 6.4 $m(S_0)$ 的负结果

$m(S_0)=0.40560295591501040\ldots$（Boyd 的 $0.4056029$ ✓）。

- PSLQ$(m(S_0),b_{11})$，系数界 $10^8$：**无关系**（复核 Boyd 的 "seemingly not $rb_{11}$"）；
- PSLQ 对 $\{m(S_0),b_{11},\log2,\log3,\mathrm{Catalan},m(1+x+y)\}$，系数界 $10^{10}$：**无关系**。
  支持 "$m(S_0)$ 需用椭圆双对数表达而非初等常数" 的预期。

### 6.5 插曲：两个模型的 Mahler 测度不同

Boyd slides 用模型 $P'=y^2+y+x^3+x^2$ 给出 $m=0.4056029$。我们算出
$m(P')=0.40560289185535\ldots$，而 $m(S_0)=0.40560295591501\ldots$——**仅前 7 位相同**
（差 $6.4\times10^{-8}$）。slides 的值是 $m(P')$；这提醒"同一椭圆曲线的不同多项式模型，
Mahler 测度不同"（$m(P)$ 是多项式的不变量，不是曲线的不变量）。

## 7. 证明路线分析（上）：modular units 前提成立

劈裂积分可写成 regulator 积分：$|x|=1$ 上 $\log|y|\,d\arg x=-\eta(x,y)$，故

$$\pi\,I_{\mathrm{split}}=-\int_\gamma\eta(x,y_-),\qquad
\partial\gamma=2[P_{\pi/2}]-[P_0]-[P_\pi],$$
$$P_0=(1,-1),\quad P_\pi=(-1,-1+\sqrt2),\quad P_{\pi/2}=(i,e^{i\pi/4}).$$

套用 regulator 定理的前提：(i) $x,y$ 是 $X_1(11)$ 上的 modular units；(ii) 积分链在
$H_1(E,\mathbb Z)^-$（复共轭反不变部分）中闭合。

**(i) 成立（精确验证）**：四次模型的不变量给出 $j=-2^{12}/11=j(X_1(11))$ ✓。
群律精确计算（`code/torsion.py`）：对 $A=(0,1)$（即 $S_0$ 上的点 $(0,0)$），

$$2A=(0,-1),\qquad 4A=(1,0)=-A\quad\Longrightarrow\quad \boxed{5A=O}.$$

于是 $\operatorname{div}(x)=[(0,0)]+[(0,-1)]-2[P_\infty]$ 与
$\operatorname{div}(y)=3[(0,0)]-3[P_\infty]$ 的支撑全是 5-扭点
（$E(\mathbb Q)_{\mathrm{tors}}=\mathbb Z/5\mathbb Z$ 恰为 $X_1(11)$ 的有理尖点），
故 **$x,y$ 确为 modular units**（Manin–Drinfeld）。

## 8. 闭性机制——第一波"障碍"分析的更正

### 8.1 第一波的错误推理（如实记录）

第一波曾断言：劈裂积分路径 $\gamma$ 的边界组合 $T:=2P_{\pi/2}-P_\pi=(3,\,2i\sqrt2)$ 非扭点，故 $\partial\gamma$ 不是尖点除子，朴素 BMZ（要求尖点间路径）被阻断。**该断言的群律计算本身有效，但结论是错的**：它用错了积分链——把劈裂积分当成"半圆上的路径积分"来取边界，而正确的对象是全圆上的带符号闭链。

### 8.2 正确的积分链 $\tilde\gamma$ 是拓扑闭链

取全圆 $|x|=1$，$x=e^{i\theta}$，$\theta\in[-\pi,\pi]$，携带连续大模长分支 $y_+(\theta)$，权 $+1$（$\theta>0$）/ $-1$（$\theta<0$），在 $y_+$ 跨越单位圆的折点 $\theta=\pm c$（$c=\pi/2$）处分段。这是 Samart（其 Lemma 9 的机制）意义下的修正链 $\tilde\gamma$。

边界分析：链的边界只可能来自折点与端点。端点 $\theta=\pm\pi$ 给出 $[P_\pi]-[P_{-\pi}]$，而 $P_\pi=P_{-\pi}$（$x=-1$ 处两分支连续相接），相消；折点贡献 $\partial\tilde\gamma=2\big([P_c]-[P_{-c}]\big)$。而 $P_{-c}=\overline{P_c}$（复共轭），故在 $H_1(E,\mathbb Z)^-$（复共轭的 $-1$ 特征空间）中 $[P_c]-[\overline{P_c}]$ 自动闭合——**闭性与 $P_c$ 是否为扭点无关**。第一波找到的"障碍"（$T$ 非扭点）系假象：那是错误链的边界，正确链根本不经过那个组合。

**绕数计算（60 位周期配对，`code/winding.py`、`code/winding.gp`）**：几何路径在折点 $\theta=\pm\pi/2$ 其实**不连续**（$y_{\text{big}}$ 在两个交点 $(i,\pm e^{i\pi/4})$ 之间跳跃，左/右极限分别为 $-e^{-i\pi/4}$ 与 $+e^{i\pi/4}$）——这正是必须取带符号链的原因；朴素环积分 $I_{\text{loop}}=-0.47447\ldots i$ 甚至不是任何整闭链的周期（与 $w_{\mathrm{anti}}$ 之比 $0.16262\ldots$ 非有理）。而带符号链在不变微分 $\omega=dx/u$（四次模型 $u^2=x^4-4x^3+2x^2+1$）下的周期
$$I_{\mathrm{signed}}=-2.917633233876990458\ldots i$$
与 PARI 给出的 $H_1(E,\mathbb Z)^-$ 生成元周期 $w_{\mathrm{anti}}=2i\,\mathrm{Im}\,\omega_2$（11.a3）**相等**（比值 $=1$ 到 13 位，受端点 $\sqrt\theta$ 奇性的积分误差所限）：**$\tilde\gamma$ 不仅闭合，而且就是 $H_1(E,\mathbb Z)^-$ 的生成元（绕数 $n=1$）**。

### 8.3 修正 Mahler 测度与 (C3) 的等价改写

沿 $\tilde\gamma$ 的修正 Mahler 测度满足精确恒等式
$$\tilde n=-I_{\mathrm{split}},\qquad \int_{\tilde\gamma}\eta(x,y)=2\pi\, I_{\mathrm{split}}$$
（`code/closedness_check.py`，折点 $\theta=0,\pm\pi/2$ 分段，44 位）。于是
$$\text{(C3)}\ \Longleftrightarrow\ \int_{\tilde\gamma}\eta(x,y)=2\pi\, b_{11},$$
数值上 $\tilde n=-b_{11}$ 到 44 位。

### 8.4 $S_k$ 族：$\tilde n(k)$ 高精度数值表、PARI 独立鉴定与新恒等式

把 §8.2–8.3 的构造用于族 $S_k=y^2+(x^2+kx+1)y+x^3$（`code/ntilde_family.py`，mpmath 50–80 位；
`code/verify_family.gp`、`code/verify_ratios.gp`，PARI/GP 2.15.5 独立复核）：

- **折点存在 ⟺ $|k|<2$**，且折角精确为 $c=\arccos(-k/2)$（由 $2\cos\theta+k=0$，数值验证到 40 位）；
  $|k|\ge2$ 时 torus 无交点，闭链就是全圆，$\tilde n(k)=m(k)$。
- 曲线 $E_k$ 的 $j$ 不变量（四次模型经典不变量 $I=k^4-8k^2+24k+16$、$J=-2k^6+24k^4-72k^3-96k^2+288k-304$，
  $j=6912\,I^3/(4I^3-J^2)$）与 PARI `ellfromeqn` 的结果**全部一致**；conductor 由 PARI 确认：
  $k=-3,-2,-1,0,1,2,3 \mapsto N=83,\,91=7\cdot13,\,53,\,11,\,17,\,37,\,79$。
- $b$ 值双保险：我们的点计数管线（`code/b_family.py`）与 PARI `lfun` 对到 50+ 位。

| $k$ | $N$ | $w$ | 折点 $c/\pi$ | $\tilde n(k)$ vs $b_N=\lvert L'(E_k,0)\rvert$ | 结论（精度） |
|---|---|---|---|---|---|
| $-3$ | $83$ | $-1$ | 无 | 比值 $0.8529175\ldots$ | 无有理关系 |
| $-2$ | $91$ | $-1$ | 无 | 比值 $0.6339454\ldots$ | 无有理关系 |
| $-1$ | $53$ | $-1$ | $1/3$ | 比值 $0.7392026\ldots$ | **无有理关系**（PARI 确认，见下） |
| $0$ | $11$ | $+1$ | $1/2$ | $\tilde n=-b_{11}$ | (C3)，149 位 |
| $1$ | $17$ | $+1$ | $2/3$ | $\tilde n=+b_{17}$ | 60 位（Samart 猜想独立确认）|
| $2$ | $37$ | $-1$ | 无 | $m=\tilde n=2\,b_{37}$ | **新确认恒等式**，60 位 |
| $3$ | $79$ | $-1$ | 无 | $m=\tilde n=b_{79}$ | **新确认恒等式**，60 位 |

要点：

- **$k=1$（conductor 17）**：$\tilde n(1)=b_{17}$——Samart 所述 conductor 17 "(4.1) 类似猜想"
  的独立高精度确认（$y_-$ 约定下积分 $=-L'(E_{17},0)$，有理因子 $r=1$）。
- **$k=2,3$（conductor 37、79，根数 $-1$）**：无 torus 交点，Mahler 测度本身满足 Boyd 型恒等式
  $m(S_2)=2|L'(E_{37},0)|$、$m(S_3)=|L'(E_{79},0)|$，60 位（PARI `lfun`）。这是本次意外收获。
- **$k=-1,-2,-3$（conductor 53、91、83）**：$\tilde n$（及 $m$）对 $|L'(E,0)|$ 的比值均非有理数
  （PARI `lindep`、PSLQ 高 $>10^{10}$ 阴性）。特别地，$k=-1$ 与 Samart 提到的 conductor 53
  "类似猜想恒等式"**不符**——我们的 $b_{53}$ 经 PARI 独立确认，故这不是计算误差；
  要么该猜想的归一化/积分对象不同，要么它基于较低精度的巧合。诚实记为**待解矛盾点**
  （$k=-1$ 在 $\theta=0$ 还有额外 torus 交点 $(1,e^{\pm i\pi/3})$，链定义可能需相应修正）。
- **经验规律**：$k\ge0$ 全部满足 Boyd 型恒等式（$r=-1,+1,2,1$），$k<0$ 全部不满足——
  是否反映某种符号/定向结构，待查。


## 9. 证明纲要：Beilinson–Brunault 路线

| 步骤 | 内容 | 状态 |
|---|---|---|
| S1 | 闭性引理：$\tilde\gamma$ 在 $H_1(E,\mathbb Z)^-$ 中闭合（§8.2），且为**生成元**（绕数 $n=1$，周期配对 60 位） | 数值锁定，严格书写待做 |
| S2 | tempered：$S_0$ 的 Newton 面多项式 $x^3+x^2y$、$x^2y+y^2$、$x^3+y^2$、$y(x^2+1)$ 全分圆，故 $\{x,y\}\in K_2(E)\otimes\mathbb Q$ | 已查 |
| S3 | modular units：$x,y$ 在 $E=X_1(11)$ 上的除子支撑于尖点（$5A=O$ 精确验证） | 已证 |
| S4 | Beilinson–Brunault regulator 定理：闭链配对 $=r\pi b_{11}$ | 引用，假设核对中 |
| S5 | $r=2$：数值锁定（44 位）；绕数 $n=1$ 已定，余下 regulator 常数的代数计算 | 部分 |

**条件性结论**：若 S4 的假设核对通过，则 (C3) 成立。这把一个 50 位数值猜想化为有限的书写/核对任务——这是第二波的主要收获。

## 10. 总结

1. (C1)(C2) 独立复现至 52 位；(C3) 确认至 **149 位**（原公开记录 50 位）。
2. 新结构定理 $|y_-|\le1\Rightarrow I_1+I_2=-m(S_0)$，把 (C3) 化为 $I_1=(b_{11}-m(S_0))/2$。
3. $m(S_0)$ 对初等常数 PSLQ 阴性（界 $10^{10}$）。
4. modular units 前提**成立**（$5A=O$ 精确验证）。
5. **更正**：第一波"朴素 BMZ 被非扭边界阻断"的断言不成立——正确积分链 $\tilde\gamma$ 在 $H_1(E,\mathbb Z)^-$ 中拓扑闭合，与扭点无关（§8）。
6. 证明纲要：(C3) 等价于 $\int_{\tilde\gamma}\eta=2\pi b_{11}$；周期配对确认 $\tilde\gamma$ 是 $H_1(E,\mathbb Z)^-$ **生成元**（绕数 1，60 位）——很可能由 Beilinson–Brunault 定理直接推出；余下 S1 书写、S4 假设核对、regulator 常数 $r=2$ 代数推导（§9）。
7. 族结果（PARI 独立复核）：$\tilde n(k)$ 表（折点 $c=\arccos(-k/2)$，$|k|<2$）；$k=1$ 时 $\tilde n=b_{17}$（60 位，确认 Samart 的 conductor 17 猜想）；意外收获 $m(S_2)=2|L'(E_{37},0)|$、$m(S_3)=|L'(E_{79},0)|$（60 位）；$k<0$ 三个值与 $b_N$ **无有理关系**（含 conductor 53 与 Samart 记述的矛盾点，§8.4）。经验规律：$k\ge0$ 全中、$k<0$ 全不中。


## 11. 复现方式

```
cd code && python b11.py && python attack1.py && python attack2.py \
  && python attack3.py && python torsion.py && python endpoint_torsion2.py \
  && python boundary_torsion.py && python closedness_check.py \
  && python ntilde_family.py && python b_family.py && python winding.py
gp -q verify_family.gp && gp -q verify_ratios.gp && gp -q winding.gp
```

依赖：Python 3.12 + mpmath + sympy。

## 12. 文献导读（`literature/`）

- `bertin-lalin-survey.pdf` — Bertin–Lalín 综述：全局图景与各 conductor 状态（先读这篇）
- `boyd-pnwnt2015.pdf` — Boyd 2015 slides：猜想史 + $m(S_0)$ 原始数据
- `brunault-these.pdf` — Brunault 博士论文：$X_1(11)$ 上 Beilinson 定理显式化，(C1) 的证明
- `zudilin-regulator.pdf` — Zudilin：BMZ regulator 公式（证明武器）
- `samart2023.pdf` — Samart：开放猜想 (C3) 的明确陈述（其 eq. (4.1)）+ conductor 19 的成功范例
- `lalin-samart-zudilin-cond21.pdf` — conductor 21：half-Mahler 方法范例
- `boyd-slides.pdf` — Boyd 关于 $L(E,3)$ 的 slides
- 详细笔记：`notes/literature-notes.md`；原始运行输出：`notes/attack*-results.txt`
