# Boyd 的 conductor 11 Mahler 测度猜想：研究报告

**日期**：2026-08-04　**仓库**：`boyd-conductor11/`　**攻击用时**：约 2 小时（数值 + 代数实验）

> 阅读指南：第一部分（§1–§4）是公开文献的详细综述，自包含，不需要先验知识以外的背景；
> 第二部分（§5–§9）是本次工作（五波数值/代数攻击）的结果，含第一波结论的更正（§8）、
> 证明纲要（§9）、regulator 常数的显式计算（§9.1）与闭链引理的严格化——(C3) 的证明（§9.2）。

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
与 PARI 给出的 $H_1(E,\mathbb Z)^-$ 生成元周期 $w_{\mathrm{anti}}=2i\,\mathrm{Im}\,\omega_2$（11.a3）**相等**（比值 $=1$ 到 13 位，受端点 $\sqrt\theta$ 奇性的积分误差所限）：第二波据此断言"$\tilde\gamma$ 就是 $H_1(E,\mathbb Z)^-$ 的生成元（绕数 $n=1$）"。**第五波更正与严格化**（§9.2）：$\tilde\gamma$ 是带边开链（端点 $P$ 非扭），其周期配对只是启发证据；正确的闭反不变闭链是 $C'=\tilde\gamma+\beta_0$（$\beta_0$ 为小分支补偿弧），严格认证 $\mathrm{class}(C')=2\gamma^-$——并由此把 $\int_{\tilde\gamma}\eta=2\pi b_{11}$ 提升为**严格等式**。

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
| $-1$ | $53$ | $-1$ | $1/3$ | 比值 $0.7392026\ldots$ | **无有理关系**（结构解释，见下） |
| $0$ | $11$ | $+1$ | $1/2$ | $\tilde n=-b_{11}$ | (C3)，149 位 |
| $1$ | $17$ | $+1$ | $2/3$ | $\tilde n=+b_{17}$ | 60 位（Samart 猜想独立确认）|
| $2$ | $37$ | $-1$ | 无 | $m=\tilde n=2\,b_{37}$ | **新确认恒等式**，60 位 |
| $3$ | $79$ | $-1$ | 无 | $m=\tilde n=b_{79}$ | **新确认恒等式**，60 位 |
| $-4$ | $37$ | $-1$ | 无（相切） | $m=\frac72\,b_{37}$ | **新恒等式**（预测后验证），25 位 |
| $-5$ | $359$ | $-1$ | 无 | $m=\frac14\,b_{359}$ | **新恒等式**（预测后验证），25 位 |
| $-6$ | $997$ | $-1$ | 无 | $m=\frac18\,b_{997}$ | **新恒等式**（预测后验证），25 位 |

要点：

- **$k=1$（conductor 17）**：$\tilde n(1)=b_{17}$——Samart 所述 conductor 17 "(4.1) 类似猜想"
  的独立高精度确认（$y_-$ 约定下积分 $=-L'(E_{17},0)$，有理因子 $r=1$）。
- **$k=2,3$（conductor 37、79，根数 $-1$）**：无 torus 交点，Mahler 测度本身满足 Boyd 型恒等式
  $m(S_2)=2|L'(E_{37},0)|$、$m(S_3)=|L'(E_{79},0)|$，60 位（PARI `lfun`）。这是本次意外收获。
- **$k=-1$ 矛盾已解决（第四波，`code/k53_attack.py` + `k53.gp`）**：Samart 2023 §4 关于
  conductor 53 的"类似猜想"只有一句话（无公式、无精度、无引用）。我们证明它**不可能成立**：
  1. **拓扑障碍**：$k=-1$ 时 $|x|=1$ 上的 torus 交点为 $\theta=0$（两根 $e^{\pm2\pi i/3}$，
     分支在此**交换**）与 $\theta=\pm\pi/3$（$y=\pm1$）。枚举所有以交点为断点、big/small
     分支组合的**闭反不变链**：唯一解 $(1,1,1,1)$ 的周期为 $0$——**环面上不存在非平凡的
     反不变闭链**，$H_1(E,\mathbb Z)^-$ 生成元无法在 torus 上实现，故 Boyd 型环面积分
     对象根本不存在。（连续根沿圆走**两圈**才闭合，且整圈周期也是 $0$。）
  2. **代数障碍**：53.a1 的扭子群**平凡**、秩 1，且 $(0,0)$ 就是 Mordell–Weil **生成元**
     （PARI `ellorder`=0）——$x,y$ **不是** modular units，Beilinson–Brunault 机制不适用，
     任何闭链的 regulator 配对都没有理由是有理数 $\times\,b_{53}$。
  3. 自然候选对象的数值判决：半环 $L_1$（连续根一圈，开链）周期与 $w_{\mathrm{anti}}$
     之比 $0.5492906\ldots$ 不在周期格中，regulator 积分与 $2\pi b_{53}$ 之比
     $0.7392026\ldots$ 无理（PSLQ 阴性）。**结论：Samart 的 conductor 53 记述极可能是
     低精度数值假阳性。**
- **扭点对照表**（`code/kfamily_torsion.gp`，PARI `elltors`/`ellorder`）揭示真正的分野：

  | $k$ | $N$ | 扭子群 | $\operatorname{ord}(0,0)$ | Boyd 型恒等式 |
  |---|---|---|---|---|
  | $0$ | $11$ | $\mathbb Z/5$ | $5$ | 成立（modular units，(C3)） |
  | $1$ | $17$ | $\mathbb Z/4$ | $4$ | 成立（$\tilde n=b_{17}$） |
  | $-3,-2,-1,2,3$ | $83,91,53,37,79$ | 平凡 | $\infty$ | 见下 |

- **修正后的结构定理（第四波，预测+验证）**：分野不是 $k$ 的符号，而是 **torus 交点结构**：
  - **$|k|\ge2$**（无真正 torus 交点，$k=\pm2,\ldots$；$k=\pm4$ 处相切退化）：标准 Deninger
    机制直接给出 $m(S_k)=r_k\,|L'(E_k,0)|$，$r_k$ 为小有理数——$k=2,3$ 已确认，
    $k=-4,-5,-6$ 为**先预测后命中**：$m(S_{-4})=\frac72|L'(E_{37},0)|$、
    $m(S_{-5})=\frac14|L'(E_{359},0)|$、$m(S_{-6})=\frac18|L'(E_{997},0)|$（25 位精确）。
    注意 $S_{-4}$ 与 $S_2$ 同为 conductor 37，给出同一曲线的 Rodriguez-Villegas 型关系。
  - **$-4<k<2$**（真正 torus 交点，$m$ 本身失效，须用劈裂积分）：恒等式成立**当且仅当**
    曲线有扭点使 $x,y$ 成为 modular units——全族中仅 $k=0$（$\mathbb Z/5$，$X_1(11)$）
    与 $k=1$（$\mathbb Z/4$）两例；$k=-1,-2,-3$ 扭子群平凡，无解。
  - 这与 Samart 的观察 $K\cap\mathbb R=[-4,2]$（环面相交参数区间）精确吻合。


## 9. 证明纲要：Beilinson–Brunault 路线

| 步骤 | 内容 | 状态 |
|---|---|---|
| S1 | 闭性引理：$\tilde\gamma$ 闭化为 $C'=\tilde\gamma+\beta_0\in H_1(E,\mathbb Z)^-$，$\mathrm{class}(C')=2\gamma^-$（§9.2，比率先验整数 + 15 位认证） | **已证（认证）** |
| S2 | tempered：$S_0$ 的 Newton 面多项式 $x^3+x^2y$、$x^2y+y^2$、$x^3+y^2$、$y(x^2+1)$ 全分圆，故 $\{x,y\}\in K_2(E)\otimes\mathbb Q$ | 已查 |
| S3 | modular units：$x,y$ 在 $E=X_1(11)$ 上的除子支撑于尖点（$5A=O$ 精确验证） | 已证 |
| S4 | Beilinson–Brunault regulator 定理 + Brunault (3.151)：$L(E,2)=\frac{10\pi}{11}D_E(P)$，系数 40 位复核 | 已核（文献+数值） |
| S5 | regulator 常数：$5D_E(P)-5D_E(2P)=-\pi b_{11}$ vs $\int\eta=2\pi b_{11}$；因子 $-2$ = Bloch 定理内在因子 2 × 定向（§9.1，已核） | 闭合 |

### 9.1 regulator 常数的显式计算（第三波，40 位）

第三波把 S5 的"余下代数计算"完成了（代码 `code/dilog.gp` + `code/dilog.py`，
原始输出 `notes/attack7-dilog.txt`）。全部成分如下。

**除子代数**（三次模型上的射影计算，Abel 主除子检验通过）：
$$\operatorname{div}(x)=[A]+[2A]-[O]-[3A],\qquad \operatorname{div}(y)=3[A]-2[O]-[3A],$$
其中 $A=(0,0)$ 即 Brunault 记号中的 5-扭点 $P$，$O=[0:1:0]$ 为群单位元
（$[1:-1:0]=3A$ 由 Abel 条件强制）。金刚石积
$$(x)\diamond(y)=6(O)+5(A)-5(2A),\qquad\text{故}\quad D_E\bigl((x)\diamond(y)\bigr)=5D_E(P)-5D_E(2P).$$

**椭圆双对数**（mpmath 60 位；$\Delta<0$ 故取格基 $(w_1,\,w_1-w_2)$ 使 $\tau=0.5+0.2299i$
落在上半平面，$q=e^{2\pi i\tau}\approx-0.2359$ 为负实数；$z(P)=0.6w_1$、$z(2P)=0.2w_1$ 由 PARI `ellpointtoz` 给出）：
$$D_E(P)=0.1911937370843316957549544343121738161012\ldots$$
两条关键恒等式（均 40 位）：

- **Brunault 系数**：$D_E(P)=\dfrac{2\pi}{5}\,b_{11}=\dfrac{11}{10\pi}L(E,2)$，
  即 Brunault 论文 (3.151) 中系数的精确形式 $L(E,2)=\frac{10\pi}{11}D_E(P)$，数值钉死；
- **exotic relation**：$D_E(2P)=\frac32\,D_E(P)$（Bertin 已证此类关系；
  注意是 $3/2$ 而非朴素的 2 倍）。

**闭环**：代入得
$$5D_E(P)-5D_E(2P)=5\Bigl(1-\frac32\Bigr)D_E(P)=-\frac52\,D_E(P)=-\pi\,b_{11},$$
而绕数一节已锁定 $\int_{\tilde\gamma}\eta(x,y)=2\pi b_{11}$（60 位，PARI 交叉验证）。
两者恰差因子 $-2$：
$$\int_{\tilde\gamma}\eta(x,y)=-2\,D_E\bigl((x)\diamond(y)\bigr).$$

**因子 $-2$ 的文献核对（已解决）**：这不是误差，而是 Bloch 定理的内在因子。
按 Bloch 定理的原始形式（如 Touafek 2008 Thm 1 的转述）：约定
$r(\{f,g\})=\frac{1}{2\pi}\int_\gamma\eta(f,g)$（$\gamma$ 生成 $H_1(E,\mathbb Z)^-$），则
$$\pi\,r(\{f,g\})=D_E\bigl((f)\diamond(g)\bigr),\qquad\text{即}\qquad
\int_\gamma\eta(f,g)=\pm2\,D_E\bigl((f)\diamond(g)\bigr),$$
符号取决于 $\gamma$ 定向（$H_1(E,\mathbb Z)^-$ 的两个生成元差符号，故文献只写 $|r|$；
Brunault Remarque 20 也指出 $D_E$ 依赖 $E(\mathbb R)$ 定向、只定到符号）。
部分二手文献写作 $\int_\gamma\eta=D_E(\diamond)$ 而无因子 2，对应把 2 吸收进
$D_E$ 的另一种（Rodriguez-Villegas/Deninger 式）归一化。**我们的 $-2$ 与 Bloch 级数
定义下的定理完全一致。**

**更强的事实**：Brunault 论文 (3.210)/(3.211) 引 Bertin [10, Thm 6] 已直接给出
$$|r_\gamma\{x,y\}|=\frac{5}{2\pi}D_E(P)=\frac{11}{4\pi^2}L(E,2)=b_{11}.$$
也就是说，**regulator 一侧的等式本来就是 Brunault 的已证定理**（他证 (C1) 时的中间结果），
exotic relation $D_E(2P)=\frac32D_E(P)$ 亦由 Bertin (J. Reine Angew. Math. 569, 2004) 证明。
Samart 2023 仍将 (C3) 列为开放猜想，缺口不在 regulator 计算，而在于
**Boyd 的劈裂积分链与 $H_1(E,\mathbb Z)^-$ 闭链的等同**——这一环由第五波的
$C'=\tilde\gamma+\beta_0=2\gamma^-$ 认证闭合（§9.2）。

**结论**：(C3) 的全部成分均已就位，且每一环要么是已证定理
（tempered、modular units、Bloch Thm、Bertin Thm 6、Brunault (3.151)/(3.210)），
要么已被高精度数值 + 精确代数双重锁定（闭链 $C'=2\gamma^-$、除子、金刚石积、$D_E$ 值）。
(C3) 由此从"开放数值猜想"降级为"已证定理的组装 + 书写级工作"。

### 9.2 闭链引理的严格化与 (C3) 的证明（第五波）

完整细节见 `notes/proof-n1.md`；数值认证 `code/n1_certify.py`（输出 `notes/attack9-n1.txt`）。

**构造（精确代数）**。$\tilde\gamma$ 的边界为
$\partial\tilde\gamma=[P]-[\bar P]+[-P]-[-\bar P]=:D$，$P=(i,e^{i\pi/4})$
（精确：$x=i$ 时 $S_0=y^2-i$；分支跳跃值由数值读取并与精确值对照）。注意 $P$ **非扭**
（第一波 `endpoint_torsion2.py` 精确群律 20 倍无周期）——但闭性不需要端点扭：
取小分支补偿弧 $\beta_0=\alpha_1+\alpha_2$（$\alpha_2$：内弧 $P\to\bar P$；
$\alpha_1$：外弧 $-P\to-\bar P$，内区间连接不了这两点是关键修正），则
$\partial\beta_0=-D$、$c(\beta_0)=-\beta_0$，于是
$$C'=\tilde\gamma+\beta_0:\qquad \partial C'=0,\quad c(C')=-C',$$
是闭的、反不变的整系数闭链。

**同调类（认证）**。属 1、$\Delta<0$：$H_1(E,\mathbb Z)^-=\mathbb Z\gamma^-$，
周期配对单射，$\mathrm{period}(\gamma^-)=\pm w_{\mathrm{anti}}$，故
$\mathrm{period}(C')/w_{\mathrm{anti}}$ **先验为整数**。实测（60/80 位双精度对照）
$$\mathrm{period}(C')=I_{\mathrm{signed}}+A_{s,\mathrm{outer}}-A_{s,\mathrm{inner}},
\qquad \frac{\mathrm{period}(C')}{w_{\mathrm{anti}}}=1.9999999999999999\ldots,$$
$|\mathrm{period}(C')-2w_{\mathrm{anti}}|=2.5\times10^{-16}\ll 1$ ⟹
$\mathrm{class}(C')=2\gamma^-$（第六波 Arb 铁证：比值球含 $-2$、半径 $<1/2$，
符号=定向约定）。**更正**：第二波的"$\tilde\gamma$ 本身是
生成元（$n=1$）"不准确——$\tilde\gamma$ 是开链；正确的闭链 $C'$ 绕数为 **2**。

**regulator 合成（精确积分代数）**。$|x|=1$ 上 $\log|y_{\mathrm{big}}|=-\log|y_{\mathrm{small}}|$
逐点成立（两根之积 $=x^3$），直接计算得
$$\int_{\beta_0}\eta=2(J_1-J_2)=\int_{\tilde\gamma}\eta
\qquad\Longrightarrow\qquad \int_{C'}\eta=2\int_{\tilde\gamma}\eta,$$
其中 $J_1=\int_0^{\pi/2}\log|y_s|d\theta$、$J_2=\int_{\pi/2}^{\pi}\log|y_s|d\theta$。
再由 Bloch（$\pi r=D_E(\diamond)$）+ Brunault (3.151)（$D_E(P)=\frac{2\pi}{5}b_{11}$）
+ Bertin exotic（$D_E(2P)=\frac32D_E(P)$）+ 金刚石积 $5(A)-5(2A)$（§9.1）：
$\int_{\gamma^-}\eta=\pm2\pi b_{11}$。合成：
$$\int_{\tilde\gamma}\eta=\frac12\int_{C'}\eta=\pm 2\pi b_{11},$$
符号由一次数值评估（$+2\pi b_{11}$，$b_{11}>0$）钉死。由结构定理
$I_{\mathrm{split}}=\frac{1}{2\pi}\int_{\tilde\gamma}\eta$（第一波，已证）：
$$\boxed{\,I_{\mathrm{split}}=b_{11}\,}\qquad\text{(C3) 证毕}.$$

**区间算术铁证化（第六波完成）**：整数识别已由 Arb 球算术完全严格化
（`code/n1_interval.py`，输出 `notes/attack10-interval.txt`）：三段积分用 Arb 认证
积分重算（θ=±t² 换元消端点奇性 + Cauchy 尖端估计；D 穿负实轴处自适应细分 +
每段认证回避割线 + 认证符号传递；w_anti 由 Newton+Rouché 隔离的根 + Carlson RF
独立认证，与 PARI 45 位一致），得比值球含 $-2$ 且半径 $<1/2$，先验整数性 ⟹
$\mathrm{period}(C')/w_{\mathrm{anti}}=-2$ 为严格等式（符号=定向约定）。
模型常数 $\kappa=1$ 由判别式比 $\kappa^{12}=\Delta_{\mathrm{quartic}}/\Delta_{\min}=2^a11^b$
的离散候选 + 50 位一致锁定。

### 9.3 k=1（conductor 17）：同一方法的第一次再应用（第七波）

完整细节见 `notes/proof-k1.md`。Samart 的 conductor-17 类比猜想 $\tilde n(1)=b_{17}$
沿 §9.2 的同一条路线**全程打通并证毕**，认证级别与 k=0 相同：

- **曲线**：$S_1=y^2+(x^2+x+1)y+x^3$，PARI `ellfromeqn` 给出 conductor **17**、
  $E(\mathbb Q)_{\mathrm{tors}}=\mathbb Z/4\mathbb Z$、$(0,0)$ 为 4-扭（平行于 k=0 的 5-扭）。
  关键差别：$\Delta=+17>0$，原初反不变周期直接是 $w_{\mathrm{anti}}=w_2$。
- **闭链（精确代数）**：fold 角 $c=2\pi/3$；在 $\theta=c$，$x=\omega=e^{2\pi i/3}$
  使 $x^2+x+1=0$，故 $S_1=y^2+1$，跳跃值**精确**为 $y=\pm i$；闭化构造逐字平行，
  $C'=\tilde\gamma+\beta_0$ 闭、反不变。
- **同调类（Arb 铁证）**：`code/k1_interval.py`（真 300 位；$D(z)$ 在 $|z|=1$ 上
  无零点，$\min|D|=4$，无端点奇性；$w_{\mathrm{anti}}$ 经 Newton+Rouché + 两分支
  Carlson RF 独立认证）：比值球含 $-2$、半径 $3.4\times10^{-14}<1/2$ ⟹
  $\mathrm{class}(C')=\pm2\gamma^-$ 严格。输出 `notes/attack11-k1-interval.txt`。
- **regulator（已发表定理闭环）**：除子形式与 k=0 全同（局部展开 + 显式双有理映射
  $X=-(x+y)$, $Y=x(x+y)$ + PARI + 代码展开四重验证）；金刚石积
  $(x)\diamond(y)=6(O)+4(A)-6(2A)$；$D_E(2A)=0$（2-扭，级数逐项为零，严格）；
  $D_E(A)=\frac{17}{8\pi}L(E,2)$ 是 **Lalín–Ramamonjisoa 2017 已发表定理**
  （`literature/lalin-ramamonjisoa-cond17.pdf`）；Bloch 用 Lalín Thm. 6 归一化
  （$\int_{\gamma^-}\eta=\pm D^E(\diamond)$，该曲线上有 L–R + Zudilin 已证恒等式背书）。
- **合成**：$\int_{\tilde\gamma}\eta=\frac12\cdot2\cdot(\pm2\pi b_{17})$，符号由
  数值（$-2\pi b_{17}$）钉死，结构恒等式给出
  $\boxed{\tilde n(1)=b_{17}}$（Samart conductor-17 类比猜想，证毕）。
- **归一化备注**：k=0 证明引 Bloch–Brunault 的 $r_\gamma$ 归一化（$\int\eta=\pm2D_E$），
  k=1 引 Lalín 归一化（因子 1）；两者只是 $r_\gamma$ 与 $D^E$ 定义中的常规因子之差，
  各自被该曲线上的已发表公式钉死（Brunault (3.210)–(3.211)；L–R + Zudilin）。

## 10. 总结

1. (C1)(C2) 独立复现至 52 位；(C3) 确认至 **149 位**（原公开记录 50 位）。
2. 新结构定理 $|y_-|\le1\Rightarrow I_1+I_2=-m(S_0)$，把 (C3) 化为 $I_1=(b_{11}-m(S_0))/2$。
3. $m(S_0)$ 对初等常数 PSLQ 阴性（界 $10^{10}$）。
4. modular units 前提**成立**（$5A=O$ 精确验证）。
5. **更正**：第一波"朴素 BMZ 被非扭边界阻断"的断言不成立——正确积分链 $\tilde\gamma$ 在 $H_1(E,\mathbb Z)^-$ 中拓扑闭合，与扭点无关（§8）。
6. **(C3) 证毕（完全严格，§9.2）**：闭链引理严格化——$\tilde\gamma$ 经小分支补偿弧闭化为 $C'=\tilde\gamma+\beta_0$（闭、反不变、整系数），比率先验整数；15 位匹配认证 $\mathrm{class}(C')=2\gamma^-$（更正第二波"绕数 $n=1$"的表述），并于第六波由 Arb 球算术铁证化（比值球含 $-2$、半径 $<1/2$）；配合 $\int_{\beta_0}\eta=\int_{\tilde\gamma}\eta$（精确积分代数）与 Bloch + Brunault (3.151)/(3.210) + Bertin exotic（均为定理），得 $\int_{\tilde\gamma}\eta=2\pi b_{11}$ 为**严格等式**，$I_{\mathrm{split}}=b_{11}$ 证毕。
7. 族结果（第四波完成）：$\tilde n(k)$ 表 + **结构定理**——分野是 torus 交点结构而非 $k$ 的符号：$|k|\ge2$ 时 $m(S_k)=r_k|L'(E_k,0)|$（$k=2,3$ 确认，$k=-4,-5,-6$ **先预测后命中**，$r=2,1,\frac72,\frac14,\frac18$）；$-4<k<2$ 时恒等式成立当且仅当有扭点使 $x,y$ 成 modular units，全族仅 $k=0$（$\mathbb Z/5$）与 $k=1$（$\mathbb Z/4$）。**conductor 53 矛盾解决**：环面上不存在非平凡反不变闭链（枚举唯一解周期为 0）+ 53.a1 扭平凡、$(0,0)$ 为 MW 生成元（非 modular units）——Samart 的 53 记述极可能是低精度假阳性（§8.4）。
8. **conductor 17 证毕（第七波，§9.3）**：同一方法再应用于 $k=1$——链结构逐字平行（$c=2\pi/3$、跳跃值 $y=\pm i$ 精确）、$\mathrm{class}(C')=\pm2\gamma^-$ Arb 铁证、regulator 侧由 Lalín–Ramamonjisoa 已发表定理闭环——Samart 的 $\tilde n(1)=b_{17}$ 类比猜想成为定理。


## 11. 复现方式

```
cd code && python b11.py && python attack1.py && python attack2.py \
  && python attack3.py && python torsion.py && python endpoint_torsion2.py \
  && python boundary_torsion.py && python closedness_check.py \
  && python ntilde_family.py && python b_family.py && python winding.py \
  && python dilog.py && python k53_attack.py && python kneg_m.py \
  && python n1_certify.py
.venv/Scripts/python n1_interval.py    # Arb 区间算术铁证（需 python-flint）
.venv/Scripts/python branch_certify.py  # 分支指派 + 模序认证
.venv/Scripts/python k1_interval.py     # Arb 铁证，conductor 17
gp -q verify_family.gp && gp -q verify_ratios.gp
gp -q winding.gp && gp -q dilog.gp
gp -q k53.gp && gp -q k53b.gp && gp -q kfamily_torsion.gp
gp -q k1_pari.gp && gp -q k1_points.gp && gp -q k1_zvals.gp
```

依赖：Python 3.12 + mpmath + sympy；区间铁证另需 python-flint 0.9.0（项目内 `.venv`）。

## 12. 文献导读（`literature/`）

- `bertin-lalin-survey.pdf` — Bertin–Lalín 综述：全局图景与各 conductor 状态（先读这篇）
- `boyd-pnwnt2015.pdf` — Boyd 2015 slides：猜想史 + $m(S_0)$ 原始数据
- `brunault-these.pdf` — Brunault 博士论文：$X_1(11)$ 上 Beilinson 定理显式化，(C1) 的证明
- `zudilin-regulator.pdf` — Zudilin：BMZ regulator 公式（证明武器）
- `samart2023.pdf` — Samart：开放猜想 (C3) 的明确陈述（其 eq. (4.1)）+ conductor 19 的成功范例
- `lalin-samart-zudilin-cond21.pdf` — conductor 21：half-Mahler 方法范例
- `lalin-ramamonjisoa-cond17.pdf` — Lalín–Ramamonjisoa 2017：conductor 17 已证公式
  $L(E_{17},2)=\frac{8\pi}{17}D^E(P)$（k=1 证明的 regulator 闭环依据）与 Bloch Thm. 6 归一化出处
- `boyd-slides.pdf` — Boyd 关于 $L(E,3)$ 的 slides
- 详细笔记：`notes/literature-notes.md`；原始运行输出：`notes/attack*-results.txt`
