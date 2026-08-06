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

（我们的 `code/b11.py` 用 $f_{11}$ 的系数和"近似函数方程"把它算到 800 位；
$L(E,2)=\sum a_n[\,\cdots\,]$ 以 $e^{-2\pi n/\sqrt{11}}$ 速度收敛，尾项界驱动截断：442 项 @ 800 dps。）

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

$m(S_0)$ 与 $b_{11}$ 之间**未知任何有理关系**：我们用 PSLQ 在 $10^8$ 系数界内未发现（Boyd 亦报告了同样的阴性结果）。
但 Boyd 同时发现：**沿 branch cut 劈开的带符号积分**仍然等于 $b_{11}$。写 $S_0=0$ 的两根

$$y_\pm(x)=-\frac{x^2+1}{2}\pm\sqrt{\frac{(x^2+1)^2}{4}-x^3},$$

$x=e^{i\theta}$ 沿上半环面走，$\theta=\pi/2$ 处正是 torus 交点 $x=i$（"分支切口"），定义

$$I_{\mathrm{split}}:=\underbrace{\frac1\pi\int_0^{\pi/2}\log|y_-(e^{i\theta})|\,d\theta}_{I_1}
\;-\;\underbrace{\frac1\pi\int_{\pi/2}^{\pi}\log|y_-(e^{i\theta})|\,d\theta}_{I_2}\ \stackrel?=\ \pm b_{11}. \tag{C3}$$

（符号取决于哪个根叫 $y_-$：$|y_+||y_-|=|x^3|=1$，换根整体变号。Samart 的 eq. (4.1) 把右端写成 $-L'(E,0)$；但按 Samart 自己的 $y_-$ 定义，左端积分是正数 $+0.1521471\ldots$，该负号疑似笔误——按我们的 $y_-$ 约定恒等式读作 $+b_{11}$，(C3) 因此记 $\pm$。）

**$y_-$ 的约定（澄清）**：$y_-$ 指 $|x|=1$ 上小模长根的**连续分支**：在 $\theta=0$ 处两根合并于 $y_-(1)=-1$；在 $(0,\pi)\setminus\{\pi/2\}$ 上由 $|y_-|<|y_+|$ 唯一确定；在第二个节点 $\theta=\pi/2$ 处按连续性延拓（两个节点处 $|y_\pm|=1$，故该处 $\log|y_-|=0$，被积函数不受取值选择的影响）。

Boyd 的原话（经 Samart 2023 §4 引用）：

> "This is in accord with our contention that in case $P$ vanishes on the torus, it is the integral
> of $\omega$ around a branch cut rather than $m(P)$, which should be rationally related to $L'(E,0)$."

Samart 2023（arXiv:2301.05390）以 Boyd 的“$=$?”猜想记号把 (C3) 记录为**开放恒等式**（其 eq. (4.1)），并指出可尝试用他在
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
  $E_1$ 为指数积分（来自 $\int_1^\infty e^{-ty}/y\,dy$）。项衰减 $\sim e^{-1.894n}$；尾项界驱动截断：442 项 @ 800 dps（旧稿"200 项够 300 位"有误，已更正）。
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

### 6.2 开放猜想 (C3) 确认到 366 位 —— **主要数值结果**（`code/attack13_c3_300.py`，存档 `notes/attack13-c3-300.txt`）

$$I_{\mathrm{split}}=0.152147141725918049486227297478634495628143589164226122809889823882023289695302776676\ldots$$
$$|I_{\mathrm{split}}-b_{11}|=9.26\times10^{-367}.$$

（800 dps 工作精度，641 秒；$b_{11}$ 的 $\eta$-级数值与 PARI/GP `lfun(E,0,1)` **330 位逐位吻合**——三条独立计算共享 310 位公共前缀。原 `attack3.py` 的 149 位结果标记 superseded。）
此前公开记录是 Boyd 的 50 位验证（Boyd 2015 slides，p. 28，即 `literature/boyd-pnwnt2015.pdf`）；本次推进到 **366 位**。

### 6.3 结构恒等式（先数值发现，后给出证明）

计算中注意到 $I_1+I_2=-m(S_0)$ 吻合到 152 位。事实上这是**定理**：

> **命题（区间算术证明）**：在 $[0,\pi]$ 上 $|y_-(e^{i\theta})|\le1\le|y_+(e^{i\theta})|$（球算术自适应二分认证：$\log|y_-|$ 在 $(0,\pi)$ 上为负，仅在折点 $\theta=\pi/2$ 与端点 $\theta=0$ 处为零；`code/branch_certify.py`——本命题是 §9.2 意义上的计算机辅助结果）。又 $|y_+y_-|=|x^3|=1$，故
> $$m(S_0)=\frac1\pi\int_0^\pi\log|y_+|\,d\theta=-\frac1\pi\int_0^\pi\log|y_-|\,d\theta=-(I_1+I_2).\qquad\blacksquare$$

**推论**：(C3) 等价于 $I_1=\dfrac{b_{11}-m(S_0)}{2}$、$I_2=-\dfrac{b_{11}+m(S_0)}{2}$。
也就是说，劈裂积分的猜想给出的是"大弧段积分"与"小弧段积分"各自的确切值。

### 6.4 $m(S_0)$ 的负结果

$m(S_0)=0.40560295591501040\ldots$（Boyd 的 $0.4056029$ ✓）。

- PSLQ$(m(S_0),b_{11})$，系数界 $10^8$：**未发现关系**（复核 Boyd 的 "seemingly not $rb_{11}$"）；
- PSLQ 对 $\{m(S_0),b_{11},\log2,\log3,\mathrm{Catalan},m(1+x+y)\}$，系数界 $10^{10}$：**未发现关系**。
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
群律精确计算（`code/torsion.py`，在四次模型上进行、以其无穷远点为群单位元）：对 $S_0$ 上的点 $A=(0,0)$（即四次坐标下的点 $(x,u)=(0,1)$），

$$2A=(0,-1),\qquad 4A=(1,0)=-A\quad\Longrightarrow\quad \boxed{5A=O}.$$

于是 $\operatorname{div}(x)=[(0,0)]+[(0,-1)]-2[P_\infty]$ 与
$\operatorname{div}(y)=3[(0,0)]-3[P_\infty]$ 的支撑全是 5-扭点。这些恰为 $X_1(11)$ 的有理尖点：
$S_0=0$ 与 $E:y^2+y=x^3-x^2=X_1(11)$ 的双有理识别是精确的（Riemann–Roch 计算复合 PARI 精确
极小化变换，往返验证，`code/kappa_exact.py`），且在标准模同构（无穷远尖点 $\mapsto O$，
$\omega_f=dx/(2y+1)$）下有理尖点 $P_v$（$v\in(\mathbb Z/11\mathbb Z)^\times/\pm1$）的像为
$$P_1=\infty\mapsto O,\quad P_2\mapsto(1,0)=3A,\quad P_3\mapsto(0,-1)=4A,\quad P_4\mapsto(0,0)=A,\quad P_5\mapsto(1,-1)=2A,$$
即恰好映满 $E(\mathbb Q)=\mathbb Z/5\mathbb Z$（Brunault Thm. 8 证明，(3.152)）。故除子支撑尖点化，
**$x,y$ 确为 modular units**（Manin–Drinfeld）。tempered 性同样精确化：Newton 面多项式
$x^3+y$、$x^3+x^2y$、$x^2y+y^2$、$y^2+y$ 全分圆给出 $\{x,y\}\in K_2(E)\otimes\mathbb Q$
（Rodríguez Villegas 判据）；更精确地，除子支撑处的 tame 符号由局部展开精确算出
（`code/bertin_diamond.py`）：$T_v\{x,y\}=\pm1$（$v\in\{O,A,2A,3A\}$），故 $2\{x,y\}$ 的
tame 符号已平凡，$\eta(x,y)$ 的实留数 $\pm\log|T_v\{x,y\}|$ 处处为零。

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

**命题（环面交，精确刻画）**：对**实** $k$，记 $x=e^{i\theta}$（$\theta\in[-\pi,\pi]$）。环面 $|x|=|y|=1$ 上的点 $(x,y)$ 落在 $S_k=0$ 上当且仅当 $\theta$ 落入下列两情形之一：
(i) $\theta=0$ 且 $|k+2|\le2$；(ii) $k+2\cos\theta=0$（可解 ⟺ $|k|\le2$）。
由此：(1) 环面与曲线相交 ⟺ $k\in[-4,2]$；(2) $-2\le k\le2$ 时出现**折点** $\theta=\pm\arccos(-k/2)$（该处 $B=x^2+kx+1=0$，$y^2=-x^3$ 的两根模长均为 1）；(3) $-4\le k\le0$ 时另有 $\theta=0$ 处的**分支交换交点**（$S_k(1,y)=y^2+(k+2)y+1$ 的两根都在单位圆上——$-4<k<0$ 时相异，$k=-4$（$y=1$）与 $k=0$（$y=-1$）处为重根即切触）；(4) $k=2$ 时折点与边界 $\theta=\pi$ 合并（$S_2(-1,y)=y^2-1$）。这与 Samart 观察到的环面相交参数区间 $K\cap\mathbb R=[-4,2]$ 精确吻合。

*证明*：记 $y=e^{i\phi}$。由 $x+x^{-1}=2\cos\theta$ 得 $B=x(k+2\cos\theta)$；方程 $e^{2i\phi}+Be^{i\phi}+e^{3i\theta}=0$ 除以 $e^{i\phi}$，并用 $e^{i\phi}+e^{i(3\theta-\phi)}=2\cos\psi\,e^{3i\theta/2}$（$\psi:=\phi-3\theta/2\in\mathbb R$），得
$$2\cos\psi\,e^{i\theta/2}=-(k+2\cos\theta),$$
其右端为实数。比较虚部与实部：
$$\cos\psi\,\sin(\theta/2)=0,\qquad 2\cos\psi\,\cos(\theta/2)=-(k+2\cos\theta).$$
第一式给出两个情形。(i) $\sin(\theta/2)=0$，即 $[-\pi,\pi]$ 上 $\theta=0$：第二式为 $2\cos\psi=-(k+2)$，可解 ⟺ $|k+2|\le2$，此时 $y=e^{i\phi}$（$\cos\phi=-(k+2)/2$）确为 $S_k(1,y)=y^2+(k+2)y+1$ 的根。(ii) $\cos\psi=0$：第二式迫使 $k+2\cos\theta=0$，可解 ⟺ $|k|\le2$；反之若 $k+2\cos\theta=0$，则 $B=0$、$S_k(e^{i\theta},y)=y^2+e^{3i\theta}$ 有两根 $y=\pm e^{i(3\theta+\pi)/2}$（模长均 1）——即折点，$\psi=\pm\pi/2$ 实现该式。两参数区间的并为 $[-4,0]\cup[-2,2]=[-4,2]$，得 (1)；(2)–(4) 直接读出（$k=2$ 时折点 $\arccos(-1)=\pi$；$k=-2$ 时折点 $\arccos(1)=0$ 与情形 (i) 的点重合）。$\blacksquare$

| $k$ | $N$ | $w$ | 折点 $c/\pi$ | $\tilde n(k)$ vs $b_N=\lvert L'(E_k,0)\rvert$ | 证据标签 |
|---|---|---|---|---|---|
| $-3$ | $83$ | $-1$ | 无 | 比值 $0.8529175\ldots$ | 未发现 (PSLQ) |
| $-2$ | $91$ | $-1$ | 无 | 比值 $0.6339454\ldots$ | 未发现 (PSLQ) |
| $-1$ | $53$ | $-1$ | $1/3$ | 比值 $0.7392026\ldots$ | 机制失效（命题，见下） |
| $0$ | $11$ | $+1$ | $1/2$ | $\tilde n=-b_{11}$ | **已证**，(C3)，366 位 |
| $1$ | $17$ | $+1$ | $2/3$ | $\tilde n=+b_{17}$ | **已证**（cond.），定理（附录 A）|
| $2$ | $37$ | $-1$ | 无 | $m=\tilde n=2\,b_{37}$ | 数值，70 位 |
| $3$ | $79$ | $-1$ | 无 | $m=\tilde n=b_{79}$ | 数值，70 位 |
| $-4$ | $37$ | $-1$ | 无（相切） | $m=\frac72\,b_{37}$ | 数值，预测，25 位 |
| $-5$ | $359$ | $-1$ | 无 | $m=\frac14\,b_{359}$ | 数值，预测，25 位 |
| $-6$ | $997$ | $-1$ | 无 | $m=\frac18\,b_{997}$ | 数值，预测，25 位 |

（"数值" 行是猜想级恒等式（见 §8.4 末与 §10）；"未发现" 指系数界 $10^8$ 的 PSLQ 搜索无有理关系。）

要点：

- **$k=1$（conductor 17）**：$\tilde n(1)=b_{17}$——Samart 所述 conductor 17 "(4.1) 类似猜想"
  的独立高精度确认（$y_-$ 约定下积分 $=-L'(E_{17},0)$，有理因子 $r=1$）。
- **$k=2,3$（conductor 37、79，根数 $-1$）**：无 torus 交点，Mahler 测度本身满足 Boyd 型恒等式
  $m(S_2)=2|L'(E_{37},0)|$、$m(S_3)=|L'(E_{79},0)|$，70 位（PARI `lfun`）。这是本次意外收获。
- **$k=-1$（conductor 53）：机制失效判定（第四波，`code/k53_attack.py` + `k53.gp`；第十波精确化 `code/k53_smith.py`）**：Samart 2023 §4 关于
  conductor 53 的"类似猜想"只有一句话（无公式、无精度、无引用）。既然待检验的精确陈述不在案，
  我们无法直接裁决它；我们能够且确实检验的是：本报告的**机制**——环面上闭反不变闭链配对
  modular-unit symbol——能否延伸到 $k=-1$。答案是否定的，有三个相互独立的原因：
  1. **拓扑障碍（精确）**：$k=-1$ 时 $|x|=1$ 上的 torus 交点为 $\theta=0$（两根 $e^{\pm2\pi i/3}$，
     分支在此**交换**）与 $\theta=\pm\pi/3$（$y=\pm1$）。以交点为断点、big/small
     分支组合的**闭反不变链**空间为 $\mathbb Z\cdot(1,1,1,1)$：反不变对称化的边界矩阵是
     $\pm1/0$ 整数矩阵，其在 $\mathbb Q$ 上的核**恰为一维**，由精确有理线性代数算出
     （`code/k53_smith.py`，非系数盒搜索）；生成元 $(1,1,1,1)$ 的周期为 $0$——而且是
     **精确**为零，非数值：$(1,1,1,1)$ 是两叶在全圆上的和，两分支
     $u=2y+B=\pm\sqrt{B^2-4x^3}$，故 $1/u_++1/u_-\equiv0$ 逐点成立，$\sum dx/u$
     恒等于零（连续根沿圆走**两圈**才闭合，该两圈回路由同一抵消而为零）。因此
     **环面上不存在非平凡的反不变闭链**，$H_1(E,\mathbb Z)^-$ 生成元无法在 torus 上实现，
     Boyd 型环面积分对象根本不存在。
  2. **代数障碍（精确除子论证）**：函数 $x,y$ 不是 modular units。事实上
     $\operatorname{div}(x)=[P]+[-P]-2[O]$、$\operatorname{div}(y)=3[P]-3[O]$，其中
     $P=(0,0)$ 是 Mordell–Weil **生成元**：53.a1 有理扭子群平凡、秩 1
     （PARI `ellorder`$(P)=0$、`ellidentify`）。该曲线是其 conductor 的 strong Weil curve，
     $X_0(53)$ 的两个尖点映到 $E$ 的扭点（Manin–Drinfeld），而扭子群平凡故都映到 $O$——
     于是 $E$ 上的 modular unit 除子只能支撑于 $\{O\}$，而 $\operatorname{div}(x)$、
     $\operatorname{div}(y)$ 的支撑含非扭点 $P$。Beilinson–Brunault 机制因此没有可作用的
     symbol，任何闭链的 regulator 配对都没有理由是有理数 $\times\,b_{53}$。
  3. **数值观察**：自然候选对象（半环 $L_1$，连续根一圈，开链）周期与 $w_{\mathrm{anti}}$
     之比 $0.5492906\ldots$ 在已算精度内未见整数性，regulator 积分与 $2\pi b_{53}$ 之比
     $0.7392026\ldots$；系数界 $10^8$ 的 PSLQ 搜索未发现有理关系。

  这些不排除经其他途径得到的恒等式；被排除的是证明 (C3) 的那个特定机制。
  在缺乏精确候选公式的情况下，conductor 53 情形**留作开放**。
- **扭点对照表**（`code/kfamily_torsion.gp`，PARI `elltors`/`ellorder`）揭示真正的分野：

  | $k$ | $N$ | 扭子群 | $\operatorname{ord}(0,0)$ | Boyd 型恒等式 |
  |---|---|---|---|---|
  | $0$ | $11$ | $\mathbb Z/5$ | $5$ | 已证（modular units，(C3)） |
  | $1$ | $17$ | $\mathbb Z/4$ | $4$ | 已证（cond.；tempered + 扭点支撑，附录 A） |
  | $-3,-2,-1,2,3$ | $83,91,53,37,79$ | 平凡 | $\infty$ | 非模单位（命题，见下） |

- **修正后的结构二分——精确结构分析 + 数值证据（第四波）**：分野不是 $k$ 的符号，而是 **torus 交点结构**：
  - **族中的模单位与 temperedness（精确）**：(i) $k=0$（$E=X_1(11)$，
    $E(\mathbb Q)_{\mathrm{tors}}=\mathbb Z/5$）时 $x,y$ 的除子支撑于有理扭点，而这些恰是模识别下
    的有理尖点（§7），故 $x,y$ 是 modular units；此情形 Boyd 型恒等式**已证**（(C3)）。
    $k=1$（$\mathbb Z/4$）时 $x,y$ 的除子支撑于有理扭点且 symbol $\{x,y\}$ tempered（Newton 面
    多项式全分圆，附录 A），故 $\{x,y\}\in K_2(E)\otimes\mathbb Q$——这是 conductor-17 证明
    唯一用到的 K 理论输入；我们**不**断言 $k=1$ 时 $x,y$ 是 modular units：在自然模模型
    $X_0(17)$（仅两个尖点）上，支撑于全部四个有理扭点的除子不可能尖点化，且未提供其他模识别。
    (ii) $k\in\{-3,-2,-1,2,3\}$ 时 $E_k$ 有理扭平凡（PARI `elltors`，精确），
    $(0,0)\ne O$ 非扭；每条 $E_k$ 是其 conductor 的 strong Weil curve，$X_0(N)$ 的尖点映到
    $E_k$ 的扭点（Manin–Drinfeld）即映到 $O$，modular unit 的除子只能支撑于 $\{O\}$，而 $x=0$
    与 $S_k$ 交于 $(0,0)$ 与 $(0,-1)$。故这些 $k$ 的 $x,y$ **不是** modular units，模单位
    regulator 机制不适用；$k=-1$ 时闭链成分也失效（上条命题）。
  - **观察到的二分（数值证据，与环面交命题联合）**：整数 $k\notin(-4,2)$（无真正 torus 交点）时
    标准 Deninger 机制预期直接适用，$m(S_k)=r_k\,|L'(E_k,0)|$（$r_k$ 为小有理数）在所有已算情形
    获数值确认——$k=2,3$ 已确认，$k=-4,-5,-6$ 为**先预测后命中**：
    $m(S_{-4})=\frac72|L'(E_{37},0)|$、$m(S_{-5})=\frac14|L'(E_{359},0)|$、
    $m(S_{-6})=\frac18|L'(E_{997},0)|$（25 位精确）；注意 $S_{-4}$ 与 $S_2$ 同为 conductor 37，
    给出同一曲线的 Rodriguez-Villegas 型关系。整数 $k$ 且 $-4<k<2$ 时，$x,y$ 为模单位的情形
    只有 $k=0$；$k=1$ 有 tempered symbol 与扭点支撑的除子（但无模单位结构，见上条）；
    $k=-1,-2,-3$（扭平凡、$\theta=0$ 处分支交换交点）在 PSLQ 界
    $10^8$ 内未发现有理关系。我们把这些记录为**猜想的证据，而非定理**。


## 9. 证明纲要：Beilinson–Brunault 路线

**主定理**：$I_{\mathrm{split}}=b_{11}$（(C3) 中的带符号劈裂积分；$b_{11}=L'(E_{11},0)$）。

证明结合精确代数与计算机辅助步骤；除以下清单外的一切步骤都是精确的（有理或符号算术）：

1. **模序与分支指派**：结构命题（§6.3）的模长排序与闭链的端点分支指派，球算术证明（§9.2）；
2. **同调类** $\mathrm{class}(C')=2\gamma^-$：由先验整数 $\mathrm{period}(C')/w_{\mathrm{anti}}$ 经球算术证明（§9.2）；
3. **regulator 积分的符号**：由严格含于半直线的认证区间包围证明（§9.2）；
4. **高精度浮点求值仅作一致性核查**（如与 $b_{11}$ 的 366 位吻合），从不作为证明步骤。

约定：$\gamma^-$ 是反不变生成元引理（§9.2）的**本原**生成元；
$\eta(x,y)=\log|x|\,d\arg y-\log|y|\,d\arg x$；$\mathbb Z[E]^-$ 取 Lalín–Ramamonjisoa 的
商约定；$D_E$ 是背景部分的椭圆双对数，按 L–R Def. 5 eq. (10) 逐字实现。

| 步骤 | 内容 | 状态 |
|---|---|---|
| S1 | 闭性引理：$\tilde\gamma$ 闭化为 $C'=\tilde\gamma+\beta_0\in H_1(E,\mathbb Z)^-$，$\mathrm{class}(C')=2\gamma^-$（§9.2，比率先验整数 + 15 位认证） | **已证（认证）** |
| S2 | tempered：$S_0$ 的 Newton 面多项式 $x^3+y$、$x^3+x^2y$、$x^2y+y^2$、$y^2+y$ 全分圆，故 $\{x,y\}\in K_2(E)\otimes\mathbb Q$ | 已查 |
| S3 | modular units：$x,y$ 在 $E=X_1(11)$ 上的除子支撑于尖点（$5A=O$ 精确验证） | 已证 |
| S4 | Beilinson–Brunault regulator 定理 + Brunault (3.151)：$L(E,2)=\frac{10\pi}{11}D_E(P)$，系数 40 位复核 | 已核（文献+数值） |
| S5 | regulator 常数：Brunault 的 Siegel 单位 regulator 定理（JNT 163 (2016) Thm. 1，已证）直接计算 $\int_{\gamma^-}\eta=\pm2\pi b_{11}$——四步精确链（尖点↔扭点、Siegel 表示、Manin 分解、$F_{\mathrm{total}}=-2f_{11}$），§9.1 锚定定理 | 闭合 |

### 9.1 regulator 常数的显式计算（第三波，40 位）

第三波把 S5 的"余下代数计算"完成了（代码 `code/dilog.gp` + `code/dilog.py`，
原始输出 `notes/attack7-dilog.txt`）。全部成分如下。

**除子代数**（三次模型上的射影计算，Abel 主除子检验通过；**更正**：此前误写
$\operatorname{div}(y)=3[A]-2[O]-[3A]$——$S_0$ 的射影闭包有**两个**无穷远点
$O_1=[0:1:0]$、$Q_\infty=[1:-1:0]$（群单位元是 $Q_\infty$），分别映到 $3A$ 与 $O$，
旧写法把二者混淆；$y$ 在 $A$ 处取值 $-1\ne0$，三重零点其实在 $2A$）：
$$\operatorname{div}(x)=[A]+[2A]-[O]-[3A],\qquad \operatorname{div}(y)=3[2A]-2[3A]-[O],$$
其中 $A=(0,0)$ 即 Brunault 记号中的 5-扭点 $P$，$O=[0:1:0]$ 为 11.a3 的群单位元。
金刚石积（修正后）
$$(x)\diamond(y)=6(O)-5(A)+5(2A),\qquad\text{故}\quad D_E\bigl((x)\diamond(y)\bigr)=-5D_E(P)+5D_E(2P).$$

**椭圆双对数**（mpmath 60 位；$\Delta<0$ 故取格基 $(w_1,\,w_1-w_2)$ 使 $\tau=0.5+0.2299i$
落在上半平面，$q=e^{2\pi i\tau}\approx-0.2359$ 为负实数；$z(P)=0.6w_1$、$z(2P)=0.2w_1$ 由 PARI `ellpointtoz` 给出）：
$$D_E(P)=0.1911937370843316957549544343121738161012\ldots$$
两条关键恒等式（均 40 位）：

- **Brunault 系数**：$D_E(P)=\dfrac{2\pi}{5}\,b_{11}=\dfrac{11}{10\pi}L(E,2)$，
  即 Brunault 论文 (3.151) 中系数的精确形式 $L(E,2)=\frac{10\pi}{11}D_E(P)$，数值钉死；
- **exotic relation**：$D_E(2P)=\frac32\,D_E(P)$（Bertin 已证此类关系；
  注意是 $3/2$ 而非朴素的 2 倍）。

**闭环**：代入 exotic relation 得
$$-5D_E(P)+5D_E(2P)=-5\Bigl(1-\frac32\Bigr)D_E(P)=\frac52\,D_E(P)=\pi\,b_{11},$$
而绕数一节已锁定 $\int_{\tilde\gamma}\eta(x,y)=2\pi b_{11}$（60 位，PARI 交叉验证）。

**regulator 的锚定（第十波重修；rev4：Siegel 单位定理直接计算）**：$\{x,y\}$ 沿
$\gamma^-$ 的 regulator 评估现在**直接**进行——用 Brunault 已证的 Siegel 单位
regulator 公式（J. Number Theory **163** (2016) 542–569，Thm. 1，
`literature/brunault-siegel.pdf`）；Bloch 菱形定理不再使用。计算所引出的因子-2 澄清
（指数引理与下文讨论）仅是文献状态说明，同样不属于证明。

**指数引理（regulator 生成元与指数）**：设 $E/\mathbb R$ 为椭圆曲线，$c$ 为复共轭
对合，$H_1(E,\mathbb Z)^\pm=\ker(c_*\mp1)$，$\{f,g\}$ 为 $E$ 上的 tempered symbol。
配对 $\gamma\mapsto\int_\gamma\eta(f,g)$ 下降到**余商**
$H_1(E,\mathbb Z)/H_1(E,\mathbb Z)^+$ 上（本文中"余商"（coinvariant quotient）专指
整同调模掉**不变子格** $H_1(E,\mathbb Z)^+$ 的这个商；我们不把 $\ker(c_*+1)$ 与该商
混用——二者相差的指数正是这里的关键）：

1. 若 $\Delta>0$，则 $H_1(E,\mathbb Z)^-=\mathbb Z b$ 同构地映上余商 $\mathbb Z[b]$：
   反不变生成元代表余商生成元（指数 $1$）；
2. 若 $\Delta<0$，则 $H_1(E,\mathbb Z)^-=\mathbb Z\gamma^-$ 在余商 $\mathbb Z[b]$ 中的
   像为 $2\mathbb Z[b]$（指数 $2$）：用反不变生成元引理（§9.2）的基，
   $\gamma^-=a-2b\equiv-2[b]$，故对每个 tempered symbol
   $$\int_{\gamma^-}\eta(f,g)=-2\int_{[b]}\eta(f,g).$$

*证明*：配对在 $H_1(E,\mathbb Z)^+$ 上平凡：若 $c(\gamma)=\gamma$，则由 $c^*\eta=-\eta$，
$$\int_\gamma\eta=\int_{c(\gamma)}\eta=\int_\gamma c^*\eta=-\int_\gamma\eta$$
（参见 Lalín–Ramamonjisoa Rem. 4），故下降到余商。$\Delta<0$ 时 $c_*:a\mapsto a$、
$b\mapsto a-b$（反不变生成元引理），故 $H_1^+=\mathbb Z a$、余商 $=\mathbb Z[b]$，
反不变类 $a-2b$ 投影为 $-2[b]$，指数 $2$；$\Delta>0$ 时 $c_*:a\mapsto a$、$b\mapsto-b$，
故 $H_1^-=\mathbb Z b$，投影 $b\mapsto[b]$ 为同构。$\blacksquare$

**锚定：Siegel 单位直接计算**。symbol $\{x,y\}$ 是一对 **modular units**（§9.1
modular units 段：其除子支撑于 $X_1(11)$ 的有理尖点）。Siegel 单位的 regulator 积分
由一条已证公式计算——既不涉及金刚石积、也不涉及同调格的选取：对 Siegel 单位
$g_u,g_v$，$u=(a,b)$、$v=(c,d)\in(\mathbb Z/N\mathbb Z)^2\setminus0$，
$$\int_0^{i\infty}\eta(g_u,g_v)
=\pi\,\Lambda^*\big(e_{a,d}e_{b,-c}+e_{a,-d}e_{b,c},\,0\big)\qquad\text{(S)}$$
其中 $e_{a,b}$ 是 Brunault Thm. 1 eq. (2) 的显式权 1、level $N^2$ Eisenstein 级数，
$\Lambda^*$ 是其 §2 的正则化完全 $L$ 值；任意模符号由线性处理（其 Rem. 2）。
（(S) 中的因子 $\pi$ 已对原文 LaTeX 源码核实——PDF 文本抽取会吞掉行内 $\pi$；我们对
公式的实现以三种独立方法 60 位复现了原文 §5.1 的 conductor-14 应用而验证，存档
`notes/attack16-siegel-anchor.txt`。）

**定理（regulator 锚定）**：取反不变生成元引理（§9.2）的 $\gamma^-$，则
$$\int_{\gamma^-}\eta(x,y)=\pm2\pi b_{11}\qquad\text{(A)}$$
（符号取决于定向）。

*证明*：全部步骤为精确有理/$q$ 展开计算或已存档（`code/siegel_anchor_step1`–`11`，
`notes/attack16-siegel-anchor.txt`）。
**(1) 尖点与除子**：在附属于 $f_{11}$ 的模参数化 $X_1(11)\to E$（无穷尖点
$\mapsto O$）下，有理尖点 $k/11\mapsto m_kA$，$(m_1,\dots,m_5)=(0,2,1,4,3)$，与
Brunault 的尖点表（(3.152)–(3.153)）一致；尖点 $i\infty$ 与 $1/11$ 在 $X_1(11)$ 上
重合（因 $\big(\begin{smallmatrix}1&0\\11&1\end{smallmatrix}\big)\in\Gamma_1(11)$）。
故 $x\circ\pi,y\circ\pi$ 的尖点阶为
$$\operatorname{ord}_{k/11}(x\circ\pi)=(-1,+1,+1,0,-1),\qquad
\operatorname{ord}_{k/11}(y\circ\pi)=(-1,+3,0,0,-2)$$
（$k=1,\dots,5$；Abel 积分 60 位计算，与本小节开头的精确除子交叉吻合）。
**(2) Siegel 表示**：记 $G_a:=\prod_{b\bmod 11}g_{a,b}$；用 Kubert–Lang 尖点阶
$\operatorname{ord}_{(r,t)}g_{a,b}=\tfrac{11}{2}B_2(\{(ar+bt)/11\})$（$X(11)$ 的 60 个
尖点上），精确有理线性代数给出
$$x\circ\pi=-\frac{G_4G_5}{G_2^2},\qquad
y\circ\pi=\frac{G_1G_5^{\,3}}{G_2^{\,3}G_3}.\qquad\text{(P)}$$
两边均 $\Gamma_1(11)$ 不变且尖点除子相同，故各比值为常数；常数为单位根（在四个
$\tau$ 点算到 70 位：$-1$ 与 $+1$），对 $\eta$ 不可见：$\log|C|=0$，且对这些行均匀
乘积 $d\arg$ 修正恒为零。
**(3) 闭链**：取模符号 $\gamma^-=\{0,\tfrac{3}{11}\}-\{0,\tfrac{8}{11}\}$：它在
$X_1(11)$ 上闭合（$3\equiv-8\bmod 11$）、在 $\tau\mapsto-\bar\tau$ 下反不变、且本原
（其周期先验为 $w_{\mathrm{anti}}$ 的非零整数倍，60 位等于 $w_{\mathrm{anti}}$），故即
反不变生成元引理的类 $\pm(a-2b)$。连分数把它分解为七个 Manin 符号：
$$\begin{aligned}
\{0,\tfrac{3}{11}\}&=\textstyle
 +\big[\begin{smallmatrix}1&0\\3&1\end{smallmatrix}\big]
 -\big[\begin{smallmatrix}1&1\\3&4\end{smallmatrix}\big]
 +\big[\begin{smallmatrix}3&1\\11&4\end{smallmatrix}\big],\\
\{0,\tfrac{8}{11}\}&=\textstyle
 +\big[\begin{smallmatrix}1&0\\1&1\end{smallmatrix}\big]
 -\big[\begin{smallmatrix}1&2\\1&3\end{smallmatrix}\big]
 +\big[\begin{smallmatrix}3&2\\4&3\end{smallmatrix}\big]
 -\big[\begin{smallmatrix}3&8\\4&11\end{smallmatrix}\big].
\end{aligned}$$
**(4) 求值**：沿这些符号对 (P) 逐项施用 (S)，分圆系数域上的精确 $q$ 展开算术
（251 个系数，远超 $M_2(\Gamma_0(11))$ 的 Sturm 界 2）证明总权 2 形式恰为
$$F_{\mathrm{total}}=-2f_{11}\ \in\ M_2(\Gamma_0(11)):$$
两个符号（$\big[\begin{smallmatrix}1&0\\3&1\end{smallmatrix}\big]$、
$\big[\begin{smallmatrix}1&2\\1&3\end{smallmatrix}\big]$）各贡献 $-f_{11}$，其余五个
恒为零。因此 $\Lambda$ 和等于 $-2\Lambda(f_{11},0)=-2b_{11}$——由
$\Lambda(f,s)=11^{s/2}(2\pi)^{-s}\Gamma(s)L(f,s)$ 与 $L(f_{11},0)=0$（根数 $+1$）有
$\Lambda(f_{11},0)=L'(f_{11},0)=b_{11}$。乘上 (S) 的因子 $\pi$，按 (3) 的定向得
$\int_{\gamma^-}\eta(x,y)=-2\pi b_{11}$，即一般情形的 (A)。$\blacksquare$
独立佐证：数值 $\Lambda$ 和 50 位等于 $-2b_{11}$（`siegel_anchor_step8.py`）；沿
$\gamma^-$ 反定向直接数值积分 $\eta(x,y)$ 得 $+0.9559686854216584787\ldots$（45 位，
`siegel_anchor_step11.py`）。

(C3) 的 regulator 一侧由此成为**已证定理**，不依赖 Bloch 菱形公式；余下的缺口——Boyd
劈裂积分链与闭链 $C'=2\gamma^-$ 的等同——由 §9.2 的认证闭合。旧稿"K₂ 秩 1 +
认证积分钉 $\lambda=1$"的循环论证已删除。

**Bloch 菱形公式的地位（不使用）**：为明确文献状态，记录本计算对 Bloch 定理的含义。
Lalín–Ramamonjisoa Thm. 6 转述的该定理称
$$\int_{\gamma}\eta(f,g)=D_E\big((f)\diamond(g)\big)\qquad\text{(B)}$$
对反不变**子群** $H_1(E,\mathbb Z)^-$ 的生成元 $\gamma$ 成立，$D_E$ 由 L–R Def. 5
eq. (10) 给出（该级数与我们实现的级数逐项相同；两条曲线上 60 位的吻合是**非严格的
数值核对**）。在本曲线（$\Delta<0$）上，(B) 的子群 factor-1 读法与锚定的已证值
不相容：$D_E((x)\diamond(y))=\frac52D_E(P)=\pi b_{11}$ 而
$\int_{\gamma^-}\eta(x,y)=\pm2\pi b_{11}$。相反，**一切**已核实的数据——上面的锚定
定理、Brunault–Bertin 对 $\{x_W,y_W\}$ 的值（字典交叉验证段）、两个生成元上的直接
数值积分（"因子 2——已解决"段）、以及 L–R 已证的 conductor-17 计算（那里子群与余商
重合，指数引理 1）——都与 (B) 以因子 1 对**余商** $H_1/H_1^+$ 的生成元 $\bar\gamma$
成立相容（余商正是 regulator 天然定义于其上的群；两端对 symbol $\mathbb Q$-线性
延拓）；由指数引理 2，子群生成元的读数则恰为两倍，与观测一致。我们记录此为
$\Delta<0$ 曲线上 Bloch 定理的显见正确表述——L–R Thm. 6 转述的子群陈述在该情形
不精确（他们的论文从未用到该情形）——但再次强调：**本文没有任何地方依赖它**。exotic relation $D_E(2P)=\frac32D_E(P)$
现同时引 Bertin 两篇：CRM Proc. Lecture Notes **36** (2004) 与
J. Reine Angew. Math. **569** (2004) 175–188（后者即 Brunault 的参考文献 [10]）；
证明出处为 CRM 版（**更正**维持：Crelle 版中此关系仍是猜想；zbMATH 书评与
Touafek–Kerada、Mellit 三处佐证 CRM 版才是证明出处）。

**归一化问题（第十波定论；第十二波由指数引理解释）**：conductor-17 的证明中 Bloch 定理用
Lalín–Ramamonjisoa Thm. 6 的归一化：$r(\{x,y\})[\gamma]=D^E((x)\diamond(y))$、
$r[\gamma]=\int_\gamma\eta$。该因子-1 陈述**在该曲线上**的正确性由自洽性论证保证
（我们自己的重构；L–R 原文中该常数 $C$ 从未被确定）：L–R §7 中实闭链的类是 $\gamma$
的先验未知整数倍 $C$，将因子-1 定理用于其 $(X)\diamond(Y)$ 并与其已证的 Corollary 2
（其 Thm. 1 eq. (5)，结合 Zudilin 恒等式即其 eq. (6)）比较，迫使 $C\cdot f=1$ 且
$C\in\mathbb Z\setminus\{0\}$，故 $f=1$；因子-2 陈述则迫使 $C=1/2\notin\mathbb Z$。
（这一自洽性论证现由指数引理 1 解释：在 $\Delta>0$ 的 conductor-17 曲线上反不变子群
**与余商重合**，故因子-1 定理不存在任何子群读法的歧义。）其 §5 的底层 $D_E$ 恒等式
我们已 60 位复现。k=0 的 $\diamond$-形式经**同一个余商生成元读法**进入：因子-1 陈述在
conductor-11 曲线上同样成立——对余商生成元 $[\gamma^-_0]=[b]$ 有
$\int_{[\gamma^-_0]}\eta=D_E((\cdot)\diamond(\cdot))$（数值验证见下段），而子群生成元
$\gamma^-=a-2b=-2[b]$ 的读数恰为其两倍（指数引理 2）——这正解释了全部三个**不同**
tempered symbol 观测到的统一因子 2（见下段）。这不影响任一证明：conductor-17 论证在
自己的曲线上自洽，conductor-11 的证明自 rev4 起完全不使用 $\diamond$-形式（regulator
由 Siegel 单位定理直接计算，见上锚定定理）；本段仅为文献状态的澄清。
旧笔记"$\pi r=D_E(\diamond)$（因子 2，Touafek 转述）"与两组数据都不符，已弃用。
**因子 2——已解决（指数现象；数据：`code/bertin_diamond.py` + `bertin_diamond.gp`，存档
`notes/attack13-bertin-diamond.txt`；新验证 `code/verify_coinvariant.gp`，存档
`notes/attack15-coinvariant.txt`）**：形式展开
$D_E((x)\diamond(y))=\frac52D_E(P)=\pi b_{11}$，而沿**子群**生成元的认证积分
$\int_{\gamma^-}\eta=2\pi b_{11}$——因子-1 归一化（L–R Thm. 6）下两边恰差 2。我们用
同一曲线上三个**独立** tempered symbol 考察此事：Weierstrass symbol、我们的 $S_0$
symbol、以及 Bertin (C1) 三次 $(X+1)(Y+1)(X+Y+1)+XY=0$ 的坐标 symbol $\{X,Y\}$（经显式
Riemann–Roch 变换映到 $E$，PARI `ellidentify` 认证像为 `11.a3`，精确变量代换
$[1,-1,-2,2]$）。$\diamond$ 值分别为 $-\frac52$、$+\frac52$、$+\frac{35}{2}$ 倍
$D_E(P)$，而沿 $\gamma^-$ 直接积分 $\eta$ 给出 $-5$、$+5$、$+35$ 倍 $D_E(P)$（前两个
如上认证；第三个由外推 Riemann 和，与 $2\pi m(C_1)=14\pi b_{11}$ 一致）：**因子在三个
情形都恰为 2**。故它是 $(E,\gamma^-)$ 的曲线级性质，与 symbol、$\diamond$ 约定（两曲线上
相同）、tame symbol（皆单位根）、$D_E$ 级数无关。

**解决**：这就是指数引理的指数现象。L–R Def. 3 的 regulator 配对在
$H_1(E(\mathbb C),\mathbb Z)^+$ 上为零（其 Remark 4），故经**余商** $H_1/H_1^+$ 分解；
Bloch 定理 (B) 因此以因子 1 评估**余商**生成元 $[\gamma^-_0]=[b]$ 上的积分，而**子群**
生成元 $\gamma^-=a-2b$ 上的积分获得 $\mathbb Z\gamma^-$ 在 $H_1/H_1^+$ 中的指数 2
（指数引理 2）。conductor-17 曲线上指数为 1（指数引理 1：$\Delta>0$、$b\mapsto-b$），
子群与余商重合，同一计算以因子 1 闭合——那里有 L–R 已证定理佐证。三个 symbol 的统一
因子 2 因此正是 $\Delta<0$ 指数，conductor-17 的因子 1 是其 $\Delta>0$ 对照。
我们还在 conductor-11 曲线本身数值验证了余商读法（`code/verify_coinvariant.gp`，存档
`notes/attack15-coinvariant.txt`）：沿平移以避开 $\eta(x_W,y_W)$ 极点的闭链积分，采用
**中心格式**离散（$\log|\cdot|$ 取子区间中点采样、乘以精确辐角增量；经验实测收敛阶
在 $N=500$ 到 $N=8000$ 的每次翻倍均为 $p=2.000000$），
$$\textstyle\int_a\eta=0,\qquad
\int_b\eta=-\pi b_{11}=D^E\big((x_W)\diamond(y_W)\big),\qquad
\int_{a-2b}\eta=2\pi b_{11},$$
其中第一个积分**谱收敛**到 $0$（$N\ge2000$ 起低于 $10^{-75}$ 的工作精度下限——配对
在 $H_1^+$ 上为零的直接体现）；第二、三个的原始误差如 $N^{-2}$ 衰减（$N=500$ 时
$-1.4\times10^{-5}$、$-6.5\times10^{-6}$，$N=8000$ 时 $-5.4\times10^{-8}$、
$-2.5\times10^{-8}$），最细一对（$N=4000/8000$）上做单步 $p=2$ Richardson 外推后，
残差为 $-1.8\times10^{-15}$（对 $-\pi b_{11}$）与 $-9.2\times10^{-16}$（对
$2\pi b_{11}$），且 $-2\int_b\eta-\int_{a-2b}\eta=4.6\times10^{-15}$（外推值上）。
**我们强调这只是数值一致性核查，不是证明**。在此前提下，因子-1 恒等式 (B) 对
$\{x_W,y_W\}$ 在余商生成元上成立，其子群生成元读数即加倍的 (A)。副产品：$D_E((X)\diamond(Y))=7\,D_E((x)\diamond(y))$ 精确成立，在 $K_2(E)$
层面解释了 Bertin $m(C_1)=7b_{11}$ 中的系数 7。Touafek 转述的形式
$\pi r=D_E(\diamond)$（其 Thm. 1）与因子-1、因子-2 两组数据都不相容，本文不用。

**与 Brunault symbol 字典的交叉验证（核查，非证明的一部分）**：我们的证明只经 conductor-17
已证实例（附录 A）使用 Bloch 定理；本段为交叉验证，不属于证明。Brunault
Thm. 3.9.3.118 对 $X_1(11)$ 上的 Weierstrass 坐标 $x,y$ 以 $D_E$ 评估
$r_{\gamma^-}\{x,y\}$，其证明终于 (3.210)–(3.211)：
$$r_{\gamma^-}\{x_W,y_W\}=\frac1{2\pi}D^E\big(8(O)+5(A)-5(2A)\big)
=-\frac{5}{2\pi}D_E(P),\qquad\text{(3.210)}$$
$$\big|r_{\gamma^-}\{x_W,y_W\}\big|=\frac{5}{2\pi}D_E(P)=b_{11},\qquad\text{(3.211)}$$
第二式把 Brunault 自己的 Cor. 3.5.101 $\zeta$ 值 $L'(E,0)=5D_E(P)$（(3.151)）代入其
(3.211)。该坐标对与我们的 $\{x_W,y_W\}$ 的识别有三重独立核查：*文本*（证明中写明
$x,y$ 是 "coordonnées de Weierstrass" 经 $j^*$ 在 "sur un modèle de Weierstrass" 上的
拉回）；*除子*（$\operatorname{div}x_W=[A]+[4A]-2[O]$、
$\operatorname{div}y_W=2[A]+[3A]-3[O]$，故 $(x_W)\diamond(y_W)\equiv8(O)+5(A)-5(2A)$
mod $\mathbb Z[E(\mathbb Q)_{\mathrm{tors}}]$）；*数值*（$5D_E(P)/(2\pi)=b_{11}$ 与 §9.2
认证值 60 位全符）。由于 Brunault 的 $\gamma$ 是**子群**生成元（其脚注 2：
$H_1(E(\mathbb C),\mathbb Z)^-=\mathbb Z\gamma$），(3.211) 读作
$\big|\int_{\gamma^-}\eta(x_W,y_W)\big|=2\pi b_{11}$——正是指数引理 2 作用于余商值
$\int_{[\gamma^-_0]}\eta(x_W,y_W)=D^E((x_W)\diamond(y_W))=-\pi b_{11}$ 的结果：计入指数
后，字典、Bloch 锚定与我们的数值三者一致。

Samart 2023 仍将 (C3) 记录为开放恒等式，缺口不在 regulator 计算，而在于
**Boyd 的劈裂积分链与 $H_1(E,\mathbb Z)^-$ 闭链的等同**——这一环由第五波的
$C'=\tilde\gamma+\beta_0=2\gamma^-$ 认证闭合（§9.2）。

**结论**：(C3) 的全部成分均已就位，且每一环要么是已证定理
（tempered、modular units、Brunault Siegel Thm. 1（regulator 锚定）），
要么已被高精度数值 + 精确代数双重锁定（闭链 $C'=2\gamma^-$、除子、$D_E$ 值）。
(C3) 由此从"开放数值猜想"降级为"已证定理的组装 + 书写级工作"。

### 9.2 闭链引理的严格化与 (C3) 的证明（第五波）

完整细节见 `notes/proof-n1.md`；数值认证 `code/n1_certify.py`（输出 `notes/attack9-n1.txt`）。

*约定（点值提升）*：每个环面分支 $\theta\mapsto y(e^{i\theta})$ 都默认视为点值提升
$\theta\mapsto(e^{i\theta},y(e^{i\theta}))\in E(\mathbb C)$；因此下文
"$y_{\mathrm{big}}(c^-)=-P$" 一类等式意指提升后的点在 $\theta=c^-$ 处等于 $E$ 上的点
$-P=(i,-e^{i\pi/4})$，而链的边界是 $E$ 上的除子。

**构造（精确代数）**。$\tilde\gamma$ 的边界为
$\partial\tilde\gamma=[P]-[\bar P]+[-P]-[-\bar P]=:D$，$P=(i,e^{i\pi/4})$
（精确：$x=i$ 时 $S_0=y^2-i$；分支跳跃值由数值读取并与精确值对照）。注意 $P$ **非扭**
（第一波 `endpoint_torsion2.py` 精确群律 20 倍无周期）——但闭性不需要端点扭：
取小分支补偿弧 $\beta_0=\alpha_1+\alpha_2$（$\alpha_2$：内弧 $P\to\bar P$；
$\alpha_1$：外弧 $-P\to-\bar P$，内区间连接不了这两点是关键修正），则
$\partial\beta_0=-D$、$c(\beta_0)=-\beta_0$，于是
$$C'=\tilde\gamma+\beta_0:\qquad \partial C'=0,\quad c(C')=-C',$$
是闭的、反不变的整系数闭链。

**同调类（认证）**。属 1、$\Delta<0$ 时周期格有基 $(w_1,w_2)$：$w_1$ 实、
$\Re(w_2/w_1)=1/2$（单实分支）；对 11.a3，PARI `E.omega` 给出
$\tau=w_2/w_1=1/2+0.2298780212\ldots i$（上半平面约定）。反不变同调此时是显式的：

**引理（反不变生成元）**：设 $E/\mathbb R$ 有 $\Delta<0$、周期基 $(w_1,w_2)$ 如上，
$(a,b)$ 为 $H_1(E,\mathbb Z)$ 的对偶基。复共轭的作用为
$$a\mapsto a,\qquad b\mapsto a-b.$$
于是 $H_1(E,\mathbb Z)^-=\mathbb Z\,(a-2b)$，生成元 $a-2b$ **本原**，其周期
$\mathrm{period}(a-2b)=w_1-2w_2=-2i\,\mathrm{Im}\,w_2$ 是周期格唯一（差符号）的纯虚周期。

*证明*：共轭固定实周期 $w_1$ 故固定类 $a$，把 $b$ 送到周期为 $\overline{w_2}=w_1-w_2$
的闭链即类 $a-b$。类 $\alpha a+\beta b$ 反不变 ⟺ $(\alpha+\beta)a-\beta b=-\alpha a-\beta b$，
即 $\beta=-2\alpha$，故 $H_1(E,\mathbb Z)^-=\mathbb Z(a-2b)$；$a-2b$ 本原因 $\gcd(1,2)=1$。
其周期 $w_1-2w_2=-2i\,\mathrm{Im}\,w_2$ 纯虚；反之纯虚周期 $mw_1+nw_2$ 须满足 $2m+n=0$，
故为 $w_1-2w_2$ 的整数倍。$\blacksquare$

周期配对 $\gamma\mapsto\mathrm{period}(\gamma)$ 把 $H_1(E,\mathbb Z)$ 等同于周期格，
故与 $\omega$ 的配对单射，$\mathrm{period}(\gamma^-)=\pm w_{\mathrm{anti}}$，于是
$\mathrm{period}(C')/w_{\mathrm{anti}}$ **先验为整数**。实测（60/80 位双精度对照）
$$\mathrm{period}(C')=I_{\mathrm{signed}}+A_{s,\mathrm{outer}}-A_{s,\mathrm{inner}},
\qquad \frac{\mathrm{period}(C')}{w_{\mathrm{anti}}}=1.9999999999999999\ldots,$$
$|\mathrm{period}(C')-2w_{\mathrm{anti}}|=2.5\times10^{-16}\ll 1$ ⟹
$\mathrm{class}(C')=2\gamma^-$（第六波 Arb 铁证：比值球含 $-2$、半径 $<1/2$，
符号=定向约定）。**更正**：开链 $\tilde\gamma$ 本身不闭，不能赋同调类；
此处所用的闭链 $C'$ 绕数为 **2**。

**同调不变性引理（与同调的配对）**：$\int_{C'}\eta(x,y)$ 只依赖于 $C'$ 在 $H_1(E,\mathbb Z)^-$ 中的同调类。
事实上，三次模型 $y^2+y=x^3-x^2$ 的面多项式全分圆，由 Rodríguez Villegas 的 tempered 判据，
symbol $\{x,y\}$ 是 tempered 的：其全部 tame 符号 $T_v(\{x,y\})$ 皆为单位根
（`code/bertin_diamond.py` 精确验证，cf. §7）。于是 $\eta$ 在支撑每点 $v$ 处的**实**留数为
$\pm\log|T_v(\{x,y\})|=0$：$\eta$ 不带 delta 质量的留数，在支撑补集上为闭形式，在支撑点处
仅有（可积的）对数奇点。由此与 $\eta$ 的配对下降到同调，再经反不变性下降到
$H_1(E,\mathbb Z)^-$（反不变生成元引理）。此处所用的链是 admissible 的：$C'$ 避开
$\operatorname{supp}\operatorname{div}(x)\cup\operatorname{supp}\operatorname{div}(y)$
（那些点 $x\in\{0,\infty\}$，而 $C'\subset\{|x|=1\}$）；折点处 $\log|y|=0$，故
$\eta=-\log|y|\,d\theta$ 沿 $C'$ 连续可积。

**regulator 合成（精确积分代数）**。$|x|=1$ 上 $\log|y_{\mathrm{big}}|=-\log|y_{\mathrm{small}}|$
逐点成立（两根之积 $=x^3$），直接计算得
$$\int_{\beta_0}\eta=2(J_1-J_2)=\int_{\tilde\gamma}\eta
\qquad\Longrightarrow\qquad \int_{C'}\eta=2\int_{\tilde\gamma}\eta,$$
其中 $J_1=\int_0^{\pi/2}\log|y_s|d\theta$、$J_2=\int_{\pi/2}^{\pi}\log|y_s|d\theta$。
再由 §9.1 的锚定（regulator 锚定定理）：$\{x,y\}$ 是 $X_1(11)$ 上的 modular units，
Brunault 的 Siegel 单位 regulator 定理（JNT 163 (2016) Thm. 1，已证）经四步精确计算
直接给出 $\int_{\gamma^-}\eta(x,y)=\pm2\pi b_{11}$（(A)；Siegel 表示 (P) + 七个 Manin
符号 + $F_{\mathrm{total}}=-2f_{11}$ 精确恒等式）。合成：
$$\int_{\tilde\gamma}\eta=\frac12\int_{C'}\eta=\pm 2\pi b_{11},$$
故由结构恒等式 $I_{\mathrm{split}}=\frac{1}{2\pi}\int_{\tilde\gamma}\eta$（§6.3，已证）得
$|I_{\mathrm{split}}|=b_{11}$。两个符号事实完成剩下的选择。**其一**，
$b_{11}=\frac{11}{4\pi^2}L(E,2)>0$：$L(E,s)$ 的 Euler 乘积在 $s=2$ 绝对收敛，且每个因子
$(1-a_p p^{-2}+p^{-1})^{-1}$ 为正，因由 Weil 界 $|a_p|\le2\sqrt p$ 有
$1-a_p p^{-2}+p^{-1}\ge1-2p^{-3/2}+p^{-1}>0$。**其二**，认证区间包围（球算术、自适应二分
加严格值域界；`code/sign_certify.py`，证书 `notes/attack14-sign-k0.txt`）给出
$$I_{\mathrm{split}}\in[\,0.1489,\ 0.1553\,]\subset(0,\infty).$$
故 $I_{\mathrm{split}}=+b_{11}$：
$$\boxed{\,I_{\mathrm{split}}=b_{11}\,}\qquad\text{(C3) 证毕}.$$

**区间算术铁证化（第六波完成）**：整数识别已由 Arb 球算术完全严格化
（`code/n1_interval.py`，输出 `notes/attack10-interval.txt`）：三段积分用 Arb 认证
积分重算（θ=±t² 换元消端点奇性 + Cauchy 尖端估计——尖端常数 $a_0$ 由局部展开手算得到
但不予轻信：脚本在 $t=\delta$ 用球算术重新评估 $f$，四个尖端逐一认证
$|f(\delta)-a_0|\le H\delta^2$，$a_0$ 若有符号或分支错误会使值移动 $2$，远超此界；
D 穿负实轴处自适应细分 +
每段认证回避割线 + 认证符号传递；w_anti 由 Newton+Rouché 隔离的根 + Carlson RF
独立认证，与 PARI 45 位一致），得比值球含 $-2$ 且半径 $<1/2$，先验整数性 ⟹
$\mathrm{period}(C')/w_{\mathrm{anti}}=-2$ 为严格等式（符号=定向约定）。
模型常数 $\kappa=1$ 由判别式比 $\kappa^{12}=\Delta_{\mathrm{quartic}}/\Delta_{\min}=2^a11^b$
的离散候选 + 50 位一致锁定。

**两项实现注记**：(i) 存档输出（如 `notes/attack10-interval.txt`）中显示的球半径遵循
Arb/python-flint 的 repr 约定——中点打印误差被折入显示半径；真实球半径更小
（显示区间包含真实球）。上文引用的所有认证界均用精确半径，而非显示半径。
(ii) 弧积分代码中，割线回避可能在某个子弧的首段选取 $i\sqrt{-D}$ 变体；现行脚本对
**每个**子弧（含首段）都认证叶选择。存档的 k=0 输出早于该逻辑的一次重构、已按其复审：
k=0 的割线穿越点位于子弧内部，其首段总用 $\sqrt D$ 变体，故认证值不受影响（取错叶会使
整段弧反号、比值以 $O(1)$ 幅度偏离 $-2$，存档输出排除了这一情形）。

**定理（认证计算）**：以下各条由精确有理算术或认证区间（球）算术证明，全部证书存档于
`notes/`、全部脚本在 `code/`（复现命令见 §11；软件：Python 3.12、mpmath、sympy、
python-flint 0.9.0/Arb、PARI/GP 2.15.5）：

1. **精确代数输入**：从 $S_0=0$ 到极小三次 $E:y^2+y=x^3-x^2$ 的双有理映射与
   $\omega_{\min}=dx/u$（$\kappa=1$）是精确多项式恒等式（`kappa_exact.py`）；$x,y$ 的除子、
   金刚石积、tame 符号、群律与扭子群 $E(\mathbb Q)=\mathbb Z/5\mathbb Z$ 在
   $\mathbb Q,\mathbb Q(\sqrt2),\mathbb Q(\zeta_8)$ 上精确（`torsion.py`、`bertin_diamond.py`）；
   链端点 $P$ 的非扭性由好素数 $17$ 与 $89$ 处的既约证明（`endpoint_torsion3.py`）。
2. **弧、定向、闭性**：§9.2 的链 $\tilde\gamma,\beta_0,C'$ 是 $|x|=1$ 上定向固定的显式
   参数化弧；端点分支指派（故 $\partial C'=0$、$c(C')=-C'$）在两种尺度
   $\varepsilon=10^{-6},10^{-9}$ 下经球算术认证（`branch_certify.py`，证书
   `notes/attack12-branch.txt`）。
3. **每段上的分支延拓**：周期积分的每个子弧上，所选平方根分支（$\sqrt D$ 或 $i\sqrt{-D}$）
   的解析性由整个子弧上的球求值认证——不限于端点附近——整体符号由认证的匹配检验跨节点
   传递（`n1_interval.py`）；结构命题的模序 $|y_{\mathrm{small}}|<1<|y_{\mathrm{big}}|$
   由自适应二分认证（`branch_certify.py`）。
4. **区间包围**：每段弧积分与本原周期 $w_{\mathrm{anti}}$ 由 Arb 认证自适应积分与
   Newton–Rouché 根隔离的 Carlson $R_F$ 包围，积分半径 $\le10^{-3}$、周期半径 $\le10^{-48}$
   （`n1_interval.py`，证书 `notes/attack10-interval.txt`）。
5. **整数识别**：类 $[C']\in H_1(E,\mathbb Z)^-=\mathbb Z\gamma^-$ 先验为整数（反不变生成元
   引理与同调类段）；认证比值球满足 $|\mathrm{ratio}-2|\le4.33\times10^{-3}<1/2$，最近整数
   唯一，故 $\mathrm{class}(C')=2\gamma^-$。
6. **符号**：最终符号由 $I_{\mathrm{split}}$ 的认证包围固定：
   $I_{\mathrm{split}}\in[0.1489,0.1553]\subset(0,\infty)$（`sign_certify.py`，证书
   `notes/attack14-sign-k0.txt`；所有值域包围——含未分离段上的凸包——完全在球算术内
   构造并验证，带程序化包含断言），而 $b_{11}=\Lambda(f_{11},2)>0$ 由 $s=2$ 处绝对收敛的
   Euler 乘积精确成立。
7. **Siegel 锚定**：Siegel 表示 (P) 是 $\Gamma_1(11)$ 不变函数的精确恒等式（Kubert–Lang
   阶的精确尖点除子匹配；常数为单位根，算到 70 位）；$\gamma^-$ 的 Manin 符号分解精确
   （连分数）；权 2 形式恒等式 $F_{\mathrm{total}}=-2f_{11}$ 由精确 $q$ 展开算术证到
   251 个系数、远超 Sturm 界 2（`siegel_anchor_step6/7/9`，存档
   `notes/attack16-siegel-anchor.txt`）；尖点-扭点对应与闭链本原性对其先验整数值
   60 位确认。

主定理证明中唯一的外部输入是 Brunault 的 Siegel 单位 regulator 公式（J. Number
Theory **163** (2016) 542–569，Thm. 1——原文以 Rankin–Selberg 计算证明），连同函数
方程求值 $\Lambda(f_{11},0)=L'(f_{11},0)=b_{11}$。Bloch 菱形定理、Bertin Thm. 6 与
Brunault 的 symbol 字典**不**用于逻辑链（后者仅作交叉验证保留，§9.1 末"与 Brunault
symbol 字典的交叉验证"段）。

## 10. 总结

1. (C1)(C2) 独立复现至 52 位；(C3) 确认至 **366 位**（$|I_{\mathrm{split}}-b_{11}|=9.26\times10^{-367}$，PARI `lfun` 330 位交叉吻合；原公开记录 50 位，Boyd 2015 slides p. 28）。
2. 新结构定理 $|y_-|\le1\Rightarrow I_1+I_2=-m(S_0)$，把 (C3) 化为 $I_1=(b_{11}-m(S_0))/2$。
3. $m(S_0)$ 对初等常数 PSLQ 阴性（界 $10^{10}$）。
4. modular units 前提**成立**（$5A=O$ 精确验证）。
5. **更正**：第一波"朴素 BMZ 被非扭边界阻断"的断言不成立——正确积分链 $\tilde\gamma$ 在 $H_1(E,\mathbb Z)^-$ 中拓扑闭合，与扭点无关（§8）。
6. **(C3) 证毕（完全严格，§9.2）**：闭链引理严格化——$\tilde\gamma$ 经小分支补偿弧闭化为 $C'=\tilde\gamma+\beta_0$（闭、反不变、整系数），比率先验整数；15 位匹配认证 $\mathrm{class}(C')=2\gamma^-$（更正第二波"绕数 $n=1$"的表述），并于第六波由 Arb 球算术铁证化（比值球含 $-2$、半径 $<1/2$）；配合 $\int_{\beta_0}\eta=\int_{\tilde\gamma}\eta$（精确积分代数）与 Brunault 的
Siegel 单位 regulator 定理（JNT 163 (2016) 542–569，Thm. 1，已证）——$\{x,y\}$ 为
modular units，四步精确计算直接给出 $\int_{\gamma^-}\eta=\pm2\pi b_{11}$（§9.1 锚定
定理；不使用 Bloch 菱形定理）——得 $\int_{\tilde\gamma}\eta=2\pi b_{11}$
为**严格等式**，$I_{\mathrm{split}}=b_{11}$ 证毕。
7. 族结果（第四波完成 + rev2/rev3 精确化）：$\tilde n(k)$ 表 + **精确结构分析 + 数值证据**框架——环面交精确刻画（命题：交 ⟺ $k\in[-4,2]$）与模单位/temperedness 结构精确刻画（模单位仅 $k=0$；$k=1$ 为 tempered + 扭点支撑，足以支撑 conductor-17 论证）；$k\notin(-4,2)$ 时 Deninger 机制预期适用、$m(S_k)=r_k|L'(E_k,0)|$ 获数值确认（$k=2,3$ 确认，$k=-4,-5,-6$ **先预测后命中**，$r=2,1,\frac72,\frac14,\frac18$）；$-4<k<2$ 时恒等式证于 $k=0$（$\mathbb Z/5$）与 $k=1$（$\mathbb Z/4$，conditional on 归一化引理），$k=-1,-2,-3$ PSLQ 阴性（界 $10^8$ 内未发现）。**conductor 53 机制失效**：闭链空间精确一维、生成元周期由 $1/u_++1/u_-\equiv0$ 恒等于 0（环面上无反不变闭链可实现）+ $x,y$ 非 modular units（精确除子论证：53.a1 为 strong Weil curve，$X_0(53)$ 尖点全映到 $O$，而 $\operatorname{div}(x),\operatorname{div}(y)$ 的支撑含非扭生成元）——机制不适用，在缺乏精确候选公式下**留作开放**（§8.4）。
8. **conductor 17（第七波，附录 A）**：同一方法再应用于 $k=1$——链结构逐字平行（$c=2\pi/3$、跳跃值 $y=\pm i$ 精确）、$\mathrm{class}(C')=\pm2\gamma^-$ Arb 铁证、regulator 侧由 Lalín–Ramamonjisoa 已发表定理闭环——Samart 的 $\tilde n(1)=b_{17}$ 类比猜想成为定理，**conditional on 附录 A 讨论的归一化引理**。
9. **rev2/rev3（第十一、十二波）严格化**：环面交/模单位的精确刻画（§8.4）+ 反不变生成元本原性引理（§9.2）+ 符号的认证区间包围（`sign_certify.py`、`k1_sign_certify.py`）+ 认证计算定理（§9.2 末六项清单）+（第十二波）因子 2 的源头定位——反不变子群在余商 $H_1/H_1^+$ 中的指数（$\Delta<0$ 时指数 2；指数引理，§9.1），Bloch 定理重述为余商生成元的因子-1 恒等式 (B)，直接锚定 (A) 取代已删除的比例引理，并新增 `verify_coinvariant.gp` 直接数值验证（$\int_b\eta=-\pi b_{11}$ 因子 1、$\int_{a-2b}\eta=2\pi b_{11}$）；Brunault Thm. 118 字典降为交叉验证；k=1 材料移入附录 A 并标注 conditional；LICENSE 与 requirements.txt 落定（§11）。
10. **rev4（第四轮审稿后）**：regulator 锚定再升级——证明不再使用 Bloch 菱形定理、Bertin Thm. 6 与 symbol 字典，改用 Brunault 已证的 Siegel 单位 regulator 公式（JNT 163 (2016) 542–569，Thm. 1）直接严格计算 $\int_{\gamma^-}\eta=-2\pi b_{11}$（§9.1 锚定定理：尖点↔扭点对应 $k/11\mapsto m_kA$、Siegel 表示 (P)、七个 Manin 符号、$F_{\mathrm{total}}=-2f_{11}$ 由 251 系数精确 $q$ 展开证明；存档 `notes/attack16-siegel-anchor.txt`，脚本 `siegel_anchor_step1`–`11`）；认证计算定理新增第 7 项记录其机器可验证部分，唯一外部输入 = Brunault Thm. 1 + 函数方程 $\Lambda(f_{11},0)=b_{11}$；`verify_coinvariant.gp` 改中心格式（实测收敛阶 $p=2$），明确标注为数值一致性核查。


## 附录 A. k=1（conductor 17）：同一方法的再应用（第七波；conditional）

**本附录的地位**：正文在逻辑上独立于本附录。这里证明的定理（Samart 的 conductor-17 类比）建立在恰好两条假设之上：

- **temperedness**：$S_1$ 的 Newton 面多项式全部是分圆多项式，因此
  $\{x,y\}\in K_2(E_1)\otimes\mathbb{Q}$——这一定义性性质在本附录中**被精确证明**（K_2 一节）。
  我们**不**使用、也**不**声称 modular-units 性质：在自然模模型 $X_0(17)$（仅两个尖点）上，
  支撑集含全部四个有理扭点的除子不可能全部由尖点支撑（§8.4 的命题）。
- **因子-1 归一化下的 Bloch 定理**（Lalín–Ramamonjisoa Thm. 6 的表述）。"这是该曲线上的正确归一化"
  来自下方归一化备注中的自洽性论证（我们自己的重构；L–R 原文并未确定该常数）；等价地，
  由指数引理 1 的指数-1 陈述：在此 $\Delta>0$ 曲线上反不变子群与余商重合，故因子-1 定理
  不存在任何子群读法的歧义。

因此本结果应读作 **conditional on 该归一化引理**；附录中其余一切与正文同一标准证明。

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
- **合成**：$\int_{\tilde\gamma}\eta=\frac12\cdot2\cdot(\pm2\pi b_{17})$，由结构恒等式得
  $|\tilde n(1)|=b_{17}$；符号与主定理同样敲定：$b_{17}>0$ 由同一 Euler 乘积论证，
  认证区间包围（`code/k1_sign_certify.py`，证书 `notes/attack14-sign-k1.txt`）给出
  $\frac1{2\pi}\int_{\tilde\gamma}\eta\in[-0.3026,-0.2961]\subset(-\infty,0)$，故
  $\boxed{\tilde n(1)=+b_{17}}$（Samart conductor-17 类比猜想，证毕）。
- **归一化备注（第八波定论；第十二波由指数引理解释）**：k=1 用的 Lalín 归一化
  （$\int_{\gamma^-}\eta=\pm D^E(\diamond)$，因子 1）经 conductor-17 自洽性论证确证为
  Bloch 定理的正确形式（我们的重构：因子-1 定理与 L–R 已证 Corollary 2 比较迫使
  $C\cdot f=1$，$C\in\mathbb Z\setminus\{0\}$，详见 §9.1；L–R 原文并未确定 $C$）——
  该论证现由指数引理 1 解释：在此 $\Delta>0$ 曲线上反不变子群与余商**重合**，因子-1
  定理不存在其他子群读法；k=0 的证明经**同一个余商生成元读法**使用 $\diamond$-形式：
  因子-1 陈述在 conductor-11 曲线上对余商生成元 $[\gamma^-_0]=[b]$ 同样成立
  （§9.1 数值验证），子群生成元读数恰为其两倍（指数引理 2），余商值经 conductor-17
  已证实例锚定于 Bloch 定理（§9.1）。


## 11. 复现方式

```
cd code && python b11.py && python attack1.py && python attack2.py \
  && python attack13_c3_300.py && python torsion.py && python endpoint_torsion2.py \
  && python boundary_torsion.py && python closedness_check.py \
  && python ntilde_family.py && python b_family.py && python winding.py \
  && python dilog.py && python k53_attack.py && python k53_smith.py && python kneg_m.py \
  && python n1_certify.py && python bertin_diamond.py \
  && python sign_certify.py && python k1_sign_certify.py
python n1_interval.py    # Arb 区间算术铁证（需 python-flint）
python branch_certify.py  # 分支指派 + 模序认证
python k1_interval.py     # Arb 铁证，conductor 17
gp -q verify_family.gp && gp -q verify_ratios.gp
gp -q winding.gp && gp -q dilog.gp && gp -q bertin_diamond.gp
gp -q k53.gp && gp -q k53b.gp && gp -q kfamily_torsion.gp
gp -q k1_pari.gp && gp -q k1_points.gp && gp -q k1_zvals.gp
gp -q verify_coinvariant.gp  # 余商 vs 子群生成元积分（∫_b=-πb11 因子 1；∫_{a-2b}=2πb11）
python siegel_anchor_step11.py  # Siegel 锚定：最终值（§9.1 锚定定理）
# 完整锚定链重建：siegel_anchor_step4.py -> step5.py -> step6.py
#   -> step7.gp -> step9.py -> step8.py（step8 约 10–20 分钟）
```

依赖（仓库根目录 `requirements.txt`）：Python ≥3.10 + mpmath + sympy；区间铁证另需
python-flint 0.9.0（Arb）。从干净 checkout 复现：`python -m pip install -r requirements.txt`——
项目内 `.venv` **不**是必需的、也不在归档中。

**Permanence（存档）**：完整研究仓库——全部脚本、全部原始输出存档（`notes/attack*.txt`）
与本报告源文件——以 git 版本控制，并作为电子补充材料随本报告发布；记录版本为标签
`rev3`（见仓库 log）。代码以 MIT license 发布（仓库根目录 `LICENSE` 文件）。
上文 `attackN` 形式的文件名均指该存档。

## 12. 文献导读（`literature/`）

- `bertin-lalin-survey.pdf` — Bertin–Lalín 综述：全局图景与各 conductor 状态（先读这篇）
- `boyd-pnwnt2015.pdf` — Boyd 2015 slides：猜想史 + $m(S_0)$ 原始数据；p. 28 载 (C3) 的 50 位验证（此前的公开纪录）
- `brunault-these.pdf` — Brunault 博士论文：$X_1(11)$ 上 Beilinson 定理显式化，(C1) 的证明
- `brunault-siegel.pdf` — Brunault：Siegel 单位的 regulator 及应用（J. Number Theory **163** (2016) 542–569，arXiv:1504.08127）：Thm. 1 的 Siegel 单位 regulator 公式——rev4 起主定理的锚定（§9.1 锚定定理）
- `zudilin-regulator.pdf` — Zudilin：BMZ regulator 公式（证明武器）
- `samart2023.pdf` — Samart：开放猜想 (C3) 的明确陈述（其 eq. (4.1)）+ conductor 19 的成功范例
- `lalin-samart-zudilin-cond21.pdf` — conductor 21：half-Mahler 方法范例
- `lalin-ramamonjisoa-cond17.pdf` — Lalín–Ramamonjisoa 2017：conductor 17 已证公式
  $L(E_{17},2)=\frac{8\pi}{17}D^E(P)$（k=1 证明的 regulator 闭环依据）与 Bloch Thm. 6 归一化出处
- `boyd-slides.pdf` — Boyd 关于 $L(E,3)$ 的 slides
- 详细笔记：`notes/literature-notes.md`；原始运行输出：`notes/attack*-results.txt`
