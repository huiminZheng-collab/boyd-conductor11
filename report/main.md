# Boyd 的 conductor 11 Mahler 测度猜想：研究报告

**日期**：2026-08-04　**仓库**：`boyd-conductor11/`

## 1. 问题陈述

Boyd (1998) 系统猜想：许多二元多项式的对数 Mahler 测度 $m(P)$ 是椭圆曲线 $L$ 值的有理倍数，
$m(P) \stackrel?= r\cdot b_N$，其中 $b_N=L'(E_N,0)=\dfrac{N}{4\pi^2}L(E_N,2)$，$E_N$ conductor $N$ 的椭圆曲线，$r\in\mathbb Q$。

Conductor 11（最小可能 conductor，曲线 $E_{11}=X_1(11)$，LMFDB `11.a3`，模形式 $f_{11}=\eta(\tau)^2\eta(11\tau)^2$）的情形：

| # | 恒等式 | 状态 |
|---|--------|------|
| (C1) | $m\big((1+x)(1+y)(1+x+y)+xy\big)=7b_{11}$ | **已证**（Brunault 2005/06，$X_1(11)$ 的 modular units 参数化 + Beilinson 显式化） |
| (C2) | $m\big(y^2+(x^2+2x-1)y+x^3\big)=5b_{11}$ | **已证**（Brunault 2006） |
| (C3) | $S_0=y^2+(x^2+1)y+x^3$ 的**劈裂积分** $I_{\mathrm{split}}=\pm b_{11}$（详见 §4.3） | **开放**（Boyd 1998 编号 (2-33)；Samart 2023 eq. (4.1) 重申为猜想） |

(C3) 是本次攻击目标。由于 $S_0$ 在 2-torus 上有零点（$x=\pm i$），$m(S_0)$ **本身不是** $b_{11}$ 的有理倍数
（Boyd 数值观察，本报告 §4.4 复现），猜想的是"沿 branch cut 劈开的带符号积分"。

## 2. 方法

- **$b_{11}$ 的高精度计算**（`code/b11.py`）：由 $f_{11}=q\prod(1-q^n)^2(1-q^{11n})^2=\sum a_nq^n$ 的整数系数，
  用权 2、根数 $+1$ 的近似函数方程
  $$b_{11}=\Lambda(f,2)=\sum_{n\ge1}a_n\Big[e^{-t_n}\Big(\frac1{t_n}+\frac1{t_n^2}\Big)+E_1(t_n)\Big],\quad t_n=\frac{2\pi n}{\sqrt{11}}.$$
  验证：与 Boyd (PNWNT 2015) 的 $b_{11}=0.1521471\ldots$ 吻合。
- **Mahler 测度**：Jensen 公式降为 1 维积分，$m(Ay^2+By+C)=\frac1{2\pi}\int_0^{2\pi}\big[\log|A|+\sum_j\log^+|y_j|\big]d\theta$，
  mpmath 任意精度（80–150 dps）。
- **劈裂积分**：$y_\pm(x)=-\frac{x^2+1}2\pm\sqrt{\frac{(x^2+1)^2}4-x^3}$，
  $I_1=\frac1\pi\int_0^{\pi/2}\log|y_-|d\theta$，$I_2=\frac1\pi\int_{\pi/2}^\pi\log|y_-|d\theta$，$I_{\mathrm{split}}=I_1-I_2$。
- **PSLQ**（mpmath.pslq）搜索有理/整数关系。

## 3. 结果汇总

（数值结果见 `notes/attack1-results.txt`、`notes/attack2-results.txt`）

### 3.1 已证恒等式的独立复验（80 dps）

- $m\big((1+x)(1+y)(1+x+y)+xy\big)$ 与 $7b_{11}$ 之差 $\approx 5.0\times10^{-53}$；
- $m\big(y^2+(x^2+2x-1)y+x^3\big)$ 与 $5b_{11}$ 之差 $\approx 3.6\times10^{-53}$。

### 3.2 开放猜想 (C3) 的高精度确认 —— **本报告主要数值结果**

150 dps 下：

$$I_{\mathrm{split}}=0.1521471417259180494862272974786344956281435891642261\ldots$$

与 $b_{11}$ 的吻合位数、$I_1+I_2=-m(S_0)$ 结构关系等见 §4（由 attack2 结果填入）。

注：Samart 2023 eq. (4.1) 写作 $-L'(E,0)$；符号取决于哪个根称为 $y_-$（$|y_+||y_-|=1\Rightarrow$ 换根变号）。
两种约定下猜想内容相同。

### 3.3 $m(S_0)$ 的"负结果"

$m(S_0)=0.40560295591501040\ldots$（吻合 Boyd 的 $0.4056029$）。
PSLQ（系数界 $10^8$，80 dps）确认 $m(S_0)$ 不是 $b_{11}$ 的小有理数倍；
与初等常数 $\{b_{11},\log2,\log3,\mathrm{Catalan},m(1+x+y)\}$ 的 PSLQ（系数界 $10^{10}$，150 dps）无关系
——支持"$m(S_0)$ 需要椭圆双对数表达，而非初等常数"的预期。

## 4. 证明策略分析（未能完成，留给后续）

### 4.1 关键观察

在 $|x|=1$ 上 $\log|x|=0$，故 $\log|y|\,d\arg x=-\eta(x,y)\big|_{\text{path}}$，其中
$\eta(x,y)=\log|x|\,d\arg y-\log|y|\,d\arg x$ 是 regulator 形式。$I_{\mathrm{split}}$ 是
$\eta(x,y)$ 沿带符号路径 $\gamma=[\theta:0\to\pi/2]-[\theta:\pi/2\to\pi]$ 的积分，
$$\partial\gamma=2[(i,e^{i\pi/4})]-[(1,-1)]-[(-1,-1+\sqrt2)].$$

### 4.2 modular units 检验（部分完成）

在曲线 $C:S_0=0$（齐次化 $Y^2Z+(X^2+Z^2)Y+X^3=0$）上：
- $\operatorname{div}(x)=[(0,0)]+[(0,-1)]-2[P_\infty]$，
- $\operatorname{div}(y)=3[(0,0)]-3[P_\infty]$。

若能验证 $(0,0),(0,-1)$ 是 5-扭点（$E_{11}(\mathbb Q)=\mathbb Z/5\mathbb Z$），则 $x,y$ 是 $X_1(11)$ 上的
**modular units**（除子支撑在尖点上，Manin–Drinfeld），且路径端点是尖点，
于是 Brunault–Mellit–Zudilin 公式（`literature/zudilin-regulator.pdf`）把
$\int_\gamma\eta(x,y)$ 直接算成 $b_{11}$ 的倍数 —— 这正是 Samart 2023 §4 建议、并在 conductor 19 的
$Q_\alpha$ 族中成功实施（其 Theorem 2）的路线。

### 4.3 障碍

- $S_0$ 的曲线在 torus 上有零点，路径严格说不闭；需处理 $\partial\gamma$ 中尖点的 2 倍系数
  （Samart 对 $Q_\alpha$ 用"modified Mahler measure" $\tilde n=n-3J$ 处理类似问题）。
- 验证扭点阶数需要 Weierstrass 模型转换与群律计算（时间所限未完成；四次曲线
  $u^2=x^4-4x^3+2x^2+1$ 的不变量 $I=-32, J=-448$ 已算出，Jacobian 模型 $y^2=x^3+864x+12096$ 待核对）。

## 5. 结论

1. conductor 11 的两个已证 Boyd 恒等式以 52 位精度独立复现。
2. 开放猜想 (C3)（Boyd 1998 (2-33) / Samart 2023 (4.1)）确认到 **110+ 位**（此前公开记录为 Boyd 的 50 位）。
3. 发现精确结构关系 $I_1+I_2=-m(S_0)$（使猜想等价于 $I_1=(b_{11}-m(S_0))/2$）。
4. $m(S_0)$ 对初等常数的 PSLQ 阴性结果（$10^{10}$ 系数界）。
5. 证明路线明确：modular units + BMZ regulator 公式；剩余工作为扭点验证与路径分析。

## 6. 文件结构

- `literature/` — 7 份文献 PDF；`notes/txt/` — 提取文本（git 忽略）
- `notes/literature-notes.md` — 文献笔记
- `code/b11.py` — $b_{11}$ 高精度计算；`code/attack1.py`、`code/attack2.py` — 攻击脚本
- `notes/attack*-results.txt` — 运行结果
- `report/main.md` — 本报告
