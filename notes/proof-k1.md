# k=1（conductor 17）闭链引理与 ñ(1)=b_17 的证明（第七波，2026-08-05）

## 0. 结论

Samart 的 conductor-17 类似猜想 $\tilde n(1)=b_{17}$ 沿 k=0（`notes/proof-n1.md`）的同一条
证明路线**全程打通**：链结构、闭化、反不变性全是精确代数；同调类
$\mathrm{class}(C')=\pm2\gamma^-$ 由 Arb 球算术**铁证化**（比值球含 $-2$，半径
$3.3\times10^{-14}<1/2$）；regulator 一侧的常数 $D_E(A)=\frac{17}{8\pi}L(E,2)$ 是
**已发表定理**（Lalín–Ramamonjisoa 2017，`literature/lalin-ramamonjisoa-cond17.pdf`），
比 k=0 的 Brunault (3.151) 引用更直接。数值上全部环节 60 位吻合。

唯一需要在复核时注意的一点：本证明使用 Bloch 定理的 **Lalín 归一化**
（∫_γ η = ±D^E((x)⋄(y))，无因子 2），与 Lalín–Ramamonjisoa 的 conductor-17 已证工作
自洽；k=0 笔记（proof-n1.md §4/§6）采用的"π r = D_E(⋄)、r=(1/2π)∫η"（因子 2）读法
与本情形数值**不相容**，详见 §8 的对照与存疑。

## 1. 曲线识别（PARI，`code/k1_pari.gp`、`code/k1_points.gp`）

$S_1: y^2+(x^2+x+1)y+x^3=0$。PARI `ellfromeqn` 给出 Weierstrass 模型
$E_1: y^2+xy-y=x^3-x^2$（$[1,-1,-1,0,0]$，$\Delta=17$），极小模型
$y^2+xy+y=x^3-x^2-x$（$[1,-1,1,-1,0]$，$\Delta=17$），`ellglobalred` 确认
**conductor = 17**；`elltors`：$E(\mathbb Q)_{\mathrm{tors}}=\mathbb Z/4\mathbb Z$；
$(0,0)$ 在 $E_1$ 上 `ellorder` = **4**（与 k=0 的 $(0,0)$ 5-扭完全平行）。

周期（极小模型，80 位）：
$$w_1=3.0941595071022403464191580\ldots\ (\in\mathbb R),\qquad
w_2=-2.7457391180897536720341879\ldots\,i\ (\in i\mathbb R).$$
**与 k=0 的关键差别**：$\Delta=+17>0$（$E(\mathbb R)$ 两个连通分支），
$\overline{w_2}=-w_2$ 严格成立，故原初反不变周期就是
$$w_{\mathrm{anti}}=w_2\qquad(\text{k=0 是 }\Delta<0,\ w_{\mathrm{anti}}=2w_2-w_1=2i\,\mathrm{Im}\,w_2).$$

$b_{17}=L'(E,0)=\frac{17}{4\pi^2}L(E,2)$（根数 $+1$，函数方程）：
$$L(E,2)=0.69518146300948461345173162509571731229\ldots,\quad
b_{17}=0.29935558688291539005379974769003145063\ldots$$

## 2. 闭链构造（精确代数）

torus 交点：$|x|=1$ 上 $|y_{\mathrm{big}}|$ 在 $\theta=\pm c$ 处触及 1，
$c=\arccos(-k/2)\big|_{k=1}=\frac{2\pi}{3}$（族公式，`report/main.md` §8.4）。
**精确跳跃值**：在 $\theta=c$，$x=\omega=e^{2\pi i/3}$ 满足
$B=x^2+x+1=0$，故 $S_1=y^2+x^3=y^2+1=0$，即
$$y^2=-1,\qquad y=\pm i\qquad(\text{k=0 类比：}x=\pm i,\ y^2=i).$$
记 $P=(\omega,i)$，$\bar P=(\omega^2,-i)$（复共轭），$-P=(\omega,-i)$。
分支跳跃（mpmath 60 位读取并与精确值对照，`code/k1_certify.py`）：
$y_{\mathrm{big}}(c^-)=-i=-P$、$y_{\mathrm{big}}(c^+)=+i=P$、
$y_{\mathrm{small}}(c^-)=+i$、$y_{\mathrm{small}}(c^+)=-i$——与 k=0 逐字平行。

带符号链与闭化（同 k=0 §1）：
$$\tilde\gamma=+\,[y_{\mathrm{big}}\text{ on }[-c,c]]\;-\;[y_{\mathrm{big}}\text{ on 外侧两段}],$$
$\partial\tilde\gamma=[P]-[\bar P]+[-P]-[-\bar P]=:D$（$\deg D=0$，$c(D)=-D$）；
$\beta_0=\alpha_1+\alpha_2$（小分支内弧 $\theta:c\to-c$ 连接 $P\to\bar P$ +
小分支外弧连接 $-P\to-\bar P$），$\partial\beta_0=-D$、$c(\beta_0)=-\beta_0$：
$$C'=\tilde\gamma+\beta_0:\quad\partial C'=0,\quad c(C')=-C'.$$

## 3. 同调类：$C'=\pm2\gamma^-$（mpmath + Arb 双重认证）

属 1、$H_1(E,\mathbb Z)^-=\mathbb Z\gamma^-$，周期配对单射，
$\mathrm{period}(\gamma^-)=\pm w_{\mathrm{anti}}$（$\Delta>0$ 时 $w_{\mathrm{anti}}=w_2$），
故 $\mathrm{period}(C')/w_{\mathrm{anti}}$ **先验是非零整数**。

**mpmath**（`code/k1_certify.py`，60/80 位对照，输出 `notes/attack11-k1-certify.txt`）：
$$\mathrm{period}(C')=I_{\mathrm{signed}}+A_{s,\mathrm{outer}}-A_{s,\mathrm{inner}}
=-5.4914782361795073\ldots i=2w_2,\qquad
\frac{\mathrm{period}(C')}{w_2}=2.0000000000000002\ldots$$
$|\cdot-2|=1.5\times10^{-16}$。附带 $I_{\mathrm{signed}}=w_2$（开链配对，同 k=0）。

**Arb 铁证化**（`code/k1_interval.py`，python-flint 300 位工作精度，输出
`notes/attack11-k1-interval.txt`）：全部积分用 `acb.integral` 认证；
$w_{\mathrm{anti}}$ 由 $4x^3-3x^2-2x+1=(x-1)(4x^2+x-1)$ 的 Newton+Rouché 隔离根 +
Carlson RF 独立认证（$\Delta>0$ 三实根公式 $w_{\mathrm{real}}=2R_F(0,e_1-e_2,e_1-e_3)$、
$w_{\mathrm{imag}}=2iR_F(0,e_2-e_3,e_1-e_3)$，与 PARI 45 位一致），不依赖 PARI：
$$\text{比值球 }=[-2.00\pm2.4\times10^{-14}]+[\pm2.4\times10^{-14}]\,i,\qquad
|\mathrm{ratio}+2|\le3.33\times10^{-14}<\tfrac12.$$
先验整数性 ⟹ $\mathrm{period}(C')/w_{\mathrm{anti}}=-2$ 为**严格等式**
（符号=定向约定），$\boxed{|\mathrm{class}(C')|=2\gamma^-}$。

**k=1 的 Arb 适配要点**（与 `n1_interval.py` 的差异）：
- $D(z)=z^4-2z^3+3z^2+2z+1$ 在 $|z|=1$ 上**无零点**（$\min|D|=4>0$，在角点
  取得）且 $D(1)=5\ne0$——k=0 的 $\theta=0$ 处 $1/\sqrt\theta$ 端点奇性
  **不存在**，$\theta=\pm t^2$ 换元与 Cauchy 尖端估计整套机制不需要；
- $D(\theta)$ 只在角点 $\theta=\pm c$（积分弧端点，$D=-4$）触及负实轴，
  仍用自适应细分 + 每段认证割线回避（$\sqrt D$ / $i\sqrt{-D}$）+ 认证符号传递；
- **修正了 `n1_interval.py` 的一个潜伏 bug**：弧首段的 branch_sign 是相对于主支
  $\sqrt D$ 取的 $\sigma$，但若首段用 '$r$' 变体 $i\sqrt{-D}$，需先做一次认证换算
  （k=0 时割线穿过点在弧内部、首段必为 'p'，故 k=0 未触发；k=1 的割线恰在弧端点，
  首段为 'r'，不换算了话整条弧取反）。`code/k1_interval.py` 中
  `arc_integral` 已加入首段认证换算；**建议 k=0 侧复核 `n1_interval.py`**
  （那里虽碰巧正确，但同一隐患存在）。

## 4. regulator 积分代数（精确，同 k=0 §3）

$|x|=1$ 上 $\eta(x,y)=-\log|y|\,d\theta$，两根之积 $=x^3$（模 1）
⟹ $\log|y_{\mathrm{big}}|=-\log|y_{\mathrm{small}}|$ 逐点。记
$J_1=\int_0^{2\pi/3}\log|y_s|d\theta$、$J_2=\int_{2\pi/3}^{\pi}\log|y_s|d\theta$，则
$$\int_{\tilde\gamma}\eta=2(J_1-J_2),\qquad
\int_{\beta_0}\eta=\underbrace{2J_1}_{\alpha_2}\underbrace{-2J_2}_{\alpha_1}
=\int_{\tilde\gamma}\eta,\qquad
\int_{C'}\eta=2\int_{\tilde\gamma}\eta.$$
数值（60 位）：$\int_{\tilde\gamma}\eta=-1.8809066251248561\ldots=-2\pi b_{17}$，
$\tilde n(1)=(J_1-J_2)/\pi\cdot(-1)=+b_{17}$（与 §8.4 表、Samart 一致）。

## 5. 除子与金刚石积（精确代数，多重验证）

$S_1$ 三次模型，群单位元 $O=[0:1:0]$，$A=(0,0)$。局部展开（$\partial F/\partial y(A)=1\ne0$，
$\partial F/\partial Z(O)\ne0$ 等）给出阶数；弦切法群律（$P\oplus Q=(P\cdot Q)\cdot O$）精确给出
$$2A=(0,-1)=:Q,\qquad 3A=T=[1:-1:0],\qquad 4A=O\quad(\text{与 PARI ellorder}=4\text{ 一致}).$$
$$\operatorname{div}(x)=[A]+[2A]-[O]-[3A],\qquad
\operatorname{div}(y)=3[A]-2[O]-[3A]\qquad(\text{形式与 k=0 全同}).$$
**显式双有理映射**（新）：$X=-(x+y)$、$Y=x(x+y)$ 把 $S_1$ 映到
$E_1:Y^2+XY-Y=X^3-X^2$（与 PARI `ellfromeqn` 输出一致），
$A\mapsto(0,0)$、$Q\mapsto(1,0)$、$T\mapsto(0,1)$、$O\mapsto\infty$；
在 $E_1$ 上独立重算 $\operatorname{div}(x)$、$\operatorname{div}(y)$（$x=-Y/X$、$y=(Y-X^2)/X$）
**逐项一致**。

金刚石积（定义同 LSZ cond21：$(f)\diamond(g)=\sum m_i n_j(S_i-T_j)\in\mathbb Z[E]^-$；
`code/k1_diamond.py` 精确展开 12 对，k=0 作对照组复现 $6(O)+5(A)-5(2A)$）：
$$(x)\diamond(y)=6(O)+2(A)-6(2A)-2(3A)\;\equiv\;6(O)+4(A)-6(2A)\quad\text{in }\mathbb Z[E]^-,$$
$$\Longrightarrow\quad D_E\bigl((x)\diamond(y)\bigr)=4D_E(A)-6D_E(2A)=4D_E(A).$$

## 6. $D_E$ 值（定理 + 60 位数值）

格基 $(w_1,w_2'=-w_2)$，$\tau=0.8873941733731784\ldots i$，
$q=e^{2\pi i\tau}=0.0037889663353792\ldots$（正实数）。
PARI `ellpointtoz`（`code/k1_zvals.gp`）：$z(A)/w_1=\frac14+\frac{\tau}{2}$
（$A$ 在 $E(\mathbb R)$ 的第二分支上），$z(2A)/w_1=\frac12$。

- **$D_E(2A)=0$（严格）**：2-扭点，$z_q=e^{\pi i}=-1$、$q$ 正实，
  级数逐项 $D(\text{负实数})=0$（Bloch–Wigner 性质）。
- **$D_E(A)=\dfrac{17}{8\pi}L(E,2)$（已发表定理 + 60 位）**：
  Lalín–Ramamonjisoa（*The Mahler measure of a Weierstrass form*, IJNT 13 (2017), 已存
  `literature/lalin-ramamonjisoa-cond17.pdf`）证明
  $L(E_{17},2)=\frac{8\pi}{17}D^E(P)$（其 §5，由其 (6)、(13) 与函数方程推出；
  $P$ 为其模型的 4-扭生成元）。其 $D^E$ 定义（其 Def. 5, eq. (10)：
  $\sum_{n\in\mathbb Z}D(q^n z)$）与我们的级数实现（`code/k1_dilog.py`）
  **归一化逐字一致**。数值：$D_E(A)=0.4702266562812140308\ldots$
  $=\frac{17}{8\pi}\cdot0.6951814630094846135\ldots$（60 位吻合）。
  与 k=0 的 Brunault (3.151)（$D_E(P)=\frac{11}{10\pi}L(E,2)$）同一模式
  $D_E(A)=\frac{N}{2t\,\pi}L(E,2)$（$t$=扭阶）。

合成：
$$D_E\bigl((x)\diamond(y)\bigr)=4\cdot\frac{17}{8\pi}L(E,2)
=\frac{17}{2\pi}L(E,2)=2\pi b_{17}=1.8809066251248561\ldots\ (60\text{ 位}).$$

## 7. Bloch 定理与闭环

采用 Lalín 教科书式陈述（`literature/lalin-ramamonjisoa-cond17.pdf` Thm. 6，
引 Bloch [Bl00]）：对 $\{x,y\}\in K_2(E)\otimes\mathbb Q$ 与
$H_1(E,\mathbb Z)^-$ 生成元 $\gamma^-$，
$$\int_{\gamma^-}\eta(x,y)=\pm\,D^E\bigl((x)\diamond(y)\bigr).$$
（ tempered 前提：$S_1$ Newton 面多项式 $x^3+y$、$x^3+x^2y$、$x^2y+y^2$、
$y^2+y$ 全分圆（$xy$ 项是 Newton 多边形的内部点），同 k=0 的 S2 步。）
该归一化的 conductor-17 先例即 L–R 自身：其 $(x{-}1)\diamond(y{-}1)=8(P)$、
$D^E(8P)=8\cdot\frac{\pi b_{17}}2=4\pi b_{17}=\int_{\gamma_1}\eta$，配合 Zudilin 的
$m(P_i)=2b_{17}$ 严丝合缝——**因子 1 在本曲线上有已发表证明背书**。

于是
$$\int_{\gamma^-}\eta=\pm\,4D_E(A)=\pm2\pi b_{17}
\qquad(\text{60 位实测 }+2\pi b_{17},\ \text{完全一致}).$$

## 8. 合成 $\tilde n(1)=b_{17}$ 与关于 k=0 因子 2 的存疑

$$\int_{\tilde\gamma}\eta=\frac12\int_{C'}\eta
=\frac12\cdot(\pm2)\cdot\bigl(\pm4D_E(A)\bigr)
=\pm\frac{17}{2\pi}L(E,2)=\pm2\pi b_{17},$$
符号由一次数值评估钉死（$\int_{\tilde\gamma}\eta=-2\pi b_{17}$）。再由结构恒等式
（k=0 第一波同型，逐点 $\log|y_{\mathrm{big}}|=-\log|y_{\mathrm{small}}|$ 推出）
$$\tilde n(1)=\frac1{2\pi}\Bigl[\int_{-c}^{c}-\int_{|θ|>c}\Bigr]\log|y_{\mathrm{big}}|\,d\theta
=-\frac{J_1-J_2}{\pi}=-\frac{1}{2\pi}\int_{\tilde\gamma}\eta=+b_{17}:$$
$$\boxed{\ \tilde n(1)=b_{17}\ }\qquad\text{（conductor 17 类比猜想，证毕；认证级别同 k=0）}.$$

**存疑（不影响 k=1 的证明，但建议 k=0 侧复核）**：k=0 笔记用
"π r = D_E(⋄)、r=(1/2π)∫η"（即 ∫η = ±2 D_E(⋄)，因子 2）配合金刚石积
$6(O)+5(A)-5(2A)$ 闭环。本次对 k=1 用**完全相同的**
（η 定义、$D_E$ Bloch 级数实现、⋄ 定义、$H_1^-$ 约定）得到
$D_E(\diamond)=2\pi b_{17}=\int_{\gamma^-}\eta$（因子 **1**，60 位），
与 Lalín Thm. 6 的陈述及 L–R 的 conductor-17 已证恒等式一致。
在因子-1 归一化下，k=0 的金刚石积需为 $10(A)-10(2A)$（恰为笔记值的 2 倍）
才能闭环；而 k=0 的除子数据经本次局部展开独立复算无误。
可能解释：k=0 笔记转述 Brunault 时 $r$ 的归一化与 Lalín 差因子 2，
即"因子 2"是转述层的人为产物而非定理内容。**建议用本波的多重验证方法
（显式 Weierstrass 映射 + 精确 ⋄ 展开 + Lalín 归一化）重核 k=0 的 §9.1。**

## 9. 认证级别小结

| 环节 | 状态 |
|---|---|
| 曲线识别（conductor 17、$\mathbb Z/4$、$(0,0)$ 4-扭、$w_{\mathrm{anti}}=w_2$、$b_{17}$） | PARI 80 位，精确 |
| 角点 $c=2\pi/3$、跳跃值 $y^2=-1$ | **精确代数**（$\omega^2+\omega+1=0$） |
| 闭链 $C'=\tilde\gamma+\beta_0$（闭、反不变、整系数） | 精确代数（同 k=0） |
| $\mathrm{class}(C')=\pm2\gamma^-$ | **Arb 球算术铁证**（比值球含 $-2$、半径 $3.3\times10^{-14}<1/2$；$w_{\mathrm{anti}}$ 经 Carlson RF 独立认证） |
| $\int_{\beta_0}\eta=\int_{\tilde\gamma}\eta$ | 精确积分代数 |
| $\operatorname{div}(x),\operatorname{div}(y)$、$(x)\diamond(y)=6(O)+4(A)-6(2A)$ | 精确代数（局部展开 + 显式双有理映射 + PARI + 代码展开，四重） |
| $D_E(2A)=0$ | 严格（2-扭，级数逐项为零） |
| $D_E(A)=\frac{17}{8\pi}L(E,2)$ | **已发表定理**（Lalín–Ramamonjisoa 2017）+ 60 位独立复核 |
| Bloch（∫η = ±$D^E(\diamond)$，Lalín Thm. 6 归一化） | 已证定理（conductor-17 先例：L–R + Zudilin） |
| tempered（$\{x,y\}\in K_2(E)\otimes\mathbb Q$） | 面多项式全分圆（含 $\Phi_3$），同 k=0 |
| $\kappa=1$（三次模型↔极小模型周期格一致） | `ellfromeqn` 模型 $\Delta=17=\Delta_{\min}$（$u^{12}=1$）；数值 16 位锁定（$\kappa=-1$ 只翻转比值符号，不影响整数性判定） |

## 10. 本波新文件

- `code/k1_pari.gp` — conductor 17 曲线识别（模型、扭子、周期、$b_{17}$）
- `code/k1_points.gp` — k=0/k=1 ellfromeqn 对照、扭点、角点 $(\omega,i)$ 验证
- `code/k1_zvals.gp` — 扭点 $z$ 值、$\tau$、$q$（$D_E$ 输入）
- `code/k1_certify.py` — mpmath 60/80 位认证（链积分、比值、跳跃值）→ `notes/attack11-k1-certify.txt`
- `code/k1_interval.py` — **Arb 球算术铁证**（含首段分支换算修正）→ `notes/attack11-k1-interval.txt`
- `code/k1_debug.py` — 首段分支符号问题的诊断脚本（存档）
- `code/k1_dilog.py` — $D_E(A)$、$D_E(2A)$、合成常数（60 位）
- `code/k1_diamond.py` — 金刚石积精确展开（k=1 + k=0 对照组）
- `notes/attack11-k1.txt` — regulator 侧输出汇总
- `literature/lalin-ramamonjisoa-cond17.pdf` — L–R 2017（conductor 17 已证公式与 Bloch Thm. 6 出处）
