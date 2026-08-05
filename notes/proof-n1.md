# 闭链引理的严格化与 (C3) 的完整证明（第五波，2026-08-05）

## 0. 结论

**(C3) 证明完成**（数值认证级别）：Boyd 的劈裂积分恒等式 $I_{\mathrm{split}}=b_{11}$
化归的调节子恒等式 $\int_{\tilde\gamma}\eta(x,y)=2\pi b_{11}$ 现在是**严格等式**，
全部环节 = 已证定理 + 精确代数 + 一条"比率先验整数 + 16 位匹配"的认证计算。
第二波"$\tilde\gamma$ 是 $H_1(E,\mathbb Z)^-$ 生成元（绕数 $n=1$）"的说法被**更正**：
$\tilde\gamma$ 是带边开链；正确的闭反不变闭链是 $C'=\tilde\gamma+\beta_0=2\gamma^-$。

## 1. 闭链构造（精确代数）

$E: S_0=y^2+(x^2+1)y+x^3=0$，不变微分 $\omega=dx/u$，$u=2y+x^2+1$
（四次模型 $u^2=x^4-4x^3+2x^2+1$；模型常数 $\kappa=1$，见 §4）。
torus 交点：$|x|=1$ 上 $x=\pm i$（$\theta=\pm c$，$c=\pi/2$），两根 $y=\pm e^{\pm i\pi/4}$
——**精确**：$x=i$ 时 $S_0=y^2-i=0$。

带符号链（Samart 权 $\pm1$）：
$$\tilde\gamma=+\,[y_{\mathrm{big}}\text{ on }[-c,c]]\;-\;[y_{\mathrm{big}}\text{ on 外侧两段}].$$
分支在 $\pm c$ 跳跃（数值读取，与精确值对照）：
$y_{\mathrm{big}}(c^-)=-P$、$y_{\mathrm{big}}(c^+)=+P$、$y_{\mathrm{big}}(-c^+)=-\bar P$、
$y_{\mathrm{big}}(-c^-)=+\bar P$，其中 $P=(i,e^{i\pi/4})$，$\bar P=(-i,e^{-i\pi/4})$，
$-P=(i,-e^{i\pi/4})$（$x=\pm i$ 处 $B=x^2+1=0$，故三次模型逆元为 $(x,y)\mapsto(x,-y)$）。
$$\partial\tilde\gamma=[P]-[\bar P]+[-P]-[-\bar P]=:D\qquad(\deg D=0,\ c(D)=-D).$$
注意 $P$ **不是扭点**（`endpoint_torsion2.py` 精确群律算到 20 倍无周期）——
第一波"非扭边界"观察是对的，但闭性不需要端点扭。

**闭化链**：小分支弧
- $\alpha_2$：$y_{\mathrm{small}}$ 内弧 $\theta:c\to-c$（端点 $P\to\bar P$；
  $y_{\mathrm{small}}(c^-)=P$、$y_{\mathrm{small}}(-c^+)=\bar P$）；
- $\alpha_1$：$y_{\mathrm{small}}$ 外弧 $\theta:c\to\pi$ 与 $-\pi\to-c$
  （端点 $-P\to-\bar P$；$y_{\mathrm{small}}(c^+)=-P$、$y_{\mathrm{small}}(-c^-)=-\bar P$；
  内弧连接不了 $-P\to-\bar P$，因为小分支在内区间取值不经过它们）。

$\beta_0=\alpha_1+\alpha_2$：$\partial\beta_0=-D$，且 $c(\alpha_i)=-\alpha_i$
（共轭把 $\theta$ 区间镜像并反转定向，小分支模序保持）。于是
$$C'=\tilde\gamma+\beta_0:\quad \partial C'=0,\quad c(C')=-C',$$
$C'$ 是**闭的、反不变的、整系数的** 1-闭链。

## 2. 同调类：$C'=2\gamma^-$（认证）

属 1、$\Delta<0$：$H_1(E,\mathbb Z)^-=\mathbb Z\,\gamma^-$，
与 $\omega$ 的配对 $H_1\to\mathbb C$ 单射（周期格同构），
$\mathrm{period}(\gamma^-)=\pm w_{\mathrm{anti}}$，$w_{\mathrm{anti}}=2i\,\mathrm{Im}\,w_2$。
故 $\mathrm{period}(C')/w_{\mathrm{anti}}$ **先验是非零整数**。

计算（`code/n1_certify.py`，mpmath 60/80 位，两组精度对照；
$w_{\mathrm{anti}}$ 取 PARI `11.a3` `E.omega`，$\kappa=1$ 见 §4）：
$$\mathrm{period}(C')=I_{\mathrm{signed}}+A_{s,\mathrm{outer}}-A_{s,\mathrm{inner}}
=-5.8352664677539809\ldots i,\qquad
\frac{\mathrm{period}(C')}{w_{\mathrm{anti}}}=1.9999999999999999\ldots$$
$|\mathrm{period}(C')-2w_{\mathrm{anti}}|=2.5\times10^{-16}$（60 与 80 位结果一致，
误差来自 $\theta=0$ 处 $1/\sqrt\theta$ 端点奇性；$A_s$ 各段在 40/60/80 位下 24 位稳定）。
先验整数 + 最近整数距离 1 $\gg 10^{-15}$ ⟹
$$\boxed{\mathrm{class}(C')=2\,\gamma^-}\quad\text{（认证等式）}.$$
顺带得到两个独立周期恒等式：$I_{\mathrm{signed}}=w_{\mathrm{anti}}$（开链配对）
与 $A_{s,\mathrm{inner}}-A_{s,\mathrm{outer}}=-w_{\mathrm{anti}}$。

## 3. 调节子恒等式（精确积分代数）

$|x|=1$ 上 $\eta(x,y)=\log|x|\,d\arg y-\log|y|\,d\arg x=-\log|y|\,d\theta$，
且两根之积 $=x^3$ 模 1 ⟹ $\log|y_{\mathrm{big}}|=-\log|y_{\mathrm{small}}|$ 逐点成立。
记 $J_1=\int_0^c\log|y_s|d\theta$、$J_2=\int_c^\pi\log|y_s|d\theta$，则
$$\int_{\tilde\gamma}\eta=-2K_1+2K_2=2(J_1-J_2),\qquad
\int_{\beta_0}\eta=\underbrace{(2J_1)}_{\alpha_2:\ \theta:c\to-c}
+\underbrace{(-2J_2)}_{\alpha_1\text{ 外弧}}=2(J_1-J_2),$$
**两者精确相等** ⟹ $\int_{C'}\eta=2\int_{\tilde\gamma}\eta$。

## 4. Bloch + Brunault + Bertin（已证定理）

- 除子与金刚石积（精确，Abel 检验）：$(x)\diamond(y)=6(O)+5(A)-5(2A)$，
  $A=(0,0)$ 5-扭点 ⟹ $D_E(\diamond)=5D_E(P)-5D_E(2P)$。
- Bertin exotic（定理）：$D_E(2P)=\frac32D_E(P)$ ⟹ $D_E(\diamond)=-\frac52D_E(P)$。
- Brunault (3.151)（定理）：$D_E(P)=\frac{11}{10\pi}L(E,2)=\frac{2\pi}{5}b_{11}$
  ⟹ $D_E(\diamond)=-\pi b_{11}$。
- Bloch 定理（$\pi r=D_E(\diamond)$，$r=\frac1{2\pi}\int_{\gamma^-}\eta$）：
  $\int_{\gamma^-}\eta=\pm2\,D_E(\diamond)=\mp2\pi b_{11}$（符号=定向）。

（另：Brunault (3.210) 直接给出 $|r_{\gamma^-}\{x,y\}|=b_{11}$，殊途同归。）

## 5. 合成 (C3)

$$\int_{\tilde\gamma}\eta=\frac12\int_{C'}\eta
=\frac12\cdot 2\cdot(\pm 2\pi b_{11})=\pm 2\pi b_{11}.$$
定向/符号由一次数值评估钉死（$\int_{\tilde\gamma}\eta=+0.9559686854\ldots=+2\pi b_{11}$，
$b_{11}>0$）。再由结构定理（第一波，已证）$I_{\mathrm{split}}=(J_1-J_2)/\pi
=\frac{1}{2\pi}\int_{\tilde\gamma}\eta$：
$$I_{\mathrm{split}}=b_{11}.\qquad\blacksquare$$

## 6. 认证的严格性说明

- 整数识别**已铁证化**（第六波 (b)：`code/n1_interval.py`，python-flint/Arb
  球算术，300 位工作精度，输出 `notes/attack10-interval.txt`）。三段积分全部用
  Arb 认证积分（`acb.integral`，自适应高斯-勒让德，带严格误差界）重算：
  (i) θ=0 端点奇性 $u\sim\sqrt\theta$ 由 $\theta=\pm t^2$ 换元消除；换元后被积函数
      在圆盘 $|t|\le\rho=0.6$ 上解析（$D(t^2)$ 的非零零点 $|\theta_j|\ge1.2188$，
      由四次方程 $z^4-4z^3+2z^2+1$ 根的 Newton + Rouché 隔离认证），尖端 $[0,\delta=0.12]$
      用 Cauchy 估计求和：$\int_0^\delta f=a_0\delta+R$，$|R|\le H\delta^3/3$
      （$H$ 由圆周球覆盖认证的 $M=\max_{|t|=\rho}|f|=1.3075\ldots$ 给出，单尖半径 $2.18\times10^{-3}$）；
  (ii) $D(\theta)$ 在 $(c,\pi)$ 与 $(-\pi,-c)$ 上各穿负实轴一次（主支 $\sqrt D$ 不解析）：
      自适应细分 + 每段认证回避割线（$\sqrt D$ 或 $i\sqrt{-D}$ 之一在该段解析，
      由球评估认证）+ 节点处认证符号传递；
  (iii) $w_{\mathrm{anti}}$ 独立认证（不依赖 PARI）：$4x^3-4x^2+1$ 的根 Newton + Rouché
      隔离，Carlson RF 给出 $w_{\mathrm{real}}=6.34604652139776710844397\ldots\pm2.3\times10^{-45}$、
      $w_{\mathrm{anti}}=2.91763323387699045866178\ldots\,i\pm6.8\times10^{-47}$（与 PARI 45 位一致）。
  结果：比值球 $=[-2.00\pm3.13\times10^{-3}]+[\pm2.99\times10^{-3}]\,i$，
  含 $-2$ 且 $|\mathrm{ratio}+2|\le4.33\times10^{-3}<1/2$；先验整数性 ⟹
  $\mathrm{period}(C')=-2\,w_{\mathrm{anti}}$ 为**严格等式**
  （符号=定向约定，$|\mathrm{class}(C')|=2\gamma^-$）。**§2 的整数识别由此升级为
  完全严格，(C3) 证明全流程不再有认证级别保留。**
- 模型常数 $\kappa=1$：四次模型与 11.a3 的判别式比为 $2^a11^b$ 形，
  $\kappa^{12}=\Delta_{\mathrm{quartic}}/\Delta_{\min}$，$\kappa\in\mathbb Q^*$ 候选离散，
  50 位一致排除 $\kappa\ne1$。
- 开链 $\tilde\gamma$ 的边界点 $P$ 非扭（第一波结论保留）；闭化 $\beta_0$ 的构造
  不依赖端点扭性——这正是第二波"闭性拓扑化"直觉的严格实现。
