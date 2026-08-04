# 文献笔记：Boyd 的 conductor 11 Mahler 测度猜想

## 记号

- $m(P)$：二元 Laurent 多项式的对数 Mahler 测度，$m(P)=\int_0^1\!\int_0^1 \log|P(e^{2\pi i u},e^{2\pi i v})|\,du\,dv$。
- $E_{11}=X_1(11)$：conductor 11 的椭圆曲线（LMFDB 11.a3），极小模型 $y^2+y=x^3-x^2$，等价写法 $y^2+y+x^3+x^2=0$（平移后）。
- 关联模形式 $f_{11}(\tau)=\eta(\tau)^2\eta(11\tau)^2=q\prod_{n\ge1}(1-q^n)^2(1-q^{11n})^2=\sum a_n q^n$，权 2，$\Gamma_0(11)$。
- $b_{11}:=L'(E_{11},0)=\dfrac{11}{4\pi^2}L(E_{11},2)$（函数方程，根数 $+1$）。数值 $b_{11}=0.1521471\ldots$

## 背景

- Deninger (1995)：$m(1+x+1/x+y+1/y)\stackrel?=L'(E_{15},0)$，源自 Bloch–Beilinson 猜想。
- Boyd (1998, Experimental Math. 7:1)：系统数值实验，给出大量 $m(P)=r\cdot b_N$ 型猜想（$r\in\mathbb{Q}$），按 conductor 分类；最小可能的 conductor 为 11, 14, 15, 17, 19, 20, 21, 24。
- Rodriguez-Villegas (1997)：此类猜想可从 Bloch–Beilinson 猜想推出；CM 情形可严格证明。

## conductor 11 的恒等式及状态

### 已证明

1. **Brunault (2005/2006，博士论文 + 论文，参数化 $X_1(11)$ by modular units，Beilinson 定理的显式版本）**：
   $$m\big((1+x)(1+y)(1+x+y)+xy\big)=\frac{77}{4\pi^2}L(E_{11},2)=7\,L'(E_{11},0)=7b_{11}.$$
   （Bertin–Lalín 综述 [Br05, Br06]；Brunault 论文集 §3.7/3.9。）

2. **Brunault (2006)**：
   $$m\big(y^2+(x^2+2x-1)y+x^3\big)=5b_{11}.$$
   （Boyd PNWNT 2015 slides, p.28。）

### 仍开放（本次攻击目标）

3. **Boyd 1998 编号 (2-33)，族 $S_k=y^2+(x^2+kx+1)y+x^3$（Boyd 记号 $S_\alpha=y^2+(x^2+\alpha x+1)y+x^3$？按 Samart 2023 §4：$S_\alpha=y^2+(x^2+\alpha x+1)y+x^3$，$K\cap\mathbb R=[-4,2]$）。**
   $k=0$：$S_0=y^2+(x^2+1)y+x^3$，零点曲线即 conductor 11 椭圆曲线 $E$。
   $S_0$ 在 2-torus 上有零点（$x=\pm i$），故 $m(S_0)$ **不是** $b_{11}$ 的有理倍数
   （Boyd：$m(S_0)=0.4056029\ldots$，"seemingly not $r b_{11}$"）。
   但**劈裂积分**（沿 branch cut 的积分）猜想为：
   $$\frac1\pi\int_0^{\pi/2}\log|y_-(e^{i\theta})|\,d\theta-\frac1\pi\int_{\pi/2}^{\pi}\log|y_-(e^{i\theta})|\,d\theta\stackrel?=-L'(E,0)=-b_{11},$$
   其中 $y_\pm(x)=-\dfrac{x^2+1}{2}\pm\sqrt{\dfrac{(x^2+1)^2}{4}-x^3}$ 为 $S_0=0$ 的两根。
   （Samart 2023, arXiv:2301.05390, eq. (4.1)：明确标注为 Boyd 的数值猜想，尚未证明；并指出可尝试用其 §3 方法 + Brunault 的 conductor 11 结果证明。）
   相关地，Boyd slides：对模型 $y^2+y+x^3+x^2=0$，带符号积分 $\frac1\pi\int_0^\pi\log|y_2(t)|dt=0.1521471\ldots=b_{11}$ 到 50 位。

4. Samart 2023 还发现 $k=1$（conductor 17）与 $k=-1$（conductor 53）的类似猜想恒等式——本课题只聚焦 conductor 11。

## 攻击思路（2 小时）

1. 高精度（50+ 位）计算 $b_{11}$：用 $f_{11}=\eta(\tau)^2\eta(11\tau)^2$ 的系数 $a_n$ 与近似函数方程
   $L(E,2)=\sum_{n\ge1}a_n e^{-2\pi n/\sqrt{11}}\big(\frac{1}{n^2}+\frac{2\pi}{n\sqrt{11}}\big)\cdot(\ldots)$（权 2、根数 +1 的标准公式），乘以 $11/(4\pi^2)$。
2. 高精度计算劈裂积分 $I_-=\frac1\pi\int_0^{\pi/2}-\frac1\pi\int_{\pi/2}^\pi$，与 $-b_{11}$ 比对（mpmath，1D 积分，快）。
3. 验证已证明恒等式 1、2（Jensen 公式降为 1D：$m(P)=\frac1{2\pi}\int(\log|A|+\sum_j\log^+|y_j|)d\theta$）。
4. PSLQ：在 $\{I_-, b_{11}\}$、$\{m(S_0), b_{11}\}$ 及更多常数上搜索小有理关系。
5. 探索证明路线：Brunault–Mellit–Zudilin regulator 公式 + 路径闭包分析（写成 proof strategy）。

## 文献（literature/）

- `bertin-lalin-survey.pdf` — Bertin–Lalín, Mahler measure of multivariable polynomials（综述，含 conductor 11 状态）
- `boyd-pnwnt2015.pdf` — Boyd, PNWNT 2015 slides（含 $m(S_0)$ 数据与 Brunault 结果）
- `brunault-these.pdf` — Brunault 博士论文（$X_1(11)$ 的 Beilinson 显式版本）
- `zudilin-regulator.pdf` — Zudilin, Regulator of modular units and Mahler measures (BMZ 公式)
- `samart2023.pdf` — Samart, Mahler measure of a nonreciprocal family（含开放猜想 (4.1)）
- `lalin-samart-zudilin-cond21.pdf` — conductor 21 的类似方法（方法参考）
- `boyd-slides.pdf` — Boyd 关于 Mahler 测度与 L(E,3) 的 slides

## 第三波补充（2026-08-04）

### k=2 (conductor 37)、k=3 (conductor 79) 恒等式的文献状态
- Lalín–Ramamonjisoa "The Mahler measure of a Weierstrass form"（umontreal.ca/~mlalin/Mahler-Weierstrass.pdf）：
  证的是 F_k=y²+kxy-x³-x 族（conductor 17 的 m(F_3)=7/2 L'(E17,0)），并综述 Rodriguez-Villegas 的
  conductor 37 恒等式 7m(W₁)=5m(W₂)（W₁,W₂ 为 E37 的**其他** Weierstrass 模型，猜测各=35L'(E37,0)）。
  **我们的 S_2=y²+(x²+2x+1)y+x³ 是 E37 的非 Weierstrass 模型，不在其列。**
- conductor 79：未检索到任何 Mahler 测度结果。
- 结论：m(S_2)=2|L'(E37,0)|、m(S_3)=|L'(E79,0)| 很可能仅见于 Boyd 1998 数值表
  （projecteuclid 被 Incapsula 拦截，未能核对原表），文献中未见证明。S_2/S_3 无 torus 交点，
  属标准 BMZ 适用区，是"可证且可能未证"的候选。

### BMZ Theorem 1 的适用范围核对（S4）
Zudilin "Regulator of modular units and Mahler measures" Thm 1（Mellit–Brunault）：
modular units g_a,g_b 沿**尖点 c/N 到 i∞** 路径的 regulator 积分 = 4L(f,2)，f 为显式 Eisenstein 组合；
若组合为尖点形式则得曲线的 L 值。**定理处理的是尖点间开路径**；我们的链 γ̃ 是闭链，
需要 Beilinson–Brunault（闭链配对版）+ γ̃ 在 H_1(E,Z)^- 中的类。H_1^- 秩 1 ⇒ γ̃ = n·γ⁻，
n 由周期配对确定。

### 绕数结果（code/winding.py + winding.gp，60 位）
- 几何路径 |x|=1 携带 y_big 在折点 θ=±π/2 **不连续**（(i,-e^{-iπ/4}) ↔ (i,+e^{iπ/4}) 跳跃）；
  朴素环积分 -0.47447i 不是任何整闭链的周期（比值 0.16262 非有理）。
- 带符号链周期 I_signed = -2.917633233876990458...i ≡ PARI 11.a3 的 w_anti = 2i·Im ω₂
  （比值 1 到 13 位，误差来自 θ=0 处 √θ 奇性）。**γ̃ = H_1(E,Z)⁻ 生成元（绕数 n=1）**。
- 附带：四次模型 dx/u 的模型常数 κ=1（Om_re = PARI w1 精确相等）。
- 含义：regulator 配对 r({x,y})[γ̃] = r({x,y})[γ⁻]，Beilinson–Brunault 给出显式有理倍数；
  数值预测 r({x,y})[γ⁻] = 2π b_11（即纲要中的 r=2），无自由同调 ambiguity。
