# Plasma benchmark 与验证边界

## 1. M16 的目标

M16 区分 code verification 与 physical validation：前者检查项目是否正确求解其声明的方程，
后者检查这些方程是否足以预测实验。解析解、时间步加密和多代码公开数据可以验证实现；只有与
有不确定度说明的实验对比，才能进一步验证物理模型。当前首先建立可在 CI 中重复运行的
full-workflow verification，不把 400 V 集成回归或一次 `SolverResult::success()` 当成物理验证。

所有 M16 时间 benchmark 都通过公开的 `FixedStepPlasmaSimulation` 运行完整的 Poisson、反应、
粒子输运、电子能量和时钟事务，而不是直接调用某个私有更新函数。空间状态选成严格均匀，使
内部 SG 通量散度为零；全部外边界由 wall flux law 接管，从而得到封闭且有解析解的常微分方程。

## 2. 均匀电离—电子能量 benchmark

考虑电子和单正离子成对生成，目标中性粒子密度和反应系数均为常数：

$$
R=k n_N n_e=\lambda n_e,
\qquad
\frac{\mathrm d n_e}{\mathrm dt}
=\frac{\mathrm d n_i}{\mathrm dt}=\lambda n_e.
$$

初态满足 $n_e=n_i=n_0$，所以电荷密度与电场严格为零。电子能量附加源令每个新生电子携带
当前平均能量：

$$
S_w=\bar\varepsilon_e R,
\qquad
w_e=n_e\bar\varepsilon_e.
$$

因此 $\bar\varepsilon_e$ 保持常数，连续解为

$$
n_e(t)=n_i(t)=n_0e^{\lambda t},
\qquad
w_e(t)=\bar\varepsilon_e n_0e^{\lambda t}.
$$

固定步显式 Euler 的离散真值则为

$$
n^{m}=n_0(1+\lambda\Delta t)^m,
\qquad
w_e^m=\bar\varepsilon_e n^m.
$$

测试同时逐项比较离散真值，并以 $\Delta t=0.1,0.05,0.025$ 对连续解计算

$$
p=\log_2(e_{\Delta t}/e_{\Delta t/2}),
$$

要求 $p>0.9$。它覆盖 Te 查表反应率、逐步反应重算、化学计量源、电子能量源、耦合原子提交和
固定时钟，而不是只验证一个标量 Euler 公式。

## 3. 壁面损失—二次发射 benchmark

两单元网格中每个单位控制体都有三个单位长度外边界面，故
$A_{\partial P}/V_P=3$。令

$$
a_e=v_eA_{\partial P}/V_P,
\quad
a_i=v_iA_{\partial P}/V_P,
\quad
a_w=v_wA_{\partial P}/V_P,
$$

离子诱导二次电子产额为 $\gamma$，发射平均能量为 $\varepsilon_\mathrm{emit}$。均匀状态满足

$$
\dot n_i=-a_i n_i,
$$

$$
\dot n_e=-a_e n_e+\gamma a_i n_i,
$$

$$
\dot w_e=-a_w w_e
+\varepsilon_\mathrm{emit}\gamma a_i n_i.
$$

当三个损失率互异时，连续解为

$$
n_i(t)=n_{i0}e^{-a_it},
$$

$$
n_e(t)=n_{e0}e^{-a_et}
+\gamma a_i n_{i0}\frac{e^{-a_it}-e^{-a_et}}{a_e-a_i},
$$

$$
w_e(t)=w_{e0}e^{-a_wt}
+\varepsilon_\mathrm{emit}\gamma a_i n_{i0}
\frac{e^{-a_it}-e^{-a_wt}}{a_w-a_i}.
$$

壁面损失会产生空间电荷和电场，而 Simulation 必须自动加入电子场功。为了让这个 benchmark
只比较上述壁面 ODE，附加能量源使用同一个粒子通量和电场重构场功，再写入其相反数。这是明确
的 manufactured-source 选择，不代表生产模型应取消 Joule heating。

测试用 $\Delta t=0.05,0.025,0.0125$ 同时比较逐步 Euler 递推和连续解析解，并要求三场最大
误差的观测阶 $p>0.85$。该判据覆盖主粒子损失、二次电子、发射能量、每步 wall 刷新以及
species/energy 的联合提交。

## 4. 公开 streamer benchmark 的接入审计

Bagheri 等人的六代码比较是合适的公开 code-verification 目标：论文和 CC0 数据集给出了三个
轴对称正 streamer 算例。Case 1 使用较高背景电子/离子密度且不含光电离，是最小候选；论文也
明确区分了 code verification 与实验 validation。

当前不能把 pemu 的平面 400 V 算例直接与该数据比较。Case 1 的普通物种边界已补齐
Dirichlet/Neumann SG 通量契约；其中 Neumann 值表示外向扩散通量，零漂移下的齐次
Neumann 即 benchmark 所需 zero-flux。`AxisymmetricMeshView` 也已提供绕 `x=r=0`
旋转后的控制体体积与面面积，并复用现有有限体积方程。仍需：

- 随局部 $|E|$ 变化的电子 mobility 与 diffusivity，而当前 `SpeciesProperties` 在 stepper
  生命周期内保存空间常数；
- 论文给出的局部场 ionization/attachment 系数及相同初态、诊断量和网格/时间收敛流程。

场强相关输运不能直接以每个面的 $|E_n|$ 查表：同一单元的径向面与轴向面会看到不同法向
分量，而 mobility/diffusivity 是局部 $|\mathbf E|$ 的标量函数。M16.4.1 因此按以下数据流
实现：先由面法向电场重构单元 $|\mathbf E|$，在单元查表，再以明确的几何插值生成面系数；
电子粒子 SG 和电子能量 SG 必须在同一时间层共同刷新，后者使用相同动态 $D_e$ 乘能量
closure factor。轴对称视图的 $r=0$ 面面积为零，场强最小二乘重构会忽略该面的零权重，
但仍要求其余正面积面的法向张成二维空间。

Case 3 还需要非局部 photoionization Helmholtz 方程，因此不作为第一个接入目标。上述剩余项
完成后，M16 才能加入 Case 1 的公开数据对比；在此之前，M16 保持“进行中”，不会用内部 golden
number 替代独立参考。

## 5. 参考资料

- B. Bagheri et al., “Comparison of six simulation codes for positive
  streamers in air,” *Plasma Sources Science and Technology* 27 (2018)
  095002, <https://doi.org/10.1088/1361-6595/aad768>；
- 对应 CC0 benchmark 数据集，DOI:
  <https://doi.org/10.17026/dans-x7r-266f>。
