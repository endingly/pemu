# 电子能量方程

## 1. M13 范围

M13 在 `physics`/`equation` 层增加独立电子能量输运，并由 simulation workflow 将它与
电子连续性、电势和反应源联合推进。
守恒变量是电子能量密度

$$
w_e = n_e\bar\varepsilon_e,
$$

其中 $n_e$ 是电子数密度，$\bar\varepsilon_e$ 是平均电子能量。厘米制网格的规范数值
单位分别为 $\mathrm{cm^{-3}}$、$\mathrm{eV}$ 和 $\mathrm{eV\,cm^{-3}}$。
`computeElectronEnergyDensity()` 与 `computeElectronMeanEnergy()` 完成两种表示间的转换；
后者在 $n_e\le n_\mathrm{floor}$ 的真空单元中定义 $\bar\varepsilon_e=0$，避免除以趋零密度。

电子能量方程本身不计算反应系数；M14 在独立的 chemistry 层以这里产生的平均能量或电场
状态计算反应率。电极/壁面能量通量仍属于 M15。

## 2. 控制方程与闭合

求解器推进

$$
\frac{\partial w_e}{\partial t}+\nabla\cdot\boldsymbol\Gamma_w=S_w,
$$

Simulation 将源项拆成内建电场功和调用者提供的附加碰撞/外部能量交换，

$$
S_w=-\boldsymbol\Gamma_e\cdot\mathbf E+R_w.
$$

其中 $oldsymbol\Gamma_e$ 由电子连续性方程使用的同一个 Scharfetter--Gummel 通量算子
计算。$R_w$ 由 `ElectronEnergyAdditionalSourceEvaluator` 提供，使反应能损模型仍保持在
equation 之外。

能量通量使用漂移—扩散闭合

$$
\boldsymbol\Gamma_w=\mathbf v_w w_e-D_w\nabla w_e.
$$

默认采用 Maxwellian、常动量传递频率和常压近似：

$$
\mu_w=\frac53\mu_e,\qquad
D_w=\frac53D_e,\qquad
\mathbf v_w=\frac53\mathbf v_e.
$$

`ExplicitElectronEnergyStepper` 允许显式传入其他正的 closure factor，但默认值固定为
$5/3$。电子面漂移速度由外部持有，求解器在每次计算前刷新能量漂移速度，因此能观察到
电静力更新后的值。

## 3. 空间离散

控制体内使用守恒有限体积散度。owner $P$ 到 neighbor $N$ 的内面通量采用现有
Scharfetter–Gummel 算子：

$$
\Gamma_{w,f}=\frac{D_w}{d_{PN}}
\left[B(-\mathrm{Pe}_f)w_P-B(\mathrm{Pe}_f)w_N\right],
\qquad
\mathrm{Pe}_f=\frac{v_{w,n,f}d_{PN}}{D_w},
$$

$$
B(x)=\frac{x}{e^x-1}.
$$

当漂移为零时它退化为中心扩散；强漂移时趋近迎风；指数平衡状态下通量严格为零。
由于 $v_w$ 与 $D_w$ 同乘 closure factor，面 Péclet 数保持与电子粒子输运一致，而整个
能量通量乘以该 factor。

普通边界接受有限、非负的能量密度 Dirichlet 值。M15 壁面可以改用线性法向通量
$\Gamma_w\cdot n=v_{w,\mathrm{loss}}w-Q_{\mathrm{in}}$；能量损失速度独立于粒子损失速度，
二次电子携带的能量由 wall assembler 根据入射粒子通量、产额和发射平均能量组装。

离散散度为

$$
(\nabla\cdot\boldsymbol\Gamma_w)_P=
\frac{1}{V_P}\sum_{f\in\partial P}s_{P,f}\Gamma_{w,f}A_f,
$$

其中 owner 的 $s_{P,f}=+1$，neighbor 的 $s_{P,f}=-1$，因此内面贡献全局严格抵消。

电场功不使用仅含漂移项的近似 $\mu_en_e|E|^2$，而由完整 SG 粒子通量的面法向分量重构：

$$
S_{E,P}=-\frac1{V_P}\sum_{f\in\partial P}
\Gamma_{e,n,f}E_{n,f}A_fd_{P,f},
$$

其中 $d_{P,f}$ 是单元中心到面中心沿全局面法向的正距离。$\Gamma_{e,n,f}E_{n,f}$ 在面法向
翻转时不变，因此 owner 和 neighbor 使用相同的面功率符号。该公式沿用 SG 算子的正交/
近正交网格假设。在厘米制裸数值中，粒子通量乘 $\mathrm{V/cm}$ 后直接得到
$\mathrm{eV/(cm^3\,s)}$，无需再乘基本电荷；若目标能量使用焦耳才需要乘 $e$。

## 4. 时间推进与稳定性

时间上采用显式 Euler：

$$
w_P^{k+1}=w_P^k+\Delta t
\left[S_{w,P}-(\nabla\cdot\boldsymbol\Gamma_w)_P^k\right].
$$

SG 离散的对角自损失率记为 $\lambda_P$。输运 CFL 上限是

$$
\Delta t_\mathrm{transport}=\min_{\lambda_P>0}\frac{1}{\lambda_P}.
$$

面向未来自适应推进，还提供不利用正源项和入流补偿的保守、状态相关非负上限

$$
\Delta t_\mathrm{positive}=
\min_P\frac{w_P}
{\lambda_Pw_P+\max(-S_{w,P},0)}.
$$

当 $w_P=0$ 且该单元没有负源时，这个状态相关上限不约束输运，因此调用者必须始终将它
与 $\Delta t_\mathrm{transport}$ 取最小值。

`computeIncrement()` 只形成未经稳定性判断的增量；`computeStableIncrement()` 同时检查
两个上限及完整候选状态，但不改变原场，方便上层先验证多个方程；`step()` 则完成同样
检查并提交。因此失败不会留下部分更新的能量场。

## 5. Simulation 联合推进

固定步长与自适应步长 Simulation 都要求显式传入 `ElectronEnergyConfiguration`，不存在
关闭电子能量方程的兼容路径。配置引用调用方持有的 $w_e$，并且必须指定电子 species、
能量 Dirichlet 边界及附加源项 evaluator。构造阶段会拒绝空 evaluator、非负电荷
或非漂移—扩散 species、不同 mesh、非法 density floor 及非法初始能量。

一个时间步按以下顺序执行：

1. 由 $n^k$ 更新电势、电场和电子漂移速度；
2. 由 $(w_e^k,n_e^k)$ 更新平均能量，随后计算反应率和物种源；
3. 使用同步后的电子 SG 粒子通量和电场组装 $-\boldsymbol\Gamma_e\cdot\mathbf E$，再由
   `ElectronEnergyAdditionalSourceEvaluator` 追加碰撞或外部源，形成 $S_w^k$；
4. 在修改状态前计算并验证全部物种与能量增量；固定步长源项若会产生负密度，则整个
   耦合状态保持不变；
5. 提交 $n^{k+1}$ 和 $w_e^{k+1}$，再刷新平均能量。

固定步长必须同时满足物种和能量限制，否则整个 workflow 进入 failed。自适应步长取
物种/能量输运限制的最小值、物种/能量非负限制的最小值，再由统一 controller 施加 safety、
增长率、用户上限和终点截断。因而 `lastTimeStepProposal()` 记录的是完整耦合系统的限制。

dump 会同步写出 `electron_energy_density` 和 `electron_mean_energy`。checkpoint schema v2
将电子 species 身份编码进能量字段 key，并与 species density、clock 及自适应步长历史在
同一次事务中保存/恢复；恢复先写临时场并验证平均能量，全部通过后才提交。

## 6. 当前限制

- $D_e$ 与 closure factor 在一个 stepper 生命周期内为空间常数；
- SG 两点通量沿用现有正交/近正交网格假设；
- 碰撞损失仍由调用者的附加源 evaluator 给出，尚无反应类型到能损的通用数据库；
- 电场功与连续性更新调用同一 SG 通量实现，但当前分别求值，尚未缓存全部物种面通量；
- 当前 wall closure 使用调用者提供的有效粒子/能量损失速度；尚未内置能量相关反射、
  thermionic emission、field emission 或 dielectric surface-charge 方程。

壁面模型、离散算子和 Simulation 生命周期仍保持单向依赖，具体数值方法见
[`plasma-wall-boundary.md`](plasma-wall-boundary.md)。

## 7. 参考模型

- [Two-Dimensional Self-Consistent Radio Frequency Plasma Simulations Relevant to the Gaseous Electronics Conference RF Reference Cell](https://pmc.ncbi.nlm.nih.gov/articles/PMC4887236/)
- [High order fluid model for streamer discharges. II](https://arxiv.org/abs/1302.4115)
- [Fluid modeling of low-temperature plasmas](https://doi.org/10.1063/5.0270714)
