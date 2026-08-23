# E/N 与 Te 相关的电子碰撞化学

## 1. M14 的模型边界

M14 将电子能量和电场状态转换为局部反应系数，并由反应网络把反应率转换为物种源。reaction
不是独立构建模块；它属于 `physics`，所有反应领域类型统一位于 `pemu::physics::reaction`：

1. `Reaction`/`ReactionNetwork` 描述反应结构和净化学计量；
2. `TabulatedRateCoefficient` 描述状态相关的局部系数；
3. `massActionReactionRate()` 与 `MassActionReactionAssembler` 计算任意体质量作用律；
4. `ReactionNetwork::accumulateSources()` 根据净化学计量把反应率汇总为物种源。

电子温度换算仍属于电子能量物理；面法向量到单元向量的重构属于 `discretization`；
Simulation 中的 evaluator 只负责把 workflow 状态适配给 reaction API。这样既合并同源的反应
概念，也避免 reaction 反向拥有离散算子或 Simulation 生命周期。M14 的表格由调用者提供；
项目不内置特定气体或截面数据库，也不把任意示例表宣称为经过验证的真实化学集。

## 2. 等效电子温度

在各向同性 Maxwellian 电子能量分布假设下，

$$
\bar\varepsilon_e=\frac32 k_\mathrm{B}T_e.
$$

低温等离子体模型常把 $k_\mathrm{B}T_e$ 本身以 eV 表示，因此代码中的 `Te in eV` 定义为

$$
T_e[\mathrm{eV}]=\frac23\bar\varepsilon_e[\mathrm{eV}].
$$

`electronTemperatureEv()` 和 `computeElectronTemperatureEv()` 实现这一转换。它们返回的是
温度对应的能量，不是 kelvin；负值或非有限值会被拒绝。该换算只有在 Maxwellian 假设
适用时才有物理意义。若反应数据直接以平均电子能量为横坐标，调用者应直接查表，不应先
做这一步转换。

## 3. 约化电场

约化电场定义为

$$
\mathcal E=\frac{|\mathbf E|}{N},
$$

其中 $N$ 是碰撞目标的数密度。对项目当前的厘米制规范，

$$
1\,\mathrm{Td}=10^{-17}\,\mathrm{V\,cm^2},
$$

所以

$$
\mathcal E[\mathrm{Td}]
=\frac{|\mathbf E|[\mathrm{V/cm}]}{N[\mathrm{cm^{-3}}]}10^{17}.
$$

`reducedElectricFieldTownsend()` 的名称显式固定了输入和输出约定，避免把 SI 裸数值与厘米
制裸数值混用。$N$ 必须严格为正；真空单元没有有限的 $E/N$，不能用零密度静默代替。

Simulation 的电场存储在面上，且只有面法向分量。discretization 中的
`reconstructCellVectorMagnitudeFromFaceNormal()` 在每个控制体内求解面积加权最小二乘问题

$$
\underset{\mathbf E_P}{\operatorname{argmin}}
\sum_{f\in\partial P}A_f
\left(\mathbf n_f\cdot\mathbf E_P-E_{n,f}\right)^2.
$$

当面法向张成网格维数时，该方法对常 Cartesian 电场严格重构；退化控制体会被明确拒绝。
这里使用 mesh 的全局 owner-to-neighbor 面法向。邻居访问同一面时无需同时翻转法向和面值，
因为存储的方程仍是同一个 $\mathbf n_f\cdot\mathbf E=E_{n,f}$。

该重构是假设单元内电场为常量的局部闭合。高非正交网格或单元内强烈变化的电场需要更高阶
梯度重构和网格收敛验证。

## 4. 速率系数查表

`TabulatedRateCoefficient` 存储严格递增的一维横坐标和非负速率系数。横坐标可以是
$E/N$（Td）、$T_e$（eV）或直接的 $\bar\varepsilon_e$（eV）；类本身保持数值通用，具体
坐标和单位由构造它的物理模型声明。

默认使用对数线性插值：

$$
\ln k(x)=(1-\theta)\ln k_i+\theta\ln k_{i+1},
\qquad
\theta=\frac{x-x_i}{x_{i+1}-x_i}.
$$

反应系数跨越多个数量级时，这通常比直接线性插值稳定。对数线性模式要求全部系数严格为
正；确实包含零系数的表必须显式选择 `RateInterpolation::linear`。

越界默认抛出异常。`RateTableBounds::clamp` 可显式选择端点钳制，但它只是一种数值策略，
并不代表端点外的碰撞物理保持不变。实现不提供线性外推，因为外推可能产生负系数或快速
增长的非物理源项。

## 5. 任意体质量作用律

反应率核心接口不固定反应物数量，而计算

$$
R_r=k_r\prod_j n_j^{\alpha_{jr}}.
$$

`Reaction::stoichiometry` 保存净生成系数 $\nu_{sr}$，`Reaction::kinetic_orders` 单独保存
$\alpha_{jr}$。两者不能互相推导。例如三体附着

$$
e+O_2+M\rightarrow O_2^-+M
$$

中第三体 $M$ 的净化学计量为零，但动力学级数为一。`MassActionReactionAssembler` 在构造时
把这些 species ID 解析成稳定场引用，之后每个时间步不再分配内存。`MassActionTerm` 也允许
显式引用外部密度场或均匀 bath-gas 密度。

若总分子数为整数 $m$，所有数密度使用 $\mathrm{cm^{-3}}$，则

$$
[k_r]=\mathrm{cm}^{3(m-1)}\mathrm{/s}.
$$

所以一体、二体和三体系数依次为 $\mathrm{s^{-1}}$、$\mathrm{cm^3/s}$ 和
$\mathrm{cm^6/s}$。assembler 根据 `kinetic_orders`、密度单位和反应率单位推导所需系数单位，
不会接受用二体系数 metadata 驱动三体反应。含 fractional empirical order 的裸数值计算受
支持；由于当前运行时单位后端只在这里处理整数幂，带 metadata 的 fractional order 会被拒绝。

kernel 先验证全部单元及完整乘积，再覆写目标场。任意零密度因子会先把反应率确定为零，
避免数学上为零的反应因其他巨大因子产生中间溢出。`binaryReactionRate()` 仅保留为委托给
通用 kernel 的便利接口，不包含独立二体实现。

数值测试覆盖 Maxwellian 能量换算、Townsend 换算、常向量精确重构、查表插值、任意三体
组装、fractional order、第三体零净源、`cm^6/s` 校验、零因子短路和失败原子性。

## 6. Simulation 接入

`PlasmaReactionRateEvaluator` 接收 `PlasmaReactionRateContext`，其中包含同一时间层的物种密度、
电势、面法向电场和电子平均能量。固定步长与自适应 Simulation 都在更新电静力并刷新平均
能量之后调用该 evaluator，因此 E/N 和 Te 化学不会读取上一时间层的派生状态。命名 context
也避免了继续扩展多个同类型位置参数。

`TabulatedElectronImpactEvaluator` 是标准的一反应实现。它引用调用者持有的目标粒子密度场，
内部复用预分配的 $|E|$、查表坐标和 $k$ 工作场，不在每个时间步分配内存。配置可选择
`reduced_electric_field_townsend` 或 `electron_temperature_ev`：前者从当前面电场重构 $|E|$
并计算 E/N，后者从当前平均能量计算 Maxwellian Te；随后都执行同一个查表和通用质量作用律
流水线。该 evaluator 只覆写指定的 `ReactionId`，Simulation 会在调用前清零完整反应率集合，
多个反应可由调用者在一个组合 evaluator 中依次组装。

标准 evaluator 默认校验厘米制/eV metadata：电子和目标密度必须为 $\mathrm{cm^{-3}}$，
面电场必须为 $\mathrm{V/cm}$，平均能量必须为 eV，反应率必须为
$\mathrm{cm^{-3}s^{-1}}$。完全无 metadata 的合成问题只有显式传入
`allow_unitless_raw_values=true` 才可运行，此时仍按上述裸数值单位解释。混合有/无 metadata
或其他单位会直接失败，而不会把 SI 数值静默送入厘米制 Townsend 公式。

## 7. 当前限制

- 查表是一维的；同时依赖 $E/N$、气体温度和组分比例的模型需要更高维数据接口；
- 当前没有 Boltzmann solver 或截面到速率系数的在线积分；
- Maxwellian 的 $T_e=2\bar\varepsilon_e/3$ 不是任意非平衡 EEDF 的普遍关系；
- 面法向电场到单元幅值采用分片常量最小二乘重构；
- 带物理 metadata 的经验 fractional reaction order 尚不支持运行时分数次单位幂；
- 当前 wall 模块处理固定有效损失速度、电子能量壁损失和固定产额二次电子发射；吸附态反应、
  表面化学网络及依赖表面状态的发射仍未实现。
