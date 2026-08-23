# 网格、边界与程序组织

## 从 Gmsh 到有限体积网格

`mesh::MoabMesh` 通过 MOAB 读取网格并一次性预计算以下面向热路径的数据：

- 稠密的 `CellId`、`FaceId`、`VertexId`，以及 face 的 owner/neighbor；
- CSR 风格的单元—面、单元—顶点邻接与顶点坐标；
- 单元中心与面积、面中心与长度、单位面法向；
- Gmsh physical group 到 `BoundaryId` 的映射。

原始 MOAB handle 只服务于构造阶段；求解阶段按稠密整数 ID 访问。公共接口
`IMesh` 抽象了这些量，计算模块不依赖 MOAB。二维网格中多边形面积和重心由鞋带
公式计算；面法向会被翻转为 owner 到 neighbor（或朝外）的约定。

## 数据与边界

`CellField<T>` 和 `FaceField<T>` 分别拥有与某一个 `IMesh` 绑定的连续存储。
算子在组合场前检查网格对象是否相同，以避免“编号碰巧相同但属于不同网格”的
错误。`BoundaryConditionSet` 是 `BoundaryId → variant<Dirichlet, Neumann>` 的
映射；重复写入同一 ID 会替换旧条件，读取未定义边界会抛出异常。

## 模块关系

```text
unit ─────────► field ◄──────── mesh
                 ▲               │
                 │               └────► boundary
physics ─────────┘

trace ─────► linalg ────────────────┐
                                    │
field + boundary + mesh + physics ──┴─► discretization ─► equation ─► simulation
                                                                        ▲
                                                                        │
                                                                      trace
```

- `mesh`：几何、拓扑与边界物理组；
- `unit`：mp-units 编译期真相源、到 LLNL `precise_unit` 的无字符串 bridge、项目级 quantity specification 及运行时单位元数据；
- `field`：单元/面标量容器、强类型 ID 场集合及场单位元数据；
- `trace`：结构化执行事件、sink concept、无领域依赖的物理测度/标量统计器与通用同步 sink；
- `output`：以 `output::dump` 和 `output::checkpoint` 分隔 ParaView 时序输出与版本化 VTKHDF 状态恢复，共享 `IMesh` 到 VTK 网格适配；它读取 `mesh/field`，只向 trace sink 发出轻量完成事件；
- `boundary`：Dirichlet、Neumann 及其集合；
- `discretization`：扩散、散度、迎风通量、泊松装配和后向欧拉装配；
- `linalg`：Eigen 稀疏矩阵类型与求解器统一接口；公共阶段统一返回可附带根因 diagnostic
  的 `SolverResult`，但不直接输出；
- `equation`：将装配与求解流程封装成稳态 `PoissonSolver`、瞬态
  `TransientDiffusionSolver`。
- `simulation`：组织电静力、反应、输运和时钟，并向注入的 trace sink 报告阶段状态。

线性后端分两阶段工作：`analyzePattern(A)` 分析稀疏结构，`factorize(A)` 进行数值
分解，最后 `solve(b,x)`。CHOLMOD 面向对称正定矩阵；UMFPACK 可处理一般非对称
稀疏矩阵。方程求解器首次调用时装配矩阵并完成前两阶段，后续仅重组右端项并复用
分解；这适合网格、系数、边界和时间步不变的多次求解。三个阶段均返回
`SolverResult`：成功时无 diagnostic，API 边界或后端失败时由 linalg 形成结构化根因，
equation 原样传播，simulation 再交给其 sink。

## 实现边界

稳态 `PoissonSolver` 的分解复用意味着如需改变扩散系数、网格或边界条件，应调用
`reset()` 或新建求解器。瞬态求解器同样假定每一步的矩阵不变，只有状态相关的 RHS
改变。

纯 Neumann 的稳态泊松问题存在常数零空间，构造 `PoissonSolver` 时必须显式传入
`PureNeumannOptions`。可选择 `PinCellGauge` 固定一个参考单元，或选择
`ZeroMeanGauge` 施加体积加权零均值约束。求解器会在每次重组 RHS 后检查离散相容性
`sum(b) = 0`；不相容时返回 `SolverStatus::IncompatibleRhs`。`PinCellGauge` 保持系统
对称正定，可使用 CHOLMOD；`ZeroMeanGauge` 会增加一个拉格朗日乘子并形成对称不定
系统，应使用 UMFPACK 等支持一般稀疏矩阵的后端。

两种策略可分别按以下方式传给直接求解器；漂移扩散 stepper 的最后一个构造参数也
接受相同配置：

```cpp
PureNeumannOptions pin{
    .gauge = PinCellGauge{.cell = reference_cell, .value = reference_phi}};

PureNeumannOptions zero_mean{.gauge = ZeroMeanGauge{}};
```

规范条件目前按单个连通计算域设计；若网格包含多个互不连通的区域，每个连通分量都
需要独立规范条件，当前接口不会自动添加这些额外约束。

物种与电子能量的 Scharfetter--Gummel 边界同时支持 Dirichlet 和 Neumann。这里的
Neumann `value` 表示外向扩散通量 $-D\nabla n\cdot\mathbf n$，总法向通量仍包含
$v_nn_P$；因此齐次 Neumann 只有在 $v_n=0$ 时才是严格 zero-flux。M15 线性 wall law
会覆盖同一物种、同一边界上的普通条件，不能与普通 Neumann 通量叠加。

## 轴对称度量视图

`AxisymmetricMeshView` 将普通二维网格解释为 $(r,z)$ 子午面，其中 `x=r\ge0`、`y=z`，
绕 $r=0$ 旋转后供现有有限体积算子直接使用。它保留拓扑、坐标、单元/面中心和子午面
单位法向，只在构造时预计算物理度量：

$$
V_P=2\pi r_P A_P,\qquad A_f^{\mathrm{axi}}=2\pi r_f L_f.
$$

这里使用 Pappus 质心定理；对于直边网格，公式也分别给出环形控制体体积与旋转面面积。
轴线面 $r_f=0$ 的物理面积自然为零。该类是非拥有视图，被包装的二维 mesh 必须比视图
活得更久，且构造后不能改变几何。它会拒绝非二维、非有限、负半径或非正底层度量，避免
把普通 XY 网格静默当作轴对称网格。
