# 网格、边界与程序组织

## 从 Gmsh 到有限体积网格

`mesh::MoabMesh` 通过 MOAB 读取网格并一次性预计算以下面向热路径的数据：

- 稠密的 `CellId`、`FaceId`，以及 face 的 owner/neighbor；
- CSR 风格的单元—面邻接；
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
mesh ──► field ───────────────┐
  │        boundary ──────────┼──► discretization ──► equation
  └───────────────────────────┘          │                 │
                                         linalg ◄───────────┘
```

- `mesh`：几何、拓扑与边界物理组；
- `field`：单元/面标量容器；
- `boundary`：Dirichlet、Neumann 及其集合；
- `discretization`：扩散、散度、迎风通量、泊松装配和后向欧拉装配；
- `linalg`：Eigen 稀疏矩阵类型与求解器统一接口；
- `equation`：将装配与求解流程封装成稳态 `PoissonSolver`、瞬态
  `TransientDiffusionSolver`。

线性后端分两阶段工作：`analyzePattern(A)` 分析稀疏结构，`factorize(A)` 进行数值
分解，最后 `solve(b,x)`。CHOLMOD 面向对称正定矩阵；UMFPACK 可处理一般非对称
稀疏矩阵。方程求解器首次调用时装配矩阵并完成前两阶段，后续仅重组右端项并复用
分解；这适合网格、系数、边界和时间步不变的多次求解。

## 实现边界

稳态 `PoissonSolver` 的分解复用意味着如需改变扩散系数、网格或边界条件，应调用
`reset()` 或新建求解器。瞬态求解器同样假定每一步的矩阵不变，只有状态相关的 RHS
改变。纯 Neumann 的稳态泊松问题存在常数零空间；项目没有为其额外施加均值约束，
使用时应至少施加足够的 Dirichlet 约束或自行处理相容性与规范条件。
