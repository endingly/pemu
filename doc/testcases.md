# 测试用例契约

本页逐一说明仓库中所有 91 个 GoogleTest 用例所守护的契约。数值阈值并非一般性
精度承诺，而是当前测试网格、双精度实现和制造解下的回归界限。`two_quads.msh`
包含两个相邻单位方形；大多数几何与算子测试以它为夹具。

## 边界条件（`submod/boundary/test/test_boundary.cpp`）

| 用例 | 保证 |
| --- | --- |
| `StoresDirichlet` | 可按 ID 保存并无歧义地取回 Dirichlet 值。 |
| `StoresNeumann` | 可保存并取回带符号的外法向 Neumann 通量。 |
| `ReplacesExistingCondition` | 同一边界 ID 的新条件覆盖旧条件，避免旧物理条件残留。 |
| `MissingBoundaryThrows` | 未配置边界不能被静默读取。 |
| `ReportsSize` | 条件集合的计数反映不同边界 ID 数。 |
| `EveryBoundaryFaceHasCondition` | 两单元网格的六个边界面都能由其 physical group 找到条件。 |
| `DetectsMissingBoundaryCondition` | 漏配第 4 号物理边界可被遍历检查发现。 |

## 网格（`submod/mesh/test/test_moab_mesh.cpp`）

| 用例 | 保证 |
| --- | --- |
| `LoadsTwoDimensionalMesh` | 输入被识别为二维网格。 |
| `HasExpectedEntityCounts` | 两四边形夹具有 2 单元、6 顶点、7 条唯一面。 |
| `HasSixBoundaryFacesAndOneInternalFace` | 边界/内面分类及边界面无 neighbor 的拓扑语义正确。 |
| `EachCellHasFourFaces` | 每个四边形连接四条面。 |
| `CellFaceConnectivityIsConsistent` | 每条面同时出现在其 owner（及内面的 neighbor）邻接表中。 |
| `CellsShareExactlyOneFace` | 两控制体只有一个公共内面。 |
| `ComputesCorrectCellCenters` | 两个中心为 `(0.5,0.5)`、`(1.5,0.5)`，且不假设 ID 顺序。 |
| `ComputesCorrectCellAreas` | 二维 `cellVolume` 是正确的单位面积。 |
| `ComputesCorrectFaceLengths` | 二维 `faceArea` 是正确的单位边长。 |
| `FaceNormalsAreUnitVectors` | 所有面法向均归一化。 |
| `InternalFaceNormalPointsFromOwnerToNeighbor` | 内面方向满足 `n·(x_N-x_P)>0`。 |
| `BoundaryFaceNormalsPointOutward` | 边界法向从控制体中心朝区域外。 |
| `InternalFaceHasCorrectGeometry` | 公共内面的中心、长度和法向与几何一致。 |
| `UsesDenseValidCellIds` | 暴露给数值模块的 cell ID 稠密且全部合法。 |
| `RejectsInvalidFaceId` | 所有面查询都会拒绝越界 ID。 |
| `RejectsInvalidCellId` | 所有单元查询都会拒绝越界 ID。 |
| `ReadsLeftBoundaryPhysicalGroup` | 左侧面读取物理组 1。 |
| `ReadsRightBoundaryPhysicalGroup` | 右侧面读取物理组 2。 |
| `ReadsBottomBoundaryPhysicalGroup` | 底部两面读取物理组 3。 |
| `ReadsTopBoundaryPhysicalGroup` | 顶部两面读取物理组 4。 |
| `InternalFaceHasNoBoundaryId` | 内面返回 `invalid_boundary`，不会伪装成边界。 |

## 场容器（`submod/field/test/test_field.cpp`）

| 用例 | 保证 |
| --- | --- |
| `CellFieldSizeMatchesMesh` | 单元场长度等于单元数。 |
| `CellFieldSupportsInitialValue` | 单元场可用常数初值构造。 |
| `CellFieldSupportsDenseIdAccess` | 稠密 `CellId` 可读写对应分量。 |
| `CellFieldFillWorks` | `fill` 覆盖全部单元自由度。 |
| `CellFieldAtRejectsInvalidCell` | 检查式单元访问会拒绝越界。 |
| `CellFieldKeepsMeshAssociation` | 单元场保留原网格身份。 |
| `FaceFieldSizeMatchesMesh` | 面场长度等于面数。 |
| `FaceFieldSupportsInitialValue` | 面场可用常数初值构造。 |
| `FaceFieldSupportsDenseIdAccess` | 稠密 `FaceId` 可读写对应分量。 |
| `FaceFieldAtRejectsInvalidFace` | 检查式面访问会拒绝越界。 |
| `FaceFieldKeepsMeshAssociation` | 面场保留原网格身份。 |

## 离散算子（`submod/discretization/test/test_operator.cpp`）

| 用例 | 保证 |
| --- | --- |
| `DivergenceConservesInternalFlux` | 单个内面通量在全域积分散度中精确相消。 |
| `DivergenceUsesOwnerNeighborOrientation` | 正法向通量对 owner 为正、对 neighbor 为负。 |
| `DivergenceSatisfiesDiscreteGaussTheorem` | 体积积分散度等于边界积分通量（此夹具为 6）。 |
| `DiffusionFluxIsExactForLinearField` | 线性场 `u=x` 的离散扩散通量精确为 `-D n_x`。 |
| `DivergenceOfLinearDiffusionFluxIsZero` | 常梯度扩散场的离散散度为零。 |
| `DivergenceRejectsDifferentMeshes` | 散度的输入/输出场不能跨网格混用。 |
| `UpwindFluxUsesOwnerForPositiveVelocity` | 内面正速度取 owner 上游值。 |
| `UpwindFluxUsesNeighborForNegativeVelocity` | 内面负速度取 neighbor 上游值。 |
| `UpwindFluxUsesInteriorStateForOutflow` | 边界出流无需外部条件，使用内部状态。 |
| `UpwindFluxUsesBoundaryStateForInflow` | 边界入流使用 Dirichlet 外部状态。 |
| `UpwindFluxRejectsMissingInflowBoundaryValue` | 入流缺少状态条件会报错。 |
| `UpwindFluxRejectsNeumannOnInflow` | Neumann 不能替代对流入流所需的状态值。 |
| `ConstantStateHasZeroAdvectionDivergence` | 常值状态与相容入流在常速度下无虚假散度。 |
| `AdvectionFluxIsGloballyConservative` | 纯内面对流的全域积分散度为零。 |
| `UpwindFluxRejectsDifferentMeshes` | 迎风算子的状态、速度、通量必须属于同一网格。 |

## 泊松有限体积装配（`submod/discretization/test/test_poisson_fvm*.cpp`）

| 用例 | 保证 |
| --- | --- |
| `AssemblesExpectedMatrixForTwoQuads` | 单位扩散、零 Dirichlet 下产生 `[[7,-1],[-1,7]]`。 |
| `AssemblesVolumeSource` | 源项按 `ρ_P V_P` 进入 RHS。 |
| `AppliesNonZeroDirichletBoundary` | 左边界值 1 以正确的边界传导系数进入左单元 RHS。 |
| `AppliesNeumannBoundaryFlux` | 正外向通量按 `-qA` 进入相邻单元 RHS。 |
| `MatrixIsSymmetric` | 混合边界下扩散离散矩阵仍为对称矩阵。 |
| `IntergrationTest` | 装配的混合边界系统可被 CHOLMOD 分析、分解和求解。 |
| `MeshFixturesHaveExpectedCellCounts` | 8/16/32/64 网格分别有 `n²` 单元。 |
| `MeshFixturesHaveExpectedBoundaryFaceCounts` | 方形 `n×n` 网格分别有 `4n` 个边界面。 |
| `HasSecondOrderSpatialConvergence` | 制造解的 L2 观测空间阶趋近 2（最细层 `2±0.05`）。 |
| `LinearSolverResidualIsSmall` | 32² 制造解的相对代数残差小于 `1e-10`。 |
| `ErrorDecreasesUnderMeshRefinement` | 逐级加密时 L2 误差单调下降。 |
| `ErrorReducesApproximatelyByFactorFour` | 二阶主导时加密一倍的误差比约为 4。 |

## 方程求解器（`submod/equation/test/test_poisson_solver.cpp`）

| 用例 | 保证 |
| --- | --- |
| `SolvesConstantDirichletSolution` | 零源和全边界常值 1 的离散解为常值 1。 |
| `WritesSolutionDirectlyIntoField` | 解会原位写回调用者给出的 `CellField`。 |
| `RejectsSolutionFieldFromDifferentMesh` | 解场与离散网格不同时拒绝求解。 |
| `RejectsNullLinearSolver` | 方程求解器不能在没有线性后端时构造。 |
| `ReusesFactorizationAcrossSolves` | 两次相同矩阵求解只分析/分解一次，但求解两次。 |

## 瞬态扩散（`submod/equation/test/test_transient_diffusion_solver.cpp`）

| 用例 | 保证 |
| --- | --- |
| `ConstantSolutionRemainsConstant` | 与边界相容的常值是后向欧拉扩散的不变解。 |
| `ZeroNeumannBoundaryConservesMass` | 零通量边界下 20 步后 `Σu_PV_P` 守恒。 |
| `DiffusionReducesCellDifference` | 扩散平滑两单元的浓度差，同时守恒总量。 |
| `HasSecondOrderSpatialConvergence` | 极小时间步隔离空间误差后，L2 空间阶约为 2。 |
| `HasFirstOrderTemporalConvergence` | 相同终止时间上，步长逐半的时间收敛阶约为 1。 |

## 线性代数（`submod/linalg/test/test_cholmod_solver.cpp`）

| 用例 | 保证 |
| --- | --- |
| `CholmodSolverTest.SolvesSpdSystem` | CHOLMOD 能解已知对称正定稀疏系统，解和相对残差均正确。 |
| `CholmodSolverTest.FactorizeWithoutAnalyzeFails` | 未分析结构时禁止数值分解。 |
| `CholmodSolverTest.SolveWithoutFactorizationFails` | 未分解时禁止求解。 |
| `CholmodSolverTest.ReusesPatternForNewMatrixValues` | 结构不变、数值改变时可复用符号分析并得到正确解。 |
| `CholmodSolverTest.ResetClearsSolverState` | `reset` 清空已分析和已分解状态。 |
| `UmfpackSolverTest.SolvesGeneralSparseSystem` | UMFPACK 能解已知一般非对称稀疏系统。 |
| `UmfpackSolverTest.FactorizeWithoutAnalyzeFails` | UMFPACK 同样强制分析先于分解。 |
| `UmfpackSolverTest.SolveWithoutFactorizationFails` | UMFPACK 同样强制分解先于求解。 |
| `UmfpackSolverTest.ReusesPatternForNewMatrixValues` | UMFPACK 在固定非对称模式下可复用符号分析。 |
| `UmfpackSolverTest.ResetClearsSolverState` | UMFPACK `reset` 清空生命周期状态。 |

## 独立五点差分泊松验证（`submod/linalg/test/test_possion.cpp`）

这组测试不经 FVM 网格模块，而是在单位正方形内点上直接组装五点差分，用于同时
验证 CHOLMOD 和经典二阶拉普拉斯离散。

| 用例 | 保证 |
| --- | --- |
| `SolvesManufacturedSolution` | `sin(πx)sin(πy)` 制造解的残差 `<1e-10`、L2 误差 `<2e-3`。 |
| `HasSecondOrderSpatialConvergence` | 8→64 内点网格的观测阶明确趋近二阶。 |
| `MatrixIsSymmetric` | 五点拉普拉斯矩阵严格对称。 |
| `MatrixCanBeCholeskyFactorized` | 该 Dirichlet 离散矩阵为可 CHOLMOD 分解的 SPD 矩阵。 |
| `HasExpectedSparsityPattern` | `n×n` 内点的非零元数为 `5n²-4n`。 |
