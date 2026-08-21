# 测试用例契约

本页逐一说明仓库中所有 206 个 GoogleTest 用例所守护的契约。数值阈值并非一般性
精度承诺，而是当前测试网格、双精度实现和制造解下的回归界限。`two_quads.msh`
包含两个相邻单位方形；大多数几何与算子测试以它为夹具。

## 日志与追踪（`submod/trace/test/test_trace.cpp`）

| 用例 | 保证 |
| --- | --- |
| `TraceSinkTest.NullSinkAcceptsStructuredEvents` | 默认空 sink 满足 trace concept，能同步接受结构化事件且不产生副作用。 |
| `PhysicalVolumeSemanticsTest.Extrudes2DMeasuresAndKeeps3DMeasures` | 二维网格的单元面积/面长度乘以面外厚度得到物理体积/面积，而三维网格保持原生几何测度。 |
| `StatisticsOptionsTest.SamplesOnlyAtConfiguredStepInterval` | 统计探针只在启用且步号满足采样间隔时运行，并拒绝零间隔等非法配置。 |
| `ScalarFieldStatisticsTest.ComputesWeightedExtremaIntegralMeanAndRms` | 标量统计器在非均匀物理权重下正确计算极值、最大绝对值、积分、L1 积分、加权均值和 RMS。 |
| `ScalarFieldStatisticsTest.DerivesIntegralUnitOutsideSampleLoop` | 统计器只在收尾边界附加并组合运行期单位，验证 $\mathrm{C/cm^3}\times\mathrm{cm^3}=\mathrm C$，而逐样本 `add` 仍只接收两个 `double`。 |
| `ScalarFieldStatisticsTest.ReportsNonFiniteValuesAndInvalidWeights` | 标量统计器显式计数负值、非有限值和非法物理权重，不把异常样本静默混入积分。 |
| `TraceSinkTest.OstreamSinkFormatsTabularStructuredEvents` | 普通事件 renderer 保持原有统一事件表格式，将步数、时间、步长、终止时间、求解状态和残差放入稳定对齐列，其余属性进入 DETAILS。 |
| `TraceSinkTest.StatisticsRendererGroupsRowsAndFormatsDedicatedColumns` | 统计 renderer 按 step 缓冲分组，以 FIELD/UNIT/MIN/MAX/MEAN/RMS/INTEGRAL 独立表格展示多个场；固定采样与体积语义元数据只输出一次，单位来自事件属性而非字段名猜测，零异常计数和 DETAILS 列均不出现。 |
| `TraceSinkTest.StatisticsRendererOnlyShowsNonZeroCountsAsDiagnostics` | 统计质量计数仅在非零时进入该 step 的 `STATISTICS DIAGNOSTICS` 小节，并保留 Warning/Diagnostic 级别；值为零的 `non_finite` 不输出。 |
| `TraceSinkTest.OstreamSinkFlushesEveryTwoSimulationStepsAndOnRunEnd` | ostream sink 在每两个 simulation 完成步后刷新文件流，并在正常或失败结束时刷新不足两个时间步的尾部记录。 |
| `TraceSinkTest.SplitSinkRoutesStatisticsSeparatelyAndFlushesThemEveryTwoSteps` | 双通道 sink 仅把显式标记的统计事件写入统计流，把普通 trace/diagnostic 写入诊断流；即使统计事件自身是 Warning/Diagnostic 也不会串流，并且统计流每两个完成步刷新一次。 |
| `EnumStringTest.DomainsAndSeveritiesHaveStableNames` | `DiagDomain` 与 `Severity` 均通过 `pemu::to_string` 映射为稳定、可读的字符串，可供过滤和机器处理使用。 |
| `TraceSinkTest.NullSinkAcceptsStructuredDiagnostics` | 同一个默认空 trace sink 能无副作用地接受标记为 diagnostic 的结构化事件。 |
| `TraceSinkTest.OstreamSinkFormatsDiagnosticEvent` | diagnostic renderer 保留根因消息与数值上下文，但与前后的普通事件共用同一个旧式事件表头；普通→diagnostic→普通切换不会插入 TRACE/DIAGNOSTIC 分段标题。 |

## 单位（`submod/unit/test/test_unit.cpp`）

| 用例 | 保证 |
| --- | --- |
| `PlasmaQuantitiesTest.DefineDimensionallyConsistentReferences` | 项目定义的数密度、数密度变化率、反应率密度、法向电场和法向漂移速度 quantity specification 只接受量纲相容的 mp-units unit。 |
| `MpUnitsBridgeTest.ProducesRuntimeMetadataWithoutStringParsing` | mp-units reference 在编译期直接桥接为可平凡复制的 `QuantityKind + precise_unit`；数密度、电荷密度、电势和法向电场的运行期单位正确，且不依赖字符串解析。 |

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
| `ExposesOrderedCellVertexConnectivity` | `IMesh` 暴露构造非结构网格所需的有序 cell→vertex 稠密连通关系。 |
| `ExposesDenseVertexCoordinates` | 每个稠密 `VertexId` 都能返回原始顶点坐标，且两单元夹具的包围盒正确。 |
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
| `FieldTest.FieldSetGroupsTypedFieldsAndPropagatesMetadata` | `FieldSet` 能按强类型稠密 ID 聚合相互独立的单元场与面场，并统一保留网格、初值、元数据、遍历和批量填充语义。 |
| `FieldTest.PhysicalQuantityIsMetadataAndRawStorageRemainsDouble` | metadata 只保存由 mp-units bridge 产生的 `QuantityKind + precise_unit`，该表示可平凡复制，而场元素、`data()` 和核心存储仍为 `double`。 |
| `FieldTest.MpUnitsBridgeCoversCanonicalPlasmaFieldMetadata` | 编译期 bridge 不经字符串解析即可把数密度、电荷密度和法向电场 reference 映射为正确的语义枚举与 LLNL `precise_unit`。 |
| `FieldTest.QuantityBoundaryConvertsToFieldStorageUnit` | `fillQuantity`/`setQuantity` 可将米输入换算为按厘米存储的裸值，`quantityAt` 能重新包装和换算输出。 |
| `FieldTest.QuantityBoundaryRejectsReferenceDifferentFromMetadata` | 边界适配器会拒绝单位不同的 reference，也会拒绝同为伏特但 quantity specification 不同的 reference。 |

## VTKHDF 场输出（`submod/output/test/test_vtkhdf_writer.cpp`）

| 用例 | 保证 |
| --- | --- |
| `MeshAdapterTest.ConvertsIMeshToPolygonalUnstructuredGrid` | `IMesh` 的顶点、cell→vertex 连通关系被无损映射为 `vtkUnstructuredGrid` polygon cells。 |
| `ParaViewReadabilityTest.OfficialVtkHdfReaderRoundTripsFieldsMetadataAndTime` | 官方 `vtkHDFWriter` 产物能由 ParaView 使用的官方 `vtkHDFReader` 读回；CellData、原始 face FieldData/拓扑、面积加权 cell-centered 可视化副本、field name、运行期 quantity/unit 以及 step/time 均保持一致。 |
| `OutputTraceTest.EmitsOnlyLightweightCompletionContextAfterWrite` | 成功输出后只向 trace 发出 `output.completed` 及 `path/step/time` 三个轻量属性，不把场数组塞入 trace。 |
| `VtkHdfWriterTest.RefusesOverwriteUnlessExplicitlyEnabled` | writer 默认拒绝覆盖已有文件，只有请求显式设置 `overwrite` 时才覆盖。 |

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
| `BernoulliTest.ValueAtZeroIsOne` | Scharfetter–Gummel 离散使用的 Bernoulli 函数满足 $B(0)=1$。 |
| `BernoulliTest.MatchesKnownValues` | $B(1)$ 与 $B(-1)$ 符合高精度参考值。 |
| `BernoulliTest.SatisfiesDifferenceIdentity` | 在正负及不同量级参数上满足 $B(-x)-B(x)=x$。 |
| `BernoulliTest.IsStableNearZero` | $x\to0$ 时采用稳定计算，并符合 $B(x)=1-x/2+O(x^2)$。 |
| `BernoulliTest.HasCorrectLargeArgumentLimits` | 大正参数下 $B(x)\to0$，大负参数下 $B(x)\sim-x$，且计算不溢出。 |
| `OperatorTest.ScharfetterGummelReducesToDiffusionAtZeroVelocity` | 速度为零时，Scharfetter–Gummel 内面通量退化为中心扩散通量。 |
| `OperatorTest.ScharfetterGummelApproachesOwnerUpwindForStrongPositiveDrift` | 强正漂移极限取 owner 状态，趋近正向迎风通量。 |
| `OperatorTest.ScharfetterGummelApproachesNeighborUpwindForStrongNegativeDrift` | 强负漂移极限取 neighbor 状态，趋近负向迎风通量。 |
| `OperatorTest.ScharfetterGummelPreservesConstantStateFlux` | 常值状态与相容边界下仅保留物理对流通量 $v_nu$。 |
| `OperatorTest.ConstantSgFluxHasZeroDivergence` | 常值状态的 Scharfetter–Gummel 通量不会产生虚假离散散度。 |
| `OperatorTest.ScharfetterGummelExactlyPreservesExponentialEquilibrium` | 对满足相邻单元指数平衡关系的状态，离散漂移扩散通量精确为零。 |

## 泊松有限体积装配（`submod/discretization/test/test_poisson_fvm*.cpp`）

| 用例 | 保证 |
| --- | --- |
| `AssemblesExpectedMatrixForTwoQuads` | 单位扩散、零 Dirichlet 下产生 `[[7,-1],[-1,7]]`。 |
| `AssemblesVolumeSource` | 源项按 `ρ_P V_P` 进入 RHS。 |
| `AppliesNonZeroDirichletBoundary` | 左边界值 1 以正确的边界传导系数进入左单元 RHS。 |
| `AppliesNeumannBoundaryFlux` | 正外向通量按 `-qA` 进入相邻单元 RHS。 |
| `DetectsPureNeumannBoundary` | 全部边界均为 Neumann 时不报告 Dirichlet 面。 |
| `DetectsMixedBoundaryAsHavingDirichlet` | 混合边界中能够识别至少一个 Dirichlet 面。 |
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
| `PureNeumannRequiresExplicitGauge` | 纯 Neumann 问题未指定规范条件时拒绝构造。 |
| `SolvesPureNeumannWithPinnedCell` | 固定参考单元后恢复线性解析解。 |
| `SolvesPureNeumannWithZeroMean` | 拉格朗日乘子约束得到体积加权零均值线性解。 |
| `RejectsIncompatiblePureNeumannRhs` | 总源项与边界通量不平衡时返回 `IncompatibleRhs`，附带 equation/poisson 根因 diagnostic，且不覆盖原解场。 |

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
| `SolverStatusTest.NamesAreStableAndHumanReadable` | 每个线性求解状态都有稳定、可读的字符串名，trace 不依赖枚举底层整数。 |
| `CholmodSolverTest.SolvesSpdSystem` | CHOLMOD 能解已知对称正定稀疏系统，解和相对残差均正确。 |
| `CholmodSolverTest.FactorizeWithoutAnalyzeFails` | 未分析结构时禁止数值分解，并返回带 linalg 域、CHOLMOD 类别和稳定失败名的 diagnostic。 |
| `CholmodSolverTest.SolveWithoutFactorizationFails` | 未分解时禁止求解，并在 API 失败出口附带 diagnostic。 |
| `CholmodSolverTest.RejectsMismatchedSolveDimensionsWithDiagnostic` | RHS 与解向量尺寸不一致时不进入后端求解，并返回 `InvalidInput` 及结构化 diagnostic。 |
| `CholmodSolverTest.ReusesPatternForNewMatrixValues` | 结构不变、数值改变时可复用符号分析并得到正确解。 |
| `CholmodSolverTest.ResetClearsSolverState` | `reset` 清空已分析和已分解状态。 |
| `UmfpackSolverTest.SolvesGeneralSparseSystem` | UMFPACK 能解已知一般非对称稀疏系统。 |
| `UmfpackSolverTest.FactorizeWithoutAnalyzeFails` | UMFPACK 同样强制分析先于分解，并附带自身类别的 diagnostic。 |
| `UmfpackSolverTest.SolveWithoutFactorizationFails` | UMFPACK 同样强制分解先于求解，并附带根因 diagnostic。 |
| `UmfpackSolverTest.RejectsMismatchedSolveDimensionsWithDiagnostic` | UMFPACK 在 API 边界拒绝不一致的 RHS/解尺寸并返回 diagnostic，不进入核心求解。 |
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

## 电场与带电粒子漂移（`submod/discretization/test/test_electric_field.cpp` 等）

| 用例 | 保证 |
| --- | --- |
| `DriftVelocityTest.PositiveSpeciesDriftsAlongElectricField` | 正电荷漂移速度与电场同向。 |
| `DriftVelocityTest.NegativeSpeciesDriftsAgainstElectricField` | 负电荷漂移速度与电场反向。 |
| `ElectrostaticDriftTest.PotentialProducesOppositeIonAndElectronDrift` | 同一电势梯度使电子与正离子产生相反漂移。 |
| `ElectricFieldTest.LinearPotentialProducesExactNormalElectricField` | 线性电势的离散面法向电场与解析值一致。 |

## 物理模型（`submod/physics/test`）

| 用例 | 保证 |
| --- | --- |
| `ChargeDensityTest.EqualOppositeSpeciesAreNeutral` | 等密度、等电荷量异号的两种粒子产生零净电荷密度。 |
| `ReactionTest.ElectronImpactIonizationComputesExpectedRate` | 电子碰撞电离率按给定电子密度、中性粒子密度和速率系数计算。 |
| `ReactionTest.ZeroElectronDensityProducesNoIonization` | 没有电子时电子碰撞电离率严格为零。 |
| `ReactionTest.IonizationCreatesElectronIonPairs` | 一次电离同时产生一个电子和一个正离子。 |
| `ReactionTest.PairProductionAccumulatesIntoExistingSource` | 反应源项累加而非覆盖已有源项。 |
| `ReactionTest.PairProductionCreatesNoNetCharge` | 成对产生电子和离子不会凭空产生净电荷。 |
| `ReactionNetworkTest.PairIonizationConservesCharge` | 反应网络能识别成对电离的电荷守恒。 |
| `ReactionNetworkTest.DetectsChargeViolatingReaction` | 违反电荷守恒的化学计量关系会被检测。 |
| `ReactionNetworkTest.AccumulatesStoichiometricSources` | 反应率依据化学计量数正确转换为各物种源项。 |
| `ReactionNetworkTest.MultipleReactionsAccumulateCorrectly` | 多个反应对同一物种的贡献正确求和。 |
| `SpeciesSetTest.AssignsDenseStableIds` | 物种 ID 稠密、稳定并可用于场数组索引。 |
| `SpeciesFieldsTest.StoresIndependentFieldsPerSpecies` | 每个物种拥有互不串扰的独立场。 |

## 固定步长与自适应步长方程推进器（`submod/equation/test`）

| 用例 | 保证 |
| --- | --- |
| `FixedStepElectrostaticDriftDiffusionStepperTest.UniformNeutralPlasmaRemainsStationary` | 固定步长下，均匀电中性等离子体是静止离散解。 |
| `FixedStepElectrostaticDriftDiffusionStepperTest.AppliedPotentialProducesOppositeSpeciesDrift` | 外加电势使电子与离子按电荷极性反向漂移。 |
| `FixedStepElectrostaticDriftDiffusionStepperTest.CflFailureDoesNotPartiallyAdvanceSpecies` | 任一物种违反 CFL 时，所有物种都保持更新前状态。 |
| `FixedStepElectrostaticDriftDiffusionStepperTest.RejectsInvalidSpeciesChargePolarity` | 电子、离子的电荷极性配置错误会在构造时被拒绝。 |
| `FixedStepElectrostaticDriftDiffusionReactionTest.IonizationProducesNeutralElectronIonPairs` | 固定步长漂移扩散与电离耦合后仍成对产生电子、离子并保持电中性。 |
| `FixedStepExplicitEulerAdvectionStepperTest.ComputesExpectedCfl` | 固定步长迎风推进器给出与网格、速度和步长一致的 CFL 数。 |
| `FixedStepExplicitEulerAdvectionStepperTest.RejectsCflGreaterThanOne` | 显式迎风 CFL 大于 1 时拒绝推进。 |
| `FixedStepExplicitEulerAdvectionStepperTest.CflOneMovesStateAcrossInternalFace` | CFL 等于 1 时，状态恰好跨过内面且总量守恒。 |
| `FixedStepExplicitSpeciesContinuityStepperTest.ConstantDensityRemainsConstant` | 与边界相容的常密度在固定步长 SG 输运下不变。 |
| `FixedStepExplicitSpeciesContinuityStepperTest.SourceIncreasesDensity` | 常源项按照 $n^{k+1}=n^k+\Delta tS$ 增加密度。 |
| `FixedStepExplicitSpeciesContinuityStepperTest.SourceChangesTotalParticlesByIntegratedSource` | 总粒子数变化等于 $\Delta t$ 乘体积积分源项。 |
| `FixedStepMultiSpeciesDriftDiffusionStepperTest.UniformNeutralPlasmaRemainsStationary` | 固定多物种推进器保持均匀中性状态。 |
| `FixedStepMultiSpeciesDriftDiffusionStepperTest.ChargeDensitySupportsMoreThanTwoSpecies` | 电荷密度对任意数量物种按 $\rho=\sum_s q_sn_s$ 求和。 |
| `FixedStepMultiSpeciesDriftDiffusionStepperTest.ImmobileSpeciesIsNotUpdated` | 标为 immobile 的物种不参与输运更新。 |
| `FixedStepMultiSpeciesDriftDiffusionStepperTest.AppliedPotentialProducesCorrectDriftForAllChargedSpecies` | 每个带电物种按迁移率和极性得到正确漂移速度。 |
| `FixedStepMultiSpeciesDriftDiffusionStepperTest.ComputesDriftForThreeTransportedChargedSpecies` | 漂移计算不隐含“只有电子和一种离子”的限制。 |
| `FixedStepMultiSpeciesDriftDiffusionStepperTest.CflFailureDoesNotPartiallyUpdateSpecies` | 多物种 CFL 预检具有原子性，不留下半更新状态。 |
| `FixedStepMultiSpeciesDriftDiffusionStepperTest.RejectsWrongSpeciesFieldCount` | 场数组物种数必须与 `SpeciesSet` 一致。 |
| `FixedStepMultiSpeciesDriftDiffusionStepperTest.RejectsFieldsFromDifferentMesh` | 多物种场不能跨网格传给推进器。 |
| `FixedStepMultiSpeciesDriftDiffusionStepperTest.MatchesLegacyTwoSpeciesSolver` | 更一般的固定多物种实现与固定双物种实现给出相同结果。 |
| `FixedStepMultiSpeciesDriftDiffusionStepperTest.PureNeumannGaugeIsAvailableFromCoupledStepper` | 纯 Neumann 规范配置可由多物种漂移扩散推进器传入内部泊松求解器，并得到零电势、零电场的中性均匀解。 |
| `AdaptiveTimeStepControllerTest.AppliesSafetyAndGrowthLimits` | 时间步提议同时服从安全系数、稳定性限制和相对上一时刻的增长上限。 |
| `AdaptiveTimeStepControllerTest.TruncatesFinalStepToRemainingTime` | 时间步提议会被剩余时间截断，保证最终时刻精确对齐。 |
| `AdaptiveTimeStepControllerTest.RejectsNoPositiveAdmissibleTimeStep` | 稳定性或正性限制为零时拒绝不存在正步长的推进。 |
| `AdaptiveStepMultiSpeciesDriftDiffusionStepperTest.UsesConfiguredMaximumTimeStep` | 稳定性允许更大步长时，自适应推进器仍服从用户给定的 `max_dt`。 |
| `AdaptiveStepMultiSpeciesDriftDiffusionStepperTest.FinalTimeStepEqualsRemainingTime` | 最后一步被截断为剩余时间，避免越过终止时刻。 |

## 固定步长与自适应步长仿真（`submod/simulation/test/test_simulation.cpp`）

| 用例 | 保证 |
| --- | --- |
| `FixedStepClockTest.AdvancesTimeFromStepIndex` | 固定时钟满足 $t=k\Delta t$，并在指定步数后结束。 |
| `FixedStepClockTest.RejectsNonPositiveTimeStep` | 固定时钟拒绝零或负步长。 |
| `FixedStepClockTest.RejectsAdvanceAfterCompletion` | 已结束的固定时钟不能继续推进。 |
| `FixedStepPlasmaSimulationTest.AdvanceOneStepEvaluatesReactionAndUpdatesSpecies` | 单步依次完成电场、反应源和固定步长输运更新。 |
| `FixedStepPlasmaSimulationTest.ReactionRateIsReevaluatedFromUpdatedStateEveryStep` | 每一步反应率都读取最新密度，而非复用旧值。 |
| `FixedStepPlasmaSimulationTest.ReactionEvaluatorSeesCurrentElectricField` | 反应模型读取的是本步泊松求解得到的电场。 |
| `FixedStepPlasmaSimulationTest.RunAdvancesUntilClockIsFinished` | `run()` 严格运行到固定时钟指定的总步数。 |
| `FixedStepPlasmaSimulationTest.EmitsOrderedTraceEventsAtEachPipelineStage` | 固定步长仿真按 run、步开始、电静力、反应率、源项和步提交的顺序发出结构化事件，并报告提交后的时间。 |
| `FixedStepPlasmaSimulationTest.EmitsReturnedDiagnosticBeforeStepFailureTrace` | 固定步长仿真消费求解结果中的底层 diagnostic，附加步数和残差上下文，在 `step.failed` 前输出且不提交失败步。 |
| `FixedStepPlasmaSimulationTest.SolvesPoissonExactlyOncePerTimeStepAndReusesFactorization` | 每步只解一次泊松方程，且固定矩阵复用分析和分解。 |
| `FixedStepPlasmaSimulationTest.RejectsClockTimeStepDifferentFromTransportTimeStep` | 仿真时钟与输运推进器的 $\Delta t$ 必须一致。 |
| `FixedStepPlasmaSimulationTest.RejectsAdvanceAfterSimulationFinished` | 固定步长仿真结束后不能再次推进。 |
| `FixedStepPlasmaSimulationTest.StatisticsRejectFieldsWithoutPhysicalMetadata` | 启用物理统计时，simulation 必须在运行前拒绝缺少预期 `QuantityKind + precise_unit` metadata 的场，避免生成无单位或错单位统计。 |
| `AdaptiveTimeClockTest.LandsExactlyOnEndTime` | 可变步长累加后时钟精确吸附到终止时刻。 |
| `AdaptiveStepPlasmaSimulationTest.RunSelectsVariableStepsAndReachesEndTime` | 自适应仿真选取可变步长、使用末步截断并准确到达终止时间。 |
| `AdaptiveStepPlasmaSimulationTest.EmitsOrderedTraceEventsWithTimeStepDiagnostics` | 自适应仿真逐阶段发出 trace，且 `timestep.selected` 和 `step.completed` 分别准确报告实际步长与提交后的时间。 |
| `AdaptiveStepPlasmaSimulationTest.EmitsReturnedDiagnosticBeforeStepFailureTrace` | 自适应仿真同样先输出返回的根因 diagnostic，再输出阶段失败 trace，并保持时钟未推进。 |
| `AdaptiveStepPlasmaSimulationTest.AdvanceOneStepEvaluatesReactionAndUpdatesSpecies` | 单个自适应步完成电静力、反应率、化学源项和物种更新，且记录实际步长。 |
| `AdaptiveStepPlasmaSimulationTest.ReactionRateIsReevaluatedFromUpdatedStateEveryStep` | 每一个自适应步都从最新粒子密度重新计算反应率。 |
| `AdaptiveStepPlasmaSimulationTest.ReactionEvaluatorSeesCurrentElectricField` | 自适应反应模型读取的电场来自同一时间层的泊松求解。 |
| `AdaptiveStepPlasmaSimulationTest.ReactionSinkLimitsTimeStepAndPreservesNonNegativeDensity` | 强反应损失会收紧正性步长限制，更新后所有物种密度保持非负。 |
| `AdaptiveStepPlasmaSimulationTest.SolvesPoissonOncePerAdaptiveStepAndReusesFactorization` | 自适应仿真每步只进行一次泊松求解，并复用固定矩阵的符号分析与数值分解。 |
| `AdaptiveStepPlasmaSimulationTest.RejectsAdvanceAfterSimulationFinished` | 自适应仿真到达终止时刻后拒绝继续推进。 |
| `AdaptiveStepPlasmaSimulation64x64Test.SolvesUniformElectronImpactIonizationAcrossMultipleSteps` | 在 `poisson_64x64.msh` 的 4096 单元上运行 5 个均匀成对电离步，检查每步步长、终止时间、末步反应率/源项、逐单元/全域密度、电中性及零电势电场。 |
| `AdaptiveStepPlasmaSimulation64x64Test.ParallelPlate400VDrivesOppositeDriftAndIonizationInCentimeterMesh` | 将 `poisson_64x64.msh` 的坐标解释为厘米，在左右端施加 $0/400\,\mathrm{V}$、上下绝缘电势边界；验证 mp-units 配置量向裸代数参数的换算、密度/电势/电场/反应率/源项的运行期 metadata 传播，以及受输运稳定性约束的多步推进、非负密度和电子相对离子向右漂移。测试启用 $1\,\mathrm{cm}$ 面外厚度统计，并验证普通 trace/diagnostic 只写入当前工作目录的 `test.diag.log`；`test.statistics.log` 按 step 输出 species、charge、potential 和 electric-field 七列统计表及单位，不含普通事件、`DETAILS` 或零值异常计数。 |
| `AdaptiveStepPlasmaSimulation64x64Test.ParallelPlate400VPreservesLowerSolverDiagnosticAfterPhysicalStep` | 使用与上述 \(0/400\,\mathrm{V}\) 厘米网格、电离反应和初始物种完全相同的配置；先以真实 CHOLMOD 完成一个物理时间步，再注入一次底层求解失败，验证 `linalg` 域 diagnostic、可读的 `SolveFailed` 状态与根因消息先于 simulation 失败事件输出。结构化失败 trace 写入当前工作目录的 `test.failure.diag.log`。 |
