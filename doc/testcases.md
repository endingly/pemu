# 测试用例契约

本页逐一说明仓库中所有 296 个 GoogleTest 用例所守护的契约。数值阈值并非一般性
精度承诺，而是当前测试网格、双精度实现和制造解下的回归界限。`two_quads.msh`
包含两个相邻单位方形；大多数几何与算子测试以它为夹具。

## 日志与追踪（`submod/trace/test/test_trace.cpp`）

| 用例 | 保证 |
| --- | --- |
| `TraceSinkTest.AnyTraceSinkForwardsNoexceptSink` | 非模板 API 使用的 owning type erasure 同时转发事件与统一的 `flush()` 契约；副本共享底层缓冲状态，空 sink 的两个操作均无副作用。 |
| `TraceSinkTest.NullSinkAcceptsStructuredEvents` | 默认空 sink 满足 trace concept，能同步接受结构化事件且不产生副作用。 |
| `PhysicalVolumeSemanticsTest.Extrudes2DMeasuresAndKeeps3DMeasures` | 二维网格的单元面积/面长度乘以面外厚度得到物理体积/面积，而三维网格保持原生几何测度。 |
| `StatisticsOptionsTest.SamplesOnlyAtConfiguredStepInterval` | 统计探针只在启用且步号满足采样间隔时运行，并拒绝零间隔等非法配置。 |
| `ScalarFieldStatisticsTest.ComputesWeightedExtremaIntegralMeanAndRms` | 标量统计器在非均匀物理权重下正确计算极值、最大绝对值、积分、L1 积分、加权均值和 RMS。 |
| `ScalarFieldStatisticsTest.DerivesIntegralUnitOutsideSampleLoop` | 统计器只在收尾边界附加并组合运行期单位，验证 $\mathrm{C/cm^3}\times\mathrm{cm^3}=\mathrm C$，而逐样本 `add` 仍只接收两个 `double`。 |
| `ScalarFieldStatisticsTest.ReportsNonFiniteValuesAndInvalidWeights` | 标量统计器显式计数负值、非有限值和非法物理权重，不把异常样本静默混入积分。 |
| `TraceSinkTest.OstreamSinkFormatsTabularStructuredEvents` | 普通事件 renderer 保持原有统一事件表格式，将步数、时间、步长、终止时间、求解状态和残差放入稳定对齐列，其余属性进入 DETAILS。 |
| `TraceSinkTest.StatisticsRendererGroupsRowsAndFormatsDedicatedColumns` | 统计 renderer 按 step 缓冲分组，以 FIELD/UNIT/MIN/MAX/MEAN/RMS/INTEGRAL 独立表格展示多个场；固定采样与体积语义元数据只输出一次，单位来自事件属性而非字段名猜测，零异常计数和 DETAILS 列均不出现。 |
| `TraceSinkTest.StatisticsRendererOnlyShowsNonZeroCountsAsDiagnostics` | 统计质量计数仅在非零时进入该 step 的 `STATISTICS DIAGNOSTICS` 小节，并保留 Warning/Diagnostic 级别；值为零的 `non_finite` 不输出。 |
| `TraceSinkTest.OstreamSinkFlushesOnlyWhenContractRequestsIt` | ostream sink 只在统一契约的显式 `flush()` 调用时刷新流，不解析 simulation 专属事件名。 |
| `TraceSinkTest.SplitSinkRoutesChannelsAndForwardsFlush` | 双通道 sink 仅按显式 channel 路由事件，并把统一的 `flush()` 同时转发给两个子 sink，不承担 simulation 刷新策略。 |
| `EnumStringTest.DomainsAndSeveritiesHaveStableNames` | `DiagDomain` 与 `Severity` 均通过 `pemu::to_string` 映射为稳定、可读的字符串，可供过滤和机器处理使用。 |
| `TraceSinkTest.NullSinkAcceptsStructuredDiagnostics` | 同一个默认空 trace sink 能无副作用地接受标记为 diagnostic 的结构化事件。 |
| `TraceSinkTest.OstreamSinkFormatsDiagnosticEvent` | diagnostic renderer 保留根因消息与数值上下文，但与前后的普通事件共用同一个旧式事件表头；普通→diagnostic→普通切换不会插入 TRACE/DIAGNOSTIC 分段标题。 |

## 单位（`submod/unit/test/test_unit.cpp`）

| 用例 | 保证 |
| --- | --- |
| `PlasmaQuantitiesTest.DefineDimensionallyConsistentReferences` | 项目定义的数密度、反应率密度、法向电场/漂移速度及电子平均能量、能量密度和能量源 quantity specification 只接受量纲相容的 mp-units unit。 |
| `MpUnitsBridgeTest.ProducesRuntimeMetadataWithoutStringParsing` | mp-units reference 在编译期直接桥接为可平凡复制的 `QuantityKind + precise_unit`；包括 $\mathrm{eV/cm^3}$ 在内的运行期单位正确，且不依赖字符串解析。 |
| `RuntimeUnitStringTest.UsesStablePlasmaEnergySpellings` | 项目 formatter 为电子能量密度、源项及积分功率输出稳定的 `eV/cm^3`、`eV/(cm^3*s)` 与 `eV/s`，未登记单位仍回退到 LLNL formatter。 |

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
| `AxisymmetricMeshViewTest.PreservesMeridionalTopologyAndCoordinates` | 轴对称视图只改变物理度量，保持底层二维拓扑、边界元数据、中心与法向。 |
| `AxisymmetricMeshViewTest.RevolvesPlanarCellsIntoAnnularVolumes` | 两个径向单位单元旋转后的总体积等于半径 2、长度 1 圆柱的 $4\pi$。 |
| `AxisymmetricMeshViewTest.RevolvesFacesIntoPhysicalSurfaceAreas` | 轴线、内部、外壁和端面的旋转面积分别满足解析几何值。 |
| `AxisymmetricMeshViewTest.RejectsOutOfRangeMeasureQueries` | 轴对称预计算度量保持 `IMesh` 的越界检查契约。 |
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
| `FieldTest.MpUnitsBridgeCoversCanonicalPlasmaFieldMetadata` | 编译期 bridge 不经字符串解析即可把数密度、电荷密度、法向电场和电子能量密度 reference 映射为正确的语义枚举与 LLNL `precise_unit`。 |
| `FieldTest.QuantityBoundaryConvertsToFieldStorageUnit` | `fillQuantity`/`setQuantity` 可将米输入换算为按厘米存储的裸值，`quantityAt` 能重新包装和换算输出。 |
| `FieldTest.QuantityBoundaryRejectsReferenceDifferentFromMetadata` | 边界适配器会拒绝单位不同的 reference，也会拒绝同为伏特但 quantity specification 不同的 reference。 |

## VTKHDF 场输出（`submod/output/test/test_vtkhdf_writer.cpp`）

| 用例 | 保证 |
| --- | --- |
| `MeshAdapterTest.ConvertsIMeshToPolygonalUnstructuredGrid` | `IMesh` 的顶点、cell→vertex 连通关系被无损映射为 `vtkUnstructuredGrid` polygon cells。 |
| `ParaViewReadabilityTest.OfficialVtkHdfReaderRoundTripsFieldsMetadataAndTime` | 官方 `vtkHDFWriter` 产物能由 ParaView 使用的官方 `vtkHDFReader` 读回；CellData、原始 face FieldData/拓扑、面积加权 cell-centered 可视化副本、field name、运行期 quantity/unit 以及 step/time 均保持一致。 |
| `ParaViewReadabilityTest.SingleVtkHdfFileRoundTripsAnOrderedTemporalSeries` | 两个不同 step/time 和场值通过一个 series 写入唯一的 VTKHDF 文件，官方 reader 能枚举并逐步读回。 |
| `OutputTraceTest.EmitsOnlyLightweightCompletionContextAfterWrite` | 成功输出后只向 trace 发出 `output.completed` 及 `path/step/time` 三个轻量属性，不把场数组塞入 trace。 |
| `VtkHdfWriterTest.UsesCellFieldSelectionNameWithoutChangingMetadata` | 写出名称可覆盖而无需复制数据，且 quantity/unit metadata 仍与原字段一致。 |
| `VtkHdfWriterTest.WritesExplicitSelectionMetadata` | selection 提供的 semantic name、quantity kind 与 unit 会写入 VTKHDF metadata，且 VTK 数组名保持独立。 |
| `VtkHdfWriterTest.RefusesOverwriteUnlessExplicitlyEnabled` | writer 默认拒绝覆盖已有文件，只有请求显式设置 `overwrite` 时才覆盖。 |

## VTKHDF checkpoint（`submod/output/test/test_checkpoint_vtkhdf.cpp`）

| 用例 | 保证 |
| --- | --- |
| `VtkHdfCheckpointTest.RoundTripsManifestCellAndFaceState` | 版本、step/time、mesh 尺寸、字段 metadata，以及 cell/face/scalar 原始状态能通过独立 checkpoint reader 完整恢复。 |
| `VtkHdfCheckpointTest.AtomicallyReplacesExistingCheckpoint` | 滚动 checkpoint 通过临时文件原子替换已有状态，恢复得到新值且目标目录不残留临时文件。 |
| `VtkHdfCheckpointTest.ValidatesAllTargetsBeforeChangingAnyField` | restore 在复制前校验全部目标；任一 metadata 不匹配时不产生部分恢复。 |
| `VtkHdfCheckpointTest.RejectsKeysThatVtkHdfWouldRewrite` | checkpoint 拒绝含 `.` 或 `/` 的 key，避免 VTKHDF 名称规范化造成恢复歧义。 |
| `VtkHdfCheckpointTest.RejectsVisualizationDumpWithoutCheckpointManifest` | checkpoint reader 拒绝没有版本化 checkpoint manifest 的普通 VTKHDF dump，确保两类文件不能误用。 |

## 离散算子（`submod/discretization/test/test_operator.cpp`）

| 用例 | 保证 |
| --- | --- |
| `DivergenceConservesInternalFlux` | 单个内面通量在全域积分散度中精确相消。 |
| `DivergenceUsesOwnerNeighborOrientation` | 正法向通量对 owner 为正、对 neighbor 为负。 |
| `DivergenceSatisfiesDiscreteGaussTheorem` | 体积积分散度等于边界积分通量（此夹具为 6）。 |
| `AxisymmetricDivergenceTest.RecoversRadialLinearFluxDivergence` | 轴对称度量下 $\boldsymbol\Gamma=r\mathbf e_r$ 的有限体积散度逐单元精确恢复为 2，并覆盖轴线零面积面。 |
| `AxisymmetricDivergenceTest.ReconstructsVectorMagnitudeAcrossZeroAreaAxisFace` | 场强重构忽略轴线零面积面的零物理权重，并由其余独立法向精确恢复常向量模。 |
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
| `OperatorTest.ScharfetterGummelNeumannAddsPrescribedDiffusiveFlux` | Neumann 值作为外向扩散通量与内部状态产生的漂移通量相加。 |
| `OperatorTest.ScharfetterGummelHomogeneousNeumannKeepsDriftFlux` | 齐次 Neumann 不会错误消除边界上的物理漂移通量。 |
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
| `ElectronFieldPowerTest.ReconstructsPositiveHeatingForElectronFluxAgainstField` | 用两单元正交网格验证完整面粒子通量与电场反向时，离散 $-\boldsymbol\Gamma_e\cdot\mathbf E$ 在两个单元中均得到正确的正加热。 |
| `ElectronFieldPowerTest.RejectsFieldsFromDifferentMeshes` | 电场功重构拒绝来自不同网格的通量、电场和目标场。 |
| `ElectronFieldPowerTest.RejectsNonFiniteFaceData` | 电场功重构在数值进入单元累加前拒绝非有限的面通量或电场。 |
| `ReactionTest.ElectronImpactIonizationComputesExpectedRate` | 电子碰撞电离率按给定电子密度、中性粒子密度和速率系数计算。 |
| `ReactionTest.ZeroElectronDensityProducesNoIonization` | 没有电子时电子碰撞电离率严格为零。 |
| `ReactionNetworkTest.PairIonizationConservesCharge` | 反应网络能识别成对电离的电荷守恒。 |
| `ReactionNetworkTest.StoresThirdBodyOrderIndependentlyFromNetStoichiometry` | 第三体可以参与速率定律而不产生净物种源。 |
| `ReactionNetworkTest.RejectsInvalidAndDuplicateKineticOrders` | 非法或重复的动力学级数在网络构造时被拒绝。 |
| `ReactionNetworkTest.DetectsChargeViolatingReaction` | 违反电荷守恒的化学计量关系会被检测。 |
| `ReactionNetworkTest.AccumulatesStoichiometricSources` | 反应率依据化学计量数正确转换为各物种源项。 |
| `ReactionNetworkTest.MultipleReactionsAccumulateCorrectly` | 多个反应对同一物种的贡献正确求和。 |
| `ReactionKineticsTest.EvaluatesArbitraryThreeBodyMassActionRate` | 通用质量作用律正确组装三个空间密度因子。 |
| `ReactionKineticsTest.AssemblerEvaluatesThreeBodyLawAndLeavesThirdBodySourceAbsent` | 三体 assembler 计算反应率，同时保持第三体净源为零。 |
| `ReactionKineticsTest.AssemblerRejectsBinaryCoefficientUnitForThreeBodyLaw` | 三体反应拒绝错误的 `cm^3/s` 二体系数单位。 |
| `WallFluxAssemblerTest.AssemblesPrimaryLossAndSecondaryParticleAndEnergyInflux` | 壁面组装器根据主损失通量同时生成二次粒子入流和发射能量入流。 |
| `WallFluxAssemblerTest.RequiresEvaluationBeforeCompleteFluxQueries` | 首次组装前或失败刷新后拒绝完整通量查询，避免混用当前主损失和过期二次发射。 |
| `WallFluxAssemblerTest.InvalidDensityLeavesPublishedInfluxUnchanged` | 非法密度导致壁面组装失败时，上一份已发布通量保持不变。 |
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
| `ExplicitElectronEnergyStepperTest.ConvertsBetweenEnergyDensityAndMeanEnergy` | $w_e=n_e\bar\varepsilon_e$ 双向转换正确，真空/低密度单元的平均能量按约定置零。 |
| `ExplicitElectronEnergyStepperTest.MaxwellianClosureScalesCompleteTransportOperator` | 默认 $5/3$ Maxwellian closure 同时缩放能量漂移速度和扩散系数，保持面 Péclet 数并按比例缩放完整 SG 算子。 |
| `ExplicitElectronEnergyStepperTest.RejectsInvalidEnergyBoundaryAtConstruction` | 能量输运构造时拒绝负 Dirichlet、非有限或缺失的边界条件。 |
| `ExplicitElectronEnergyStepperTest.HomogeneousNeumannPreservesConstantZeroDriftEnergy` | 零漂移下齐次 Neumann 保持常电子能量密度及全域积分。 |
| `ExplicitElectronEnergyStepperTest.PrescribedNeumannOutflowLimitsPositiveTimeStep` | 正向外流能量 Neumann 通量作为状态无关耗散进入电子能量 positivity 上限。 |
| `ExplicitElectronEnergyStepperTest.ConstantEnergyDensityRemainsConstant` | 零漂移、相容边界和零源下，常电子能量密度保持不变。 |
| `ExplicitElectronEnergyStepperTest.SourceChangesIntegratedEnergyByIntegratedSource` | 常值状态下总能量变化等于 $\Delta t$ 乘能量源的控制体积分。 |
| `ExplicitElectronEnergyStepperTest.DiffusionRedistributesSymmetricEnergyConservatively` | 对称边界下能量由高值单元扩散到低值单元，同时保持全域积分。 |
| `ExplicitElectronEnergyStepperTest.RejectsTransportCflViolationWithoutChangingEnergy` | 超过 SG 显式输运 CFL 时拒绝推进，能量场保持事务性不变。 |
| `ExplicitElectronEnergyStepperTest.RejectsNegativeSourceUpdateWithoutChangingEnergy` | 保守非负步长上限识别过强能量 sink，失败推进不留下部分更新。 |
| `FixedStepMultiSpeciesDriftDiffusionStepperTest.UniformNeutralPlasmaRemainsStationary` | 固定多物种推进器保持均匀中性状态。 |
| `FixedStepMultiSpeciesDriftDiffusionStepperTest.ChargeDensitySupportsMoreThanTwoSpecies` | 电荷密度对任意数量物种按 $\rho=\sum_s q_sn_s$ 求和。 |
| `FixedStepMultiSpeciesDriftDiffusionStepperTest.ImmobileSpeciesIsNotUpdated` | 标为 immobile 的物种不参与输运更新。 |
| `FixedStepMultiSpeciesDriftDiffusionStepperTest.AppliedPotentialProducesCorrectDriftForAllChargedSpecies` | 每个带电物种按迁移率和极性得到正确漂移速度。 |
| `FixedStepMultiSpeciesDriftDiffusionStepperTest.ComputesDriftForThreeTransportedChargedSpecies` | 漂移计算不隐含“只有电子和一种离子”的限制。 |
| `FixedStepMultiSpeciesDriftDiffusionStepperTest.CflFailureDoesNotPartiallyUpdateSpecies` | 多物种 CFL 预检具有原子性，不留下半更新状态。 |
| `FixedStepMultiSpeciesDriftDiffusionStepperTest.SourcePositivityFailureDoesNotPartiallyUpdateSpecies` | 固定步长先计算全部物种增量；任一强负源会产生负候选密度时拒绝整步且不修改其他物种。 |
| `FixedStepMultiSpeciesDriftDiffusionStepperTest.RejectsWrongSpeciesFieldCount` | 场数组物种数必须与 `SpeciesSet` 一致。 |
| `FixedStepMultiSpeciesDriftDiffusionStepperTest.RejectsFieldsFromDifferentMesh` | 多物种场不能跨网格传给推进器。 |
| `FixedStepMultiSpeciesDriftDiffusionStepperTest.MatchesLegacyTwoSpeciesSolver` | 更一般的固定多物种实现与固定双物种实现给出相同结果。 |
| `FixedStepMultiSpeciesDriftDiffusionStepperTest.PureNeumannGaugeIsAvailableFromCoupledStepper` | 纯 Neumann 规范配置可由多物种漂移扩散推进器传入内部泊松求解器，并得到零电势、零电场的中性均匀解。 |
| `FixedStepExplicitSpeciesContinuityStepperTest.NeumannOutflowContributesToFixedAndAdaptiveLossRates` | 普通 Neumann 边界保留漂移出流；固定、自适应估计分别计入对角损失与状态无关的规定通量 sink。 |
| `AdaptiveTimeStepControllerTest.AppliesSafetyAndGrowthLimits` | 时间步提议同时服从安全系数、稳定性限制和相对上一时刻的增长上限。 |
| `AdaptiveTimeStepControllerTest.TruncatesFinalStepToRemainingTime` | 时间步提议会被剩余时间截断，保证最终时刻精确对齐。 |
| `AdaptiveTimeStepControllerTest.RejectsNoPositiveAdmissibleTimeStep` | 稳定性或正性限制为零时拒绝不存在正步长的推进。 |
| `AdaptiveStepMultiSpeciesDriftDiffusionStepperTest.UsesConfiguredMaximumTimeStep` | 稳定性允许更大步长时，自适应推进器仍服从用户给定的 `max_dt`。 |
| `AdaptiveStepMultiSpeciesDriftDiffusionStepperTest.FinalTimeStepEqualsRemainingTime` | 最后一步被截断为剩余时间，避免越过终止时刻。 |
| `AdaptiveStepMultiSpeciesDriftDiffusionStepperTest.PrescribedNeumannOutflowLimitsSelectedTimeStep` | 多物种自适应 controller 将正向规定 Neumann 通量纳入 positivity limit，并用安全系数选出可提交步长。 |
| `LinearBoundaryFluxTest.AggregatesSeveralWallFacesOnOneOwnerCell` | 同一控制体上的多个壁面损失按 $\sum_f v_fA_f/V$ 聚合。 |
| `LinearBoundaryFluxTest.RejectsAggregateLossRateOverflow` | 单个壁面贡献均有限但控制体总损失率溢出时明确拒绝。 |
| `LinearWallContinuityTest.WallFluxChangesIntegratedParticlesExactly` | 线性壁面通量引起的全域粒子变化严格等于边界面积积分。 |
| `AdaptiveMultiSpeciesWallTransportTest.CouplesIonLossToSecondaryElectronFluxAndConservesFaceBalance` | 自适应多物种算子将离子壁面损失耦合为二次电子入流，并保持逐面收支。 |
| `FixedMultiSpeciesWallTransportTest.UsesTheSameSecondaryEmissionFluxContract` | 固定步长路径复用同一二次发射通量契约。 |
| `ElectronEnergyWallTest.AppliesIndependentEnergyLossAndSecondaryEmissionFlux` | 电子能量壁面使用独立损失速度并计入二次电子携带的能量。 |
| `ElectronEnergyWallTest.WallLossContributesToPositivityLimit` | 壁面能量耗散进入显式 positivity timestep。 |

## 固定步长与自适应步长仿真（`submod/simulation/test/test_simulation.cpp`）

| 用例 | 保证 |
| --- | --- |
| `FixedStepClockTest.AdvancesTimeFromStepIndex` | 固定时钟满足 $t=k\Delta t$，并在指定步数后结束。 |
| `FixedStepClockTest.RejectsNonPositiveTimeStep` | 固定时钟拒绝零或负步长。 |
| `FixedStepClockTest.RejectsAdvanceAfterCompletion` | 已结束的固定时钟不能继续推进。 |
| `FixedStepClockTest.RestoresInitialStepAndDerivedTime` | 固定时钟可从经校验的 checkpoint step 重建，并继续保持 $t=k\Delta t$。 |
| `FixedStepPlasmaSimulationTest.AdvanceEvaluatesReactionAndUpdatesSpecies` | `advance()` 单步依次完成电场、反应源和固定步长输运更新。 |
| `FixedStepPlasmaSimulationTest.ReactionRateIsReevaluatedFromUpdatedStateEveryStep` | 每一步反应率都读取最新密度，而非复用旧值。 |
| `FixedStepPlasmaSimulationTest.ReactionEvaluatorSeesCurrentElectricField` | 反应模型读取的是本步泊松求解得到的电场。 |
| `FixedStepPlasmaSimulationTest.RunAdvancesUntilClockIsFinished` | `run()` 严格运行到固定时钟指定的总步数。 |
| `FixedStepPlasmaSimulationTest.StateMachineSupportsStartPauseResumeAndStop` | workflow 在 `ready/running/paused/stopped` 间执行合法转换，暂停时拒绝推进，停止后不可再次运行；刷新策略在每两个完成步及 pause/stop 生命周期边界调用统一 sink 契约。 |
| `FixedStepPlasmaSimulationTest.StopFailureTransitionsWorkflowToFailed` | dump series 关闭失败时 `stop()` 不会伪装为成功停止，而是将 workflow 提交为 `failed`。 |
| `FixedStepPlasmaSimulationTest.EmitsOrderedTraceEventsAtEachPipelineStage` | 固定步长仿真按 run、步开始、电静力、反应率、源项和步提交的顺序发出结构化事件，并报告提交后的时间。 |
| `FixedStepPlasmaSimulationTest.WritesInitialPeriodicAndFinalFieldSnapshots` | 固定步 simulation 将同步状态在初始、周期及终态交给 output writer，并发出轻量完成 trace。 |
| `FixedStepPlasmaSimulationTest.EmitsReturnedDiagnosticBeforeStepFailureTrace` | 固定步长仿真消费求解结果中的底层 diagnostic，附加步数和残差上下文，在 `step.failed` 前输出且不提交失败步。 |
| `FixedStepPlasmaSimulationTest.SolvesPoissonExactlyOncePerTimeStepAndReusesFactorization` | 每步只解一次泊松方程，且固定矩阵复用分析和分解。 |
| `FixedStepPlasmaSimulationTest.RejectsClockTimeStepDifferentFromTransportTimeStep` | 仿真时钟与输运推进器的 $\Delta t$ 必须一致。 |
| `FixedStepPlasmaSimulationTest.RejectsAdvanceAfterSimulationFinished` | 固定步长仿真结束后不能再次推进。 |
| `FixedStepPlasmaSimulationTest.StatisticsRejectFieldsWithoutPhysicalMetadata` | 启用物理统计时，simulation 必须在运行前拒绝缺少预期 `QuantityKind + precise_unit` metadata 的场，避免生成无单位或错单位统计。 |
| `FixedStepPlasmaSimulationTest.CheckpointRejectsReorderedSpeciesWithoutChangingDensity` | checkpoint 将 species ID/name 纳入稳定 key；物种顺序不一致时恢复失败，且目标 density、clock 和状态保持不变。 |
| `FixedStepPlasmaSimulationTest.CheckpointRoundTripsCallerEnergyMetadata` | checkpoint 恢复临时场沿用调用方 energy metadata；字段名与 transport 内部展示名不同时仍可恢复自身写出的状态。 |
| `FixedStepPlasmaSimulationTest.AdvanceUpdatesMandatoryElectronEnergyAndMeanEnergy` | 固定步 workflow 用显式能量源推进必需的电子能量密度，并在物种与能量共同提交后刷新平均电子能量。 |
| `FixedStepPlasmaSimulationTest.SpeciesPositivityFailureLeavesCoupledStateUnchanged` | 固定步强反应 sink 失败时，species density、electron energy 和 clock 均保持在原时间层，workflow 进入 failed。 |
| `FixedStepPlasmaSimulationTest.RejectsElectronEnergyUnitDifferentFromTransportMetadata` | 即使 quantity kind 相同，调用方能量场使用与 transport 不同的存储单位也会在构造阶段被拒绝。 |
| `FixedStepPlasmaSimulationTest.RejectsMissingElectronEnergyAdditionalSourceEvaluator` | Simulation 构造阶段拒绝缺失电子能量附加源 evaluator 的配置，防止静默跳过碰撞/外部能量闭合。 |
| `AdaptiveTimeClockTest.LandsExactlyOnEndTime` | 可变步长累加后时钟精确吸附到终止时刻。 |
| `AdaptiveStepPlasmaSimulationTest.RunSelectsVariableStepsAndReachesEndTime` | 自适应仿真选取可变步长、使用末步截断并准确到达终止时间。 |
| `AdaptiveStepPlasmaSimulationTest.CheckpointWorkflowRestoresClockDensityAndTimeStepHistory` | simulation 自动保存自适应 checkpoint，并在新 workflow 中恢复 density、step/time 与 `previous_dt`，续算结果和未中断运行一致。 |
| `AdaptiveStepPlasmaSimulationTest.WritesScheduledSnapshotsAtAdaptiveStateTimes` | 自适应 simulation 按已完成状态的 step/time 调度初始、周期及终态输出。 |
| `AdaptiveStepPlasmaSimulationTest.EmitsOrderedTraceEventsWithTimeStepDiagnostics` | 自适应仿真逐阶段发出 trace，且 `timestep.selected` 和 `step.completed` 分别准确报告实际步长与提交后的时间。 |
| `AdaptiveStepPlasmaSimulationTest.EmitsReturnedDiagnosticBeforeStepFailureTrace` | 自适应仿真同样先输出返回的根因 diagnostic，再输出阶段失败 trace，并保持时钟未推进。 |
| `AdaptiveStepPlasmaSimulationTest.AdvanceEvaluatesReactionAndUpdatesSpecies` | 单个 `advance()` 完成电静力、反应率、化学源项和物种更新，且记录实际步长。 |
| `AdaptiveStepPlasmaSimulationTest.ReactionRateIsReevaluatedFromUpdatedStateEveryStep` | 每一个自适应步都从最新粒子密度重新计算反应率。 |
| `AdaptiveStepPlasmaSimulationTest.ReactionEvaluatorSeesCurrentElectricField` | 自适应反应模型读取的电场来自同一时间层的泊松求解。 |
| `AdaptiveStepPlasmaSimulationTest.ReactionSinkLimitsTimeStepAndPreservesNonNegativeDensity` | 强反应损失会收紧正性步长限制，更新后所有物种密度保持非负。 |
| `AdaptiveStepPlasmaSimulationTest.ElectronEnergyDepletionLimitsTheCoupledTimeStep` | 强电子能量耗散会进入统一自适应步长 proposal，并在物种本身允许更大步长时收紧耦合步长、保持能量非负。 |
| `AdaptiveStepPlasmaSimulationTest.SolvesPoissonOncePerAdaptiveStepAndReusesFactorization` | 自适应仿真每步只进行一次泊松求解，并复用固定矩阵的符号分析与数值分解。 |
| `AdaptiveStepPlasmaSimulationTest.RejectsAdvanceAfterSimulationFinished` | 自适应仿真到达终止时刻后拒绝继续推进。 |
| `FixedStepWallSimulationTest.AdvancesParticleAndElectronEnergyWallFluxInOneAtomicWorkflow` | Simulation 在同一步中耦合提交壁面粒子损失、二次电子和电子能量通量。 |
| `AdaptiveStepWallSimulationTest.AdvancesParticleAndElectronEnergyWallFluxInOneAtomicWorkflow` | 自适应 Simulation 复用相同壁面耦合，并把粒子与电子能量边界损失纳入统一步长和原子提交。 |
| `AdaptiveStepPlasmaSimulation64x64Test.SolvesUniformElectronImpactIonizationAcrossMultipleSteps` | 在 `poisson_64x64.msh` 的 4096 单元上运行 5 个均匀成对电离步，检查每步步长、终止时间、末步反应率/源项、逐单元/全域密度、电中性及零电势电场。 |
| `AdaptiveStepPlasmaSimulation64x64Test.ParallelPlate400VDrivesOppositeDriftAndIonizationInCentimeterMesh` | 将 `poisson_64x64.msh` 的坐标解释为厘米，在左右端施加 $0/400\,\mathrm{V}$、上下绝缘电势边界；验证 mp-units 配置量向裸代数参数的换算、密度/电势/电场/反应率/源项/电子能量的运行期 metadata 传播，以及受物种与能量联合稳定性约束的多步推进、非负密度和电子相对离子向右漂移。电子能量源由同一 SG 粒子通量形成的电场功与 $15.76\,\mathrm{eV}$ 每次电离的附加损失共同组成，并逐单元验证有限且为正。测试启用 $1\,\mathrm{cm}$ 面外厚度统计，并验证普通 trace/diagnostic 只写入当前工作目录的 `test.diag.log`；`test.statistics.log` 按 step 输出 species、charge、potential、electric-field、电子能量密度、平均能量和非零能量源统计，能量复合单位稳定显示为 `eV/cm^3` 与 `eV/(cm^3*s)`。真实 `VtkHdfWriter` 将包括电子能量密度和平均能量在内的 step 0 至终态同步状态写入 `<build>/submod/simulation/output/parallel-plate-400v/parallel-plate-400v.vtkhdf`，测试再用官方 reader 校验时间步数与终态 step。 |
| `PlasmaSimulation64x64CheckpointTest.ParallelPlate400VCheckpointRestoresAndContinuesSimulation` | 在关闭 dump/trace 的 $0/400\,\mathrm{V}$ 厘米网格固定步仿真中，由 simulation workflow 在第 3 步自动保存 checkpoint；新建 simulation 通过 `restoreCheckpoint()` 自动恢复物种密度、电子能量密度与 clock 后续跑至第 8 步，最终能量、density、电势、电荷和全部 face 电场/漂移速度与未中断运行逐项一致。 |
| `AdaptiveStepPlasmaSimulation64x64Test.ParallelPlate400VPreservesLowerSolverDiagnosticAfterPhysicalStep` | 使用与上述 \(0/400\,\mathrm{V}\) 厘米网格、电离反应和初始物种完全相同的配置；先以真实 CHOLMOD 完成一个物理时间步，再注入一次底层求解失败，验证 `linalg` 域 diagnostic、可读的 `SolveFailed` 状态与根因消息先于 simulation 失败事件输出。结构化失败 trace 写入当前工作目录的 `test.failure.diag.log`。 |

## Plasma benchmark（`submod/simulation/test/test_plasma_benchmarks.cpp`）

| 用例 | 保证 |
| --- | --- |
| `PlasmaBenchmarkTest.HomogeneousIonizationEnergyHasFirstOrderTemporalConvergence` | 完整 workflow 对均匀成对电离和等平均能量电子源同时匹配显式 Euler 离散真值，并对连续解析解表现出一阶时间收敛。 |
| `PlasmaBenchmarkTest.WallLossAndSecondaryEmissionHaveFirstOrderTemporalConvergence` | 粒子壁损失、二次电子和发射能量匹配三变量线性 ODE 的离散递推，并随时间步减半一阶收敛到连续解析解。 |
