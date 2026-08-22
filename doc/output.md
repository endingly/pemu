# VTKHDF 场输出

`submod/output` 是独立于 `trace` 的持久化模块。第一阶段只实现面向 ParaView 的
VTKHDF field dump；`trace` 不保存网格或场数据，只在写入成功后同步接收一个
`output.completed` 事件及 `path`、`step`、`time` 三个属性。

## 数据映射

`toVtkUnstructuredGrid` 使用 `IMesh::vertex` 和 `IMesh::cellVertices` 构造 VTK
points 与 polygon cells。当前 MOAB 后端仍只接受平面二维 polygon mesh，但输出端
不依赖 MOAB handle。

| pemu 数据 | VTKHDF 位置 | 语义 |
| --- | --- | --- |
| `CellField<double>` | `CellData` | 数组名取自 `FieldMetadata::name`。 |
| `FaceField<double>` 原始值 | `FieldData` | 保持稠密 `FaceId` 顺序，不伪装成 `CellData`。 |
| face owner/neighbor/center/area/boundary | `FieldData` | 使原始 face 数组保持可解释、可重建的数据关联。 |
| 可选 face 可视化副本 | `CellData` | 对每个 cell 的相邻 face 做面积加权平均，并使用独立数组名。 |

每个场同时在四个并行的一维字符串数组中登记 association、原始 name、
`QuantityKind` 与 LLNL `precise_unit` 字符串：

- `pemu_field_association`
- `pemu_field_name`
- `pemu_field_quantity_kind`
- `pemu_field_unit`

官方 writer 会按 VTKHDF 规范替换数据数组名中的 `.` 与 `/`；metadata 的 name 列
始终保存未修改的原始 field name，因此语义名称不会丢失。

`pemu_step` 与 `pemu_time` 位于每个时间步的 `FieldData`，time 还写入 VTK 的标准
`DATA_TIME_STEP` information。`IFieldOutputWriter::openSeries` 建立一次输出会话，按
step/time 递增地 `append` 快照，最后 `finish` 成一个 `.vtkhdf`/`.hdf` 时间序列文件。
单快照 `write` 是上述 series 生命周期的便利封装。

## Simulation 调度

`FixedStepPlasmaSimulation` 与 `AdaptiveStepPlasmaSimulation` 通过构造参数
`simulation::FieldOutputOptions` 启用输出。配置持有一个非拥有的
`IFieldOutputWriter`，writer 必须比 simulation 活得更久；默认的空 writer 指针保持
输出关闭，因此不会改变既有数值推进或产生文件。

启用后，simulation 会在电势、电场与当前密度同步的时刻采集：可选初始状态、每
`every_steps` 个已完成状态一次，以及可选终态。终态会额外进行一次电静力更新，确保
终态密度与输出的电势/电场属于同一时刻。全部被选中的状态按时间顺序进入同一个文件，
路径固定为 `<directory>/<file_stem>.vtkhdf`；默认目录为 `output`，默认 stem 为
`plasma`。若同一路径已存在，仍由 writer 的 `overwrite` 选项决定是否允许覆盖。

```cpp
pemu::output::VtkHdfWriter writer;
pemu::simulation::FieldOutputOptions output_options{
    .writer = &writer,
    .directory = "results",
    .file_stem = "discharge",
    .every_steps = 10,
    .write_initial = true,
    .write_final = true,
};

// output_options 是 Fixed/AdaptiveStepPlasmaSimulation 构造函数的最后一个参数。
// writer 必须在 simulation 完成前保持存活。
```

快照包含每个物种的数密度、空间电荷密度、电势、法向电场和各物种法向漂移速度。
物种字段使用带 SpeciesId 的稳定输出名称，例如 `species_0_e_number_density`，以避免
同类物种字段共用 metadata 名称时产生 VTK 数组重名。整个时间序列成功完成后，若配置了
trace sink，会发出一次轻量的 `output.completed` 事件，属性中的 step/time 对应末个快照。

## 转换与内存模型

每个 series 只调用一次 `toVtkUnstructuredGrid`：MOAB 网格的 points/connectivity 由会话
持有并在所有时间步复用。结束写出时，temporal source 对这份 topology 做 shallow copy，
动态数值数组则通过 `vtkDoubleArray::SetArray(..., save=1)` 直接包装快照缓冲区；因此不会在
每个时间步重新构造完整 `vtkUnstructuredGrid`，也不会在已有快照缓冲区与 VTK 数组之间再
复制一次。共享 topology 的 mesh MTime 保持不变，`vtkHDFWriter` 因而只把静态 VTK 网格
几何和连通关系写入一次。

当前实现为了使用 VTK 9.4 的 `WriteAllTimeSteps` 管线，会先把被选择的动态场快照保存在内存，
到 simulation 终态再完成文件。这把每步成本限制为动态场复制与可选 face→cell 可视化计算，
但内存随快照数增长。后续大规模优化应在保持现有 series API 的前提下改为在线 HDF5 append，
并进一步研究静态 face 辅助数组布局、chunk/compression 及异步 I/O；无需再次改变 simulation
调度接口。

## FaceField 约束

VTK unstructured grid 没有 face-data association，因此原始 `FaceField` 只能作为
FieldData 保存。若 `FaceFieldSelection::include_cell_centered_visualization` 为真，
writer 额外产生一个 cell-centered `CellData` 数组；这个数组明确标记为
`cell_visualization`，不能代替原始 face 数据参与 restart 或守恒计算。

`CellFieldSelection::name` 与 `FaceFieldSelection::name` 可以仅覆盖写出名称，不复制
数值数组，也不改变字段本身的物理量 metadata。这让 simulation 能为同质的物种字段
提供唯一文件名，同时保留其单位和 quantity kind。

## Checkpoint 边界

`IFieldOutputWriter` 与模板接口 `ICheckpointWriter<State>` 相互独立。当前没有任何
checkpoint writer 或 restart reader 实现，也不把 VTKHDF 可视化文件声明为
checkpoint。后续确定完整 simulation state 与一致性协议后，可在不改变 field dump
接口的前提下实现 checkpoint/restart。

## 可读性验证

`ParaViewReadabilityTest` 写出真实 VTKHDF 文件，再用官方 `vtkHDFReader` 完整读回并
验证网格、全部数据关联和 metadata；时间序列测试还逐 step 验证不同场值、`pemu_step`
及 time。ParaView 的 VTKHDF reader 使用同一 VTK I/O 实现，因此这些 round-trip 同时
作为无 GUI 的 ParaView 可读性回归测试。
