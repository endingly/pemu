# Output：dump 与 checkpoint

`submod/output` 是独立于 `trace` 的持久化模块，公共接口按用途分层：

- `pemu/output/common`：dump/checkpoint 共用的 stamp、record 与 mesh→VTK adapter；
- `pemu/output/dump` 与 `pemu::output::dump`：面向 ParaView 的时序场输出；
- `pemu/output/checkpoint` 与 `pemu::output::checkpoint`：可校验、可恢复的状态快照。

旧的顶层 output 头文件已经移除，所有调用方均使用上述分层入口。`trace` 不保存网格或
场数据；dump 写入成功后仍只同步接收 `output.completed` 及 `path`、`step`、`time`
三个属性。

## 数据映射

`toVtkUnstructuredGrid` 使用 `IMesh::vertex` 和 `IMesh::cellVertices` 构造 VTK
points 与 polygon cells。当前 MOAB 后端仍只接受平面二维 polygon mesh，但输出端
不依赖 MOAB handle。

| pemu 数据 | VTKHDF 位置 | 语义 |
| --- | --- | --- |
| `CellField<double>` | `CellData` | 默认数组名取自 `FieldMetadata::name`，selection 可覆盖。 |
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
保存未修改的输出语义名称（selection 提供 `meta_data` 时为 `meta_name`），因此语义名称
不会丢失。

`pemu_step` 与 `pemu_time` 位于每个时间步的 `FieldData`，time 还写入 VTK 的标准
`DATA_TIME_STEP` information。`output::dump::IWriter::openSeries` 建立一次输出会话，按
step/time 递增地 `append` 快照，最后 `finish` 成一个 `.vtkhdf`/`.hdf` 时间序列文件。
单快照 `write` 是上述 series 生命周期的便利封装。

## Simulation 调度

`FixedStepPlasmaSimulation` 与 `AdaptiveStepPlasmaSimulation` 通过构造参数
`simulation::FieldOutputOptions` 启用输出。配置持有一个非拥有的
`output::dump::IWriter`，writer 必须比 simulation 活得更久；默认的空 writer 指针保持
输出关闭，因此不会改变既有数值推进或产生文件。

启用后，simulation 会在电势、电场与当前密度同步的时刻采集：可选初始状态、每
`every_steps` 个已完成状态一次，以及可选终态。终态会额外进行一次电静力更新，确保
终态密度与输出的电势/电场属于同一时刻。全部被选中的状态按时间顺序进入同一个文件，
路径固定为 `<directory>/<file_stem>.vtkhdf`；默认目录为 `output`，默认 stem 为
`plasma`。若同一路径已存在，仍由 writer 的 `overwrite` 选项决定是否允许覆盖。

```cpp
pemu::output::dump::VtkHdfWriter writer;
pemu::simulation::FieldOutputOptions output_options{
    .writer = &writer,
    .directory = "results",
    .file_stem = "discharge",
    .every_steps = 10,
    .write_initial = true,
    .write_final = true,
};

// output_options 位于末尾的 checkpoint_options 之前。
// writer 必须在 simulation 完成前保持存活。
```

快照包含每个物种的数密度、空间电荷密度、电势、法向电场和各物种法向漂移速度。
物种字段使用带 SpeciesId 的稳定输出名称，例如 `species_0_e_number_density`，以避免
同类物种字段共用 metadata 名称时产生 VTK 数组重名。整个时间序列成功完成后，若配置了
trace sink，会发出一次轻量的 `output.completed` 事件，属性中的 step/time 对应末个快照。

simulation 的运行状态、`start/advance/pause/stop/run` 语义以及自动 checkpoint 调度见
[Simulation workflow](simulation-workflow.md)。

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

如需为某个输出数组声明不同于源字段的语义，可在 selection 的 `meta_data` 中提供非空
`meta_name`、`quantity_kind` 和 `unit`。writer 会将它们写入 `pemu_field_*` metadata
数组，而不会修改源字段或该数组的 VTK 名称；空 `meta_name` 保持源字段 metadata 的既有
行为。

## VTKHDF checkpoint

checkpoint 使用 VTKHDF，而不是直接调用 HDF5 C API。纯 HDF5 只定义容器，仍需项目自行
发明 group/dataset 协议；VTKHDF 已提供稳定的网格与数据关联结构，并且当前依赖和
round-trip 路径已经覆盖它。checkpoint 与 dump 虽共享容器和 mesh adapter，但 API、
manifest 和语义互不替代：普通 dump 文件没有 checkpoint marker，reader 会拒绝恢复。

`output::checkpoint::VtkHdfWriter` 将以下内容写入单快照 `.vtkhdf`/`.hdf`：

- `pemu.vtkhdf.checkpoint` marker、格式版本、step/time；
- 完整 points/cells，以及按 FaceId 排列的 owner、neighbor、center、area、normal 和
  boundary id 校验数组；
- 所选 CellField/FaceField 的原始 double 值；
- 可选 scalar workflow 状态（例如自适应推进的 `previous_dt`）；
- association、稳定 key、原字段 metadata name、quantity kind、unit 与值数量。

覆盖已有 checkpoint 时，writer 先在目标目录写入并校验临时 VTKHDF，再以文件系统原子
rename 提交；写入失败会清理临时文件并保留上一份有效 checkpoint。

VTKHDF 会改写数组名中的 `.` 和 `/`，所以 checkpoint key 明确禁止这两个字符，避免恢复
时出现名称歧义。`VtkHdfReader::inspect` 可只读取并验证 manifest；`restore` 在写入任何
目标字段之前统一校验格式版本、完整 mesh/FaceId 顺序、key、尺寸和物理 metadata。默认
要求调用方提供所有已保存字段，因此恢复要么全部成功，要么保持目标状态不变。

```cpp
using namespace pemu::output::checkpoint;

VtkHdfWriter writer;
writer.write({
    .mesh = &mesh,
    .path = "state.vtkhdf",
    .stamp = {.step = step, .time = time},
    .cell_fields = cell_sources,
    .face_fields = face_sources,
});

VtkHdfReader reader;
const auto manifest = reader.inspect("state.vtkhdf");
reader.restore("state.vtkhdf", {
    .mesh = &mesh,
    .cell_fields = cell_targets,
    .face_fields = face_targets,
});
```

checkpoint 保存由调用方明确选择的数值状态；外部配置、reaction network、边界条件和求解器
实例仍由应用构造。底层 reader/writer 仍可独立使用；simulation workflow 会额外写入 workflow
schema，并将 species ID 与名称编码进密度 key，从而拒绝物种缺失或顺序不一致的恢复；它还会
恢复 clock，并为自适应推进保存/恢复步长历史。

## 可读性验证

`ParaViewReadabilityTest` 写出真实 VTKHDF 文件，再用官方 `vtkHDFReader` 完整读回并
验证网格、全部数据关联和 metadata；时间序列测试还逐 step 验证不同场值、`pemu_step`
及 time。ParaView 的 VTKHDF reader 使用同一 VTK I/O 实现，因此这些 round-trip 同时
作为无 GUI 的 ParaView 可读性回归测试。
