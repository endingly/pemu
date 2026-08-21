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

`pemu_step` 与 `pemu_time` 位于 `FieldData`，time 还写入 VTK 的标准
`DATA_TIME_STEP` information。一次 `write` 生成一个 `.vtkhdf`/`.hdf` 文件，调用层
可按 step 或 time 选择路径并调度输出。

## FaceField 约束

VTK unstructured grid 没有 face-data association，因此原始 `FaceField` 只能作为
FieldData 保存。若 `FaceFieldSelection::include_cell_centered_visualization` 为真，
writer 额外产生一个 cell-centered `CellData` 数组；这个数组明确标记为
`cell_visualization`，不能代替原始 face 数据参与 restart 或守恒计算。

## Checkpoint 边界

`IFieldOutputWriter` 与模板接口 `ICheckpointWriter<State>` 相互独立。当前没有任何
checkpoint writer 或 restart reader 实现，也不把 VTKHDF 可视化文件声明为
checkpoint。后续确定完整 simulation state 与一致性协议后，可在不改变 field dump
接口的前提下实现 checkpoint/restart。

## 可读性验证

`ParaViewReadabilityTest` 写出真实 VTKHDF 文件，再用官方 `vtkHDFReader` 完整读回并
验证网格、全部数据关联和 metadata。ParaView 的 VTKHDF reader 使用同一 VTK I/O
实现，因此该 round-trip 同时作为无 GUI 的 ParaView 可读性回归测试。
