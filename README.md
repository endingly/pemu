# pemu 

`pemu` 是一个面向二维静电低温等离子体的 C++23 CPU 仿真内核。它以有限体积法
（finite-volume method，FVM）为核心：在任意多边形控制体上守恒地离散通量，
再将所得稀疏线性系统交给 SuiteSparse 求解。

当前实现覆盖：

- 稳态泊松/扩散方程 $-∇·(ε∇φ)=ρ$
- 标量对流通量的一阶迎风离散；
- 瞬态扩散方程的后向欧拉时间推进；
- 电子能量密度的显式 Scharfetter–Gummel 输运、稳定性与非负性约束；
- 多物种反应网络、E/N 或 Te 查表反应率和通用质量作用律；
- 壁面粒子/能量损失与离子诱导二次电子发射；
- 固定/自适应 Simulation workflow、结构化诊断、VTKHDF dump 与 checkpoint 恢复；
- Gmsh 网格的 MOAB 读取、边界物理组及几何/拓扑预处理；
- 基于官方 `vtkHDFWriter` 的独立 VTKHDF 场输出与 ParaView 可读性验证。

## 阅读顺序

1. [数学模型与离散化](doc/mathematics-and-discretization.md)：方程、符号、面通量和离散精度；
2. [网格、边界与程序组织](doc/mesh-boundary-and-architecture.md)：数据约定及模块边界；
3. [诊断、日志与中间状态追踪](doc/diagnostics-and-tracing.md)：结构化事件、sink 与仿真耦合；
4. [Output 持久化](doc/output.md)：VTKHDF dump、metadata、版本化 checkpoint 与恢复；
5. [Simulation workflow](doc/simulation-workflow.md)：轻量状态机、运行控制与 checkpoint 调度；
6. [电子能量方程](doc/electron-energy-equation.md)：M13 模型、SG 离散与显式稳定性；
7. [壁面模型](doc/plasma-wall-boundary.md)：粒子、能量损失与二次发射；
8. [Plasma benchmark](doc/plasma-benchmarks.md)：M16 解析验证与公开 benchmark 接入边界；
9. [测试用例契约](doc/testcases.md)：全部测试各自保证的性质。

## 当前适用范围

`MoabMesh` 后端读取平面 XY 二维网格；`cellVolume` 在二维中表示面积，`faceArea`
表示边长。`AxisymmetricMeshView` 可将 `x=r\ge0` 的二维网格转换为绕轴控制体度量。
扩散系数在现有 API 中仍为非负的空间常数。对内面的两点
法使用控制体中心连线在面法向上的投影距离，因而最自然地适用于正交或接近正交
的网格；非正交修正、各向异性扩散和三维几何尚未实现。

边界物理组约定为：左/右/下/上分别对应 `BoundaryId` 1/2/3/4。测试网格
`data/meshfiles/two_quads.msh` 是两个相邻的单位四边形，收敛性测试使用
`poisson_{8,16,32,64}x{8,16,32,64}.msh`。

## 构建与运行

项目使用 CMake 3.30、C++23、Ninja 与 vcpkg manifest 模式。依赖为 Eigen3、MOAB、
SuiteSparse（CHOLMOD、UMFPACK）、GoogleTest 和 fmt；版本基线记录于
[`vcpkg.json`](vcpkg.json)。配置预设要求环境变量 `VCPKG_ROOT` 指向 vcpkg 根目录。

```bash
cmake --preset ci-linux-gcc-debug
cmake --build --preset ci-linux-gcc-debug
ctest --preset ci-linux-gcc-debug --output-on-failure
```

也可以用一条命令执行同一 Debug 工作流：

```bash
cmake --workflow --preset ci-linux-gcc-debug
```

Release 预设将命令中的 `debug` 替换为 `release`。测试数据目录在 CMake 配置时以
`PEMU_MESH_TEST_DATA_DIR` 传入各测试目标，因此不依赖执行时的当前工作目录。

按组件筛选测试可使用：

```bash
ctest --test-dir out/build/ci-linux-gcc-debug -R 'Poisson|Transient' --output-on-failure
```

完整的测试意图、方程背景与阈值见 [测试用例契约](doc/testcases.md)

在使用 vcpkg 之前需要一些系统组件：

```bash
apt install libtool libtool-bin autoconf autoconf-archive automake libtool
```
