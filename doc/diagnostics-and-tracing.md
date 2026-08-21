# 诊断、日志与中间状态追踪

## 1. 模块职责

`pemu::trace` 是一个独立的小型模块。它定义结构化 trace 事件、按模块分类的诊断事件、
严重级别、sink concept、通用标量统计器和少量通用 sink；不直接理解网格、场、物种、
反应或求解器。simulation 的轻量适配层负责把这些领域对象逐项送入统计器。

这个边界有意区分三类功能：

- **log**：供人阅读的运行消息，例如一步开始、一步结束或求解失败；
- **trace**：可由程序消费的连续运行状态，例如当前步数、时间、残差和自适应步长限制；
- **diagnostics**：需要关注的离散诊断，例如输入不合法、求解失败或数值状态异常；
- **result output**：电势、密度、电场等大规模网格场的快照与文件格式写出。

当前模块用同一种结构化事件和 sink 实现 trace 与 diagnostics，并提供将事件格式化为
日志行的最小 ostream sink；它不是通用日志系统。结果输出不应把整个场塞进事件属性；
后续更适合建立独立的 `pemu::output` 或 `pemu::io` 模块，由它负责快照频率、网格关联
和 VTK/HDF5 等格式。

## 2. 依赖方向

```text
┌───────────────┐
│  pemu::trace  │  定义事件、诊断域与 sink concept
└───────▲───────┘
        │ 公共事件类型
      linalg ─────► discretization ─────► equation ─────► simulation
 API 失败处形成诊断                              传播结果       消费并交给 sink
                                                               │
                                                               ▼
                                                          application
                                                        选择过滤与落盘策略
```

依赖规则如下：

1. `trace` 不依赖其他 pemu 模块；
2. `linalg` 公有依赖 `pemu::trace`，因为 `SolverResult` 可以携带统一的诊断事件；
3. `equation` 原样传播底层诊断，也可以为自己检测到的可恢复失败形成诊断；
4. `simulation` 公有依赖 `pemu::trace`，消费返回诊断并交给注入的 sink；
5. `field`、`physics`、`mesh` 和 `unit` 不因普通参数校验而依赖 trace；
6. 最终应用决定输出到 `std::clog`、内存、测试收集器或未来的异步后端。

这样不会出现底层数值算子偷偷访问全局 logger 的情况，也不会令一个线性求解器同时
承担数值计算和 I/O 策略。

## 3. 结构化事件模型

`TraceEvent` 包含：

- `kind`：`Trace` 表示运行轨迹，`Diagnostic` 表示需要关注的诊断；
- `domain`：由 `DiagDomain` 表示的来源模块，例如 `simulation`；
- `category`：模块内的子系统，例如 `fixed_step`；
- `name`：事件名称，例如 `step.completed`；
- `message`：可选的人类可读说明，普通 trace 通常为空；
- `severity`：`Trace`、`Debug`、`Info`、`Warning`、`Error` 或 `Critical`；
- `attributes`：由名称和强类型标量值构成的只读视图。

属性当前支持布尔值、有符号/无符号整数、`double` 和 `std::string_view`。事件及属性是
同步、非拥有视图：sink 必须在调用返回前消费它们；需要长期保存的 sink 必须复制名称
和字符串值。仿真侧用定长栈数组构造属性，不为每个事件分配堆内存。

`TraceSink` 是 concept，而不是虚函数接口：

```cpp
template <typename Sink>
concept TraceSink = requires(Sink& sink, const TraceEvent& event) {
  { sink(event) } noexcept -> std::same_as<void>;
};
```

sink 必须是 `noexcept`，诊断失败不能改变数值推进的控制流。`OstreamTraceSink` 捕获流
异常并记录自身的失败状态；其实现使用 fmt 在内部缓冲区中完成格式化，再写入目标流。
它在首次事件前输出一次表头，将 `step`、`time`、`dt`、`end_time`、求解器状态及残差
放入对齐列；一个事件缺少的属性对应空白单元格，其他属性和诊断消息进入 `DETAILS` 列。
为支持长时间仿真的 `tail -f` 预览，文件流默认在每两个 simulation 时间步完成后刷新；
`run.completed` 或 `run.failed` 会无条件刷新尚未刷新的尾部事件。此策略只刷新 C++ 流
缓冲区以便观察，不提供断电耐久性的 `fsync` 保证。
所有枚举的文本转换统一使用 `pemu::to_string()`；求解状态因此显示稳定字符串（如
`Success`、`SolveFailed`），
不输出 `enum class SolverStatus` 的底层整数。
fmt 只出现在实现文件中，并作为 `pemu::trace` 的私有构建依赖，不会泄漏到公共头文件。
自定义网络或文件 sink 也应在内部处理重试、丢弃或错误计数。

## 4. 结构化诊断模型

`diag.hpp` 将诊断域定义集中在一处。`DiagDomain` 的枚举项直接采用模块名：
`linalg`、`mesh`、`unit`、`field`、`trace`、`boundary`、`discretization`、`physics`、
`equation` 和 `simulation`。`diagDomainName()` 返回相同的稳定字符串，使 sink 可以按
模块过滤或生成机器可读记录。

诊断不定义第二种事件或第二套 sink，而是使用 `EventKind::Diagnostic` 标记
`TraceEvent`。其中 `category` 和 `name` 合起来构成稳定诊断码，例如：

```cpp
#include <pemu/trace/ostream_trace_sink.hpp>

pemu::trace::TraceEvent diagnostic{
    .kind = pemu::trace::EventKind::Diagnostic,
    .domain = pemu::trace::DiagDomain::equation,
    .category = "poisson",
    .name = "not_converged",
    .message = "linear solve did not converge",
    .severity = pemu::trace::Severity::Warning,
};

pemu::trace::OstreamTraceSink sink(std::clog);
sink(diagnostic);
```

它仍通过 `TraceSink` 输出为
`[Warning] equation.poisson.not_converged: linear solve did not converge`。应用可统一按照
`kind`、`domain` 和 `severity` 过滤，不需要管理两套 sink，也不会产生不一致的级别、
属性和输出策略。

### 4.1 返回值携带诊断的边界

`SolverResult` 可以携带一个可选的 `TraceEvent`。成功结果不携带诊断；失败结果由最先
明确知道根因的模块创建 `EventKind::Diagnostic`，上层必须原样传播，不能把
`linalg.cholmod.factorize.not_positive_definite` 抹平成含义更弱的
`simulation.step.failed`。

返回诊断使用 `makeDiagnosticEvent()` 构造。其域、类别、名称和消息必须引用静态存储，
返回时不附带属性 span；这样复制 `SolverResult` 不会产生悬空视图。simulation 消费时
再在栈上附加 `step`、`solver_status`、`residual_norm` 和 `relative_residual`，同步调用
sink。根因诊断之后仍可发出 `step.failed`，前者说明“为什么失败”，后者说明“哪个上层
阶段因此终止”。

### 4.2 什么情况下不应返回诊断

返回值附带诊断只用于调用者可能恢复、改用其他策略或正常终止的失败。构造参数非法、
越界访问、跨网格误用等编程契约错误仍抛出异常；它们不应为了统一形式而改成容易被忽略
的返回码。`linalg` 坚持 API-boundary / failure-path instrumentation：只在
`analyzePattern()`、`factorize()`、`solve()` 的入口检查和后端失败出口形成诊断，不在
矩阵遍历、分解或回代核心路径中发事件，也不持有 sink。

## 5. 与仿真的耦合

固定步长和自适应仿真都增加了第二个模板参数 `TraceSink`，默认值为
`pemu::trace::NullTraceSink`。sink 由构造函数注入并按值保存；有共享状态需求时，sink
可以内部持有引用或指针。默认 sink 是空类型，通过 `[[no_unique_address]]` 保存，不引入
虚调用，也不改变原有调用代码；`if constexpr` 同时跳过默认路径上的事件和属性构造。

例如，将自适应仿真的事件直接写到标准日志流：

```cpp
#include <pemu/trace/ostream_trace_sink.hpp>

pemu::trace::OstreamTraceSink trace(std::clog);

pemu::simulation::AdaptiveStepPlasmaSimulation simulation(
    density, reactions, transport, evaluator,
    pemu::simulation::AdaptiveTimeClock(end_time), trace,
    pemu::trace::StatisticsOptions{
        true, 10, 1.0 * cm, isq::length[cm]});

simulation.run();
```

输出为一行一个事件，例如：

```text
[Debug] simulation.adaptive_step.timestep.selected step=2 dt=0.02 transport_limit=inf positivity_limit=inf stability_limit=inf
```

当前事件顺序如下：

| 事件 | 发生位置 | 主要属性 |
| --- | --- | --- |
| `run.started` | `run()` 进入时 | `step`、`time`、`end_time` |
| `step.started` | 读取本步状态之前 | `step`、`time`，以及固定 `dt` 或自适应 `remaining_time` |
| `electrostatics.completed` | 电荷密度、泊松方程和电场更新成功后 | `solver_status`、残差 |
| `step.failed` | 电静力求解返回失败状态时 | `solver_status`、残差 |
| `reaction_rates.completed` | 反应率求值后 | 反应率场数 |
| `sources.completed` | 化学计量累积为物种源项后 | 源项场数 |
| `physics.species.statistics` | 输运更新前，每个已采样物种一次 | 极值、负值/非有限值计数、体积均值、粒子总数、RMS |
| `physics.charge.statistics` | 输运更新前 | 空间电荷极值、净电荷、绝对电荷及相对不平衡度 |
| `field.potential.statistics` | 输运更新前 | 电势极值、体积均值、RMS |
| `field.electric_field_normal.statistics` | 输运更新前 | 面法向电场极值、最大绝对值、面积均值、RMS |
| `timestep.selected` | 自适应推进完成步长选择与状态更新后 | `dt`、输运/正性/稳定性限制 |
| `step.completed` | 状态更新成功且时钟提交后 | 新 `step`、新 `time`、实际 `dt`、相对残差 |
| `run.completed` | 时钟到达终止条件后 | `step`、`time`、`end_time` |
| `run.failed` | `run()` 收到失败的 `SolverResult` 后 | `step`、`time`、`end_time` |

这里需要特别注意时间层：四类 statistics 事件位于电静力、反应率和源项均完成之后，
但在输运更新之前，因此全部描述同一个 $k$ 层状态。`step.completed` 中的密度已经是
$n^{k+1}$，而电势、电场、反应率和源项仍是用于本次推进的 $k$ 层工作量。

若电静力返回诊断，事件顺序为“底层根因 diagnostic → `step.failed`”。固定步长和自适应
仿真都遵循这一顺序，并且失败步不会提交时钟。

## 6. 物理体积与场统计

### 6.1 二维物理体积语义

二维网格接口中的 `cellVolume()` 实际返回单元面积 $A_c$，`faceArea()` 返回面长度
$\ell_f$。为使单位为 $\mathrm{cm}^{-3}$ 的数密度能够积分为粒子数，统计器把二维网格
解释成面外厚度为 $L_z$ 的平板：

$$
V_c=A_cL_z,\qquad A_f=\ell_fL_z.
$$

`StatisticsOptions` 的第三项是带 mp-units 单位的 $L_z$，第四项是网格坐标 reference；厘米
网格传入 `1.0 * cm, isq::length[cm]`。构造后只保存换算过的面外厚度裸值与桥接得到的
`precise_unit`。三维网格直接使用原生单元体积与面面积，并忽略厚度数值。事件中的
`volume_semantics=planar_extrusion` 或 `native_3d` 明示所用约定。

### 6.2 通用标量统计

`ScalarFieldStatisticsAccumulator` 逐样本累积，不分配与网格规模相关的临时数组，并用
补偿求和降低大范围数值相加的舍入误差。对值 $u_i$ 和物理权重 $w_i$，它计算

$$
I=\sum_iw_iu_i,\qquad
\bar u=\frac{I}{\sum_iw_i},\qquad
u_{\mathrm{rms}}=\sqrt{\frac{\sum_iw_iu_i^2}{\sum_iw_i}},
$$

以及最小值、最大值、最大绝对值、$L^1$ 积分、负值数、非有限值数和非法权重数。
species 的 `total_number` 是数密度的体积积分；charge 的 `net_charge` 是空间电荷密度
的体积积分，`relative_imbalance=|Q|/\int|\rho|\,dV`。

统计量同时携带运行期单位。若字段单位为 $[u]$，几何权重单位为 $[w]$，则 min、max、
mean 与 RMS 使用 $[u]$，权重和使用 $[w]$，integral 与 $L^1$ integral 使用
$[u][w]$。单位只在统计器收尾和 trace 事件组装阶段附加并相乘；遍历 cell/face 的热循环
仍只接收 `double value, double weight`。启用统计时，simulation 构造函数还会验证数密度、
电荷密度、电势和法向电场都具有期望的 `QuantityKind`，缺失或错误 metadata 会立即拒绝。

`TraceAttribute` 可附带 `precise_unit`，ostream 输出形如
`minimum=3.125 [V]`、`net_charge=0 [C]`。LLNL Units 可以把等价单位规范化显示，例如
`cm^3` 显示成 `mL`；无量纲量统一显示为 `[1]`。

正常统计作为 `Debug/Trace` 事件输出。物种密度出现负值时升级为 `Warning/Diagnostic`；
任意场出现 NaN、无穷值、非法物理权重或无有效统计量时升级为 `Error/Diagnostic`。

统计扫描的成本为每个采样步 $O(N_sN_c+N_f)$，因此默认关闭。启用后可用
`sample_every_steps` 独立控制采样频率；`NullTraceSink` 路径通过 `if constexpr` 完全
跳过事件和统计构造。

大规模场结果则走独立 output 管线。建议由 simulation 在“初始状态”“每隔若干步”及
“终止状态”提供一致的 snapshot 观察点，output 模块负责深拷贝或同步写出。日志过滤、
结果采样频率和数值推进步长应彼此独立。

## 7. 当前限制

- `OstreamTraceSink` 是同步输出，逐步打印大量事件会影响长时间仿真的吞吐；生产环境
  应使用带级别过滤、步数抽样或后台队列的自定义 sink；
- 当前 linalg 和泊松相容性检查会形成返回诊断；反应 evaluator 或输运更新直接抛出的
  异常仍原样向上传播，尚未定义可恢复的通用推进结果；
- 当前没有全局 logger、运行时 sink 注册表或跨线程排序，避免在单线程数值原型阶段
  提前固定并发模型。
