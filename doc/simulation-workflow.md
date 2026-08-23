# Simulation workflow

`pemu/simulation/workflow.hpp` 聚合固定步、自适应步、状态、dump 与 checkpoint 配置。
`FixedStepPlasmaSimulation` 与 `AdaptiveStepPlasmaSimulation` 是非模板的轻量 façade：公开头文件只保留
workflow API，反应率计算和 trace 回调分别以 `PlasmaReactionRateEvaluator`、`trace::AnyTraceSink`
注入；两个类分别实现在对应的 `.cpp` 中，状态转换、dump 生命周期和 checkpoint 调度由私有公共
workflow 组件统一实现。
simulation 是同步的轻量状态机，不创建后台线程，也不依赖第三方状态机库：

| 状态 | 含义 | 可用转换 |
| --- | --- | --- |
| `ready` | 已完成配置，尚未运行；允许恢复 checkpoint | `start()`、`advance()`、`run()`、`stop()` |
| `running` | 可以执行完整物理步 | `advance()`、`run()`、`pause()`、`stop()` |
| `paused` | 保留全部内存状态，等待继续 | `start()`、`run()`、`stop()` |
| `completed` | 已到达时间终点 | 只读查询及显式保存 checkpoint |
| `stopped` | 用户主动终止 | 终态 |
| `failed` | 数值推进或 workflow 任务失败 | 终态 |

`start()` 只做同步状态转换；`advance()` 执行一个完整的电静力、反应、物种/电子能量源项、输运和
时间提交事务。首次在 `ready` 状态调用 `advance()` 会自动 `start()`。`run()` 从 `ready`
启动或从 `paused` 恢复，然后同步运行到 `completed`；需要交互式控制时，调用方应使用
`start()` / `advance()` / `pause()`，在物理步之间协作式暂停。`stop()` 会关闭已经采集的
dump series，但不会把未达到时间终点的状态伪装成 final snapshot。

```cpp
#include <pemu/simulation/workflow.hpp>

simulation.start();
auto result = simulation.advance();
simulation.pause();

// 继续直到时间终点。
simulation.run();
```

## 必需的电子能量配置

两个 Simulation 构造函数都要求在 reaction evaluator 与 clock 之间传入
`ElectronEnergyConfiguration`，不提供关闭电子能量方程的旧重载。配置引用调用方持有的
电子能量密度场，并指定 electron species、能量边界、density floor、输运 closure factor 和
必填的 `ElectronEnergyAdditionalSourceEvaluator`。workflow 先以连续性方程相同的 SG 电子
粒子通量自动组装 $-\boldsymbol\Gamma_e\cdot\mathbf E$；additional evaluator 随后在同一
时间层的电场、反应率和物种源均已更新后写入碰撞或外部源。附加源目标场先填充 NaN，因此
漏写任何单元都会使当前 step 失败，而不会静默当作零源。

生产 workflow 默认要求 density、potential、electric field、drift velocity、energy 和 source
具有完整且代数一致的运行期 metadata；simulation 在构造阶段一次性比较实际字段与 transport
声明的 `QuantityKind + precise_unit`，并由电势单位到平均能量存储单位计算电场功换算系数。
纯数值单元测试若有意不携带 metadata，必须显式设置 `allow_unitless_raw_values=true`；此时约定
电势裸值为 V、平均电子能量裸值为 eV。该开关不能绕过部分或互相矛盾的 metadata。

固定步长在提交前先计算并验证全部物种与能量候选状态，任一物种的强负源失败都不会留下
部分 density/energy 更新。自适应步长将物种和能量的输运/正性限制分别
取最小值后统一选步。`electronEnergyDensity()`、`electronMeanEnergy()` 与
`electronEnergySource()` 提供只读 workflow 视图；能量密度仍是调用方传入的同一个状态对象。

## Checkpoint workflow

`CheckpointOptions` 将 checkpoint 调度纳入 simulation。每次成功提交物理状态后，
simulation 根据 `every_steps` 和 `write_final` 保存同一个可覆盖的 VTKHDF checkpoint；覆盖通过
同目录临时文件原子提交，写入失败不会破坏上一份成功状态。
checkpoint 写入失败会使 workflow 进入 `failed`。固定步和自适应步都自动生成稳定的物种
密度 key，并以 electron species ID/name 绑定电子能量密度。自适应步还保存 `previous_dt`，
使恢复后的 growth limiter 与未中断运行一致。

```cpp
output::checkpoint::VtkHdfWriter checkpoint_writer;

simulation::CheckpointOptions checkpoint_options{
    .writer = &checkpoint_writer,
    .path = "checkpoint/plasma.vtkhdf",
    .every_steps = 100,
    .write_final = true,
    .overwrite = true,
};

// checkpoint_options 是 Fixed/AdaptiveStepPlasmaSimulation 构造函数的末个参数。
```

恢复由 simulation 自己组装字段目标、校验时间并重建 clock：

```cpp
output::checkpoint::VtkHdfReader checkpoint_reader;
simulation.restoreCheckpoint(checkpoint_reader, "checkpoint/plasma.vtkhdf");
simulation.run();
```

`restoreCheckpoint()` 只允许在 `ready` 状态调用。reader 只加载一次文件并先恢复到临时状态；
workflow schema、species ID/name、mesh、metadata、电子能量/平均能量和时间校验全部成功后，
才提交 density、energy、clock 及时间积分历史。schema v2 的 checkpoint 必含电子能量场；
reaction network、边界条件、能量源 evaluator 及求解器配置仍由应用构造，这些外部配置必须与
checkpoint 对应的算例一致。
