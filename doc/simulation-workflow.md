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

`start()` 只做同步状态转换；`advance()` 执行一个完整的电静力、反应、源项、输运和
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

## Checkpoint workflow

`CheckpointOptions` 将 checkpoint 调度纳入 simulation。每次成功提交物理状态后，
simulation 根据 `every_steps` 和 `write_final` 保存同一个可覆盖的 VTKHDF checkpoint；覆盖通过
同目录临时文件原子提交，写入失败不会破坏上一份成功状态。
checkpoint 写入失败会使 workflow 进入 `failed`。固定步和自适应步都自动生成稳定的物种
密度 key。自适应步还保存 `previous_dt`，使恢复后的 growth limiter 与未中断运行一致。

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
workflow schema、species ID/name、mesh、metadata 和时间校验全部成功后，才提交 density、clock
及时间积分历史。reaction network、边界条件及求解器配置仍由应用构造，这些外部配置必须与
checkpoint 对应的算例一致。
