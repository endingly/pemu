# 单位系统：mp-units 与场元数据

## 1. 设计目标

本项目使用 mp-units 处理数值模块边界上的量纲和单位换算，但不把
`mp_units::quantity` 放入有限体积循环、稀疏矩阵或线性求解器。其基本分层是：

$$
\text{带单位的输入}
\xrightarrow{\text{mp-units 检查与换算}}
\texttt{Field<double> 中的约定单位数值}
\xrightarrow{\text{核心数值运算}}
\text{带单位的输出视图}.
$$

这样做同时满足两点：

- 配置、测试、I/O 和后处理处可以写出明确的物理单位，并由编译器检查换算是否合法；
- 核心代数仍处理连续的 `double` 数组，不改变 Eigen、SuiteSparse 和有限体积算子的标量类型。

mp-units 通过 vcpkg 引入，当前锁定基线提供 2.5.0。项目通过纯 interface 目标
`pemu::unit` 统一引入 `mp-units::mp-units`；`field` 通过 `pemu::unit` 使用单位能力，
而 `unit` 不依赖 `field`、`mesh` 或 `physics`。

## 2. `FieldMetadata`

`CellField<T>` 与 `FaceField<T>` 都拥有一个 `FieldMetadata`，其中可记录：

- 字段实例名称，例如 `electron density`；
- 可选的 `pemu::unit::PhysicalQuantityMetadata`，它对 mp-units 的 quantity specification 和 unit
  进行类型擦除，并缓存由 mp-units 自动生成的单位符号。

物理量不再以 `"electric potential"` 之类的字符串重复记录。构造元数据时直接传入完整的
mp-units reference：

```cpp
field::makeFieldMetadata("potential", isq::electric_potential[V]);
```

其中 `isq::electric_potential` 是 quantity specification，`V` 是 unit。表达式
`isq::electric_potential[V]` 由 mp-units 检查物理量与单位是否相容；元数据中的单位符号也由
`mp_units::unit_symbol` 生成，调用者没有第二份可以写错的物理量名称或单位符号。

没有指定元数据的旧构造方式仍然有效，因此无量纲制造解和旧测试不必伪装成有量纲问题。
元数据按“每个场一份”存储，不随单元或面重复；访问
`field[cell]`、`data()`、`span()` 和 `fill(double)` 时不会查询单位。

尤其有：

```cpp
static_assert(std::same_as<field::CellField<double>::value_type, double>);
```

这保证单位没有进入场的核心标量类型。

字段及其通用集合统一由 `field` 模块提供：`FieldSet<Field, Id>` 持有一组同类型场，
`CellFieldSet<T, Id>` 与 `FaceFieldSet<T, Id>` 分别是单元场和面场的便捷别名。集合负责
网格关联、统一初始化、元数据传播、批量 `fill` 和按稠密强类型 ID 访问。`physics` 中的
`SpeciesCellFields`、`SpeciesFaceFields` 与 `ReactionRateFields` 现在只是把物理 ID 绑定到
这些通用容器的语义别名，不再各自实现或存储场容器。

## 3. 带单位边界适配器

`quantity_io.hpp` 提供三个边界函数：

- `setQuantity`：将 mp-units 量转换为字段元数据指定单位下的裸数值；
- `fillQuantity`：用一个带单位的量填充整个字段；
- `quantityAt`：把某个裸数值重新包装为带单位的量。

例如，字段按厘米存储，但输入使用米：

```cpp
using namespace mp_units;
using namespace mp_units::si::unit_symbols;

constexpr auto centimetre_length = isq::length[cm];
field::CellField<double> x(
    mesh, 0.0, field::makeFieldMetadata("x", centimetre_length));

field::setQuantity(x, mesh::CellId{0}, 2.5 * m, centimetre_length);

// 核心存储值是 250.0，单位由 metadata 说明为 cm。
const double raw = x[0];

// 在边界重新获得 quantity，并可继续转换为 m。
const auto value = field::quantityAt(x, mesh::CellId{0}, centimetre_length);
const double metres = value.numerical_value_in(m);
```

调用边界适配器时必须给出字段实际的 mp-units reference。如果元数据缺失，或者请求的
quantity specification 或 unit 与元数据不一致，接口会抛出 `std::invalid_argument`。因此，
即便两个物理量恰好使用相同单位，也不能在边界处静默互换。单位转换发生在传入的 quantity
与 reference 的存储单位之间；核心数组中不会保存逐元素 quantity 对象。

## 4. 厘米制等离子体场约定

`field::centimetrePlasmaFieldMetadata()` 给出当前厘米网格算例的一组一致元数据。其类型
`field::PlasmaFieldMetadata` 也位于 `field` 模块，避免场的单位目录散落在物理和求解器模块：

| 字段 | 存储单位 |
| --- | --- |
| 粒子数密度 $n_s$ | $\mathrm{cm^{-3}}$ |
| 粒子源项 $S_s$ | $\mathrm{cm^{-3}\,s^{-1}}$ |
| 反应率 $R_r$ | $\mathrm{cm^{-3}\,s^{-1}}$ |
| 电荷密度 $\rho$ | $\mathrm{C\,cm^{-3}}$ |
| 电势 $\phi$ | $\mathrm{V}$ |
| 电场 $E_n$ | $\mathrm{V/cm}$ |
| 漂移速度 $v_n$ | $\mathrm{cm/s}$ |
| 输运损失率 $\lambda$ | $\mathrm{s^{-1}}$ |

ISQ 已定义的电荷密度、电势和频率直接使用 `mp_units::isq` 中的 quantity specification。
粒子数密度、粒子数密度变化率、反应率密度、法向电场和法向漂移速度没有完全对应的标量
ISQ 类型，因此在 `pemu::unit::plasma_quantity` 中定义为项目级 mp-units quantity
specifications；它们仍由 mp-units 负责量纲组合与单位相容性检查。

现代固定步长和自适应多物种漂移扩散推进器会把这组元数据传给其内部电荷密度、电势、
电场、漂移速度和工作场。仿真再把相同约定传给反应率场与物种源项场。调用方创建初始
`SpeciesCellFields` 时，应把 `number_density` 元数据传入，从而使整条仿真链路具有一致说明。

400 V 平行板测试中的物理常数和模型参数首先写成 mp-units quantity，例如

$$
\varepsilon_0=8.8541878128\times10^{-14}\,\mathrm{F/cm},
\qquad
\mu_e=10^3\,\mathrm{cm^2/(V\,s)},
$$

随后仅在构造现有求解器参数时通过 `numerical_value_in(...)` 提取厘米制裸值。这正是单位
系统和核心代数之间的边界。

## 5. 有意保留的边界

当前实现不在每次加法、乘法或散度计算中动态检查单位。核心算子仍依靠其数学契约，例如
扩散通量函数的输入必须已经使用同一套约定单位。原因是这些调用位于高频数值路径，且其
量纲关系已经由离散方程固定。

单位系统主要防止以下错误：

- 在配置阶段把 $\mathrm{m^2/(V\,s)}$ 的迁移率误当成 $\mathrm{cm^2/(V\,s)}$ 裸值；
- 给一个按厘米存储的场写入米数值而忘记比例换算；
- I/O 或后处理把 $\mathrm{V/cm}$ 标成 $\mathrm{V/m}$；
- 把速度字段误当作电场字段通过单位边界接口读取。

尚未覆盖的工作包括：网格自身的长度元数据、边界条件值的单位元数据、序列化单位信息，
以及为每一种物种分别保存更具体的字段名称。这些都可以沿同一边界适配方式扩展，无需改变
核心矩阵和场数据的表示。
