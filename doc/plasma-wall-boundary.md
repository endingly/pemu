# Plasma wall / electrode 边界模型

## 1. M15 的分层

M15 将壁面交换视为边界法向通量，而不是体反应源。电极电势继续使用已有的
`boundary::Dirichlet`；它已经能表达接地、固定偏压和平行板电极，不需要再定义一套 electrode
voltage 类型。新增的 `pemu::physics::wall` 只描述等离子体粒子和能量与材料表面的交换，随后由
equation 层把这些通量离散进边界控制体。

实现依次分为壁面通量契约、二次发射组装、物种连续性离散、电子能量离散、数值验证和
Simulation 接入。这样 wall physics 不依赖时间循环，也不会把 Scharfetter–Gummel 的网格细节
放进物理参数类型。

## 2. 线性壁面通量契约

对边界面外法向 $\mathbf n$ 和物种 $s$，采用

$$
\boldsymbol\Gamma_s\cdot\mathbf n
=v_{s,\mathrm{loss}}n_s-\Gamma_{s,\mathrm{in}}.
$$

$v_{s,\mathrm{loss}}\ge0$ 是有效损失速度，第一项表示吸收产生的向外通量；
$\Gamma_{s,\mathrm{in}}\ge0$ 表示壁面发射回等离子体，因此带负号。这个接口故意不把热速度、
反射率、Bohm 速度或 sheath closure 硬编码为唯一公式：调用者先根据所选模型计算有效损失
速度，离散层只处理明确的线性通量契约。这样可以在不改变连续性算子的情况下替换壁面闭合。

电子或其他带能粒子的边界能量通量使用同样结构：

$$
\boldsymbol\Gamma_{w,s}\cdot\mathbf n
=v_{w,s,\mathrm{loss}}w_s-Q_{s,\mathrm{in}}.
$$

粒子损失速度和能量损失速度彼此独立。它们的比值取决于电子能量分布和采用的壁面近似，不能
由粒子通量接口擅自推断。

## 3. 二次发射

入射物种 $i$ 的主损失通量为

$$
\Gamma_{i,\mathrm{loss}}
=v_{i,\mathrm{loss}}n_i.
$$

固定产额 $\gamma_{i\rightarrow s}$ 产生

$$
\Gamma_{s,\mathrm{in}}
=\sum_i\gamma_{i\rightarrow s}\Gamma_{i,\mathrm{loss}},
$$

若每个发射粒子的平均能量为 $\varepsilon_{i\rightarrow s}$，则

$$
Q_{s,\mathrm{in}}
=\sum_i\varepsilon_{i\rightarrow s}
\gamma_{i\rightarrow s}\Gamma_{i,\mathrm{loss}}.
$$

`WallFluxAssembler` 在构造时把 boundary ID、物种 ID 和面解析完成；每次 `evaluate()` 只遍历
实际二次发射通道，并用当前边界 owner cell 的密度刷新发射通量，不扫描全部物种和体网格，也
不进行内存分配。候选缓存仅按不同的 `(face, emitted_species)` 目标分配。所有候选通量验证成功后
才发布，因此一次失败不会留下半更新的 wall state；首次求值前查询完整粒子/能量通量会被拒绝，
避免把尚未组装的二次发射静默当作零。

显式稳定性中的壁面损失率统一由离散层计算
$\lambda_{w,P}=\sum_{f\in\partial P}v_{\mathrm{loss},f}A_f/V_P$。fixed 和 adaptive 路径复用同一
实现，且会检查同一 owner cell 上多个壁面的聚合结果，防止逐面有限但总和溢出。

该模型与常见流体等离子体壁面条件保持相同的“粒子损失减发射”结构。具体的热运动、迁移修正
和反射系数属于 $v_\mathrm{loss}$ 的上游闭合；当前阶段也不求解 dielectric surface charge、
thermionic emission 或 field emission。

## 4. 符号与单位

- 正通量始终表示从 owner cell 指向壁外；
- `particle_loss_velocity` 使用长度/时间，乘数密度后得到粒子面通量；
- `energy_loss_velocity` 使用长度/时间，乘能量密度后得到能量面通量；
- `yield` 是每个入射粒子产生的发射粒子数；
- `emitted_mean_energy` 与电子能量方程采用同一能量单位，当前厘米制算例通常为 eV。

M15 后续离散必须把 $v_\mathrm{loss}A/V$ 加入显式损失率和 positivity timestep，不能只在最终
面通量上覆盖数值，否则自适应时间步会遗漏壁面耗散。

## 5. 参考模型

本实现的通量分解参考了 Hagelaar、de Hoog 与 Kroesen 提出的气体放电流体边界条件，以及
工程等离子体模型中电子随机热损失、离子轰击二次电子发射和发射电子能量通量的常见写法。
本项目选择有效损失速度作为接口，是为了显式暴露模型假设，而不是把某一套热速度定义固定成
不可替换的常数。

- G. J. M. Hagelaar, F. J. de Hoog, G. M. W. Kroesen, “Boundary conditions
  in fluid models of gas discharges,” *Physical Review E* 62, 1452 (2000),
  DOI: <https://doi.org/10.1103/PhysRevE.62.1452>；
- COMSOL Plasma Module, “The Wall Boundary Condition,”
  <https://doc.comsol.com/6.3/doc/com.comsol.help.plasma/plasma_ug_drift_diffusion.07.22.html>。
