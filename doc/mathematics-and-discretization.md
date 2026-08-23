# 数学模型、数值方法与离散化

本项目处理定义在二维区域 $\Omega$ 上的标量守恒问题。未知量可以是温度、浓度、电势或其他标量场；程序不限定其具体物理意义。它们的共同点是：一个小区域中某守恒量的变化，等于穿过边界的净通量与区域内源项的平衡。

项目采用单元中心有限体积法（finite-volume method, FVM）。它先在每个控制体上写积分守恒律，再近似控制体各面的通量；不是直接在节点上近似微分。因此同一条内面的通量对相邻两单元必为一进一出，离散后仍严格守恒。

## 1. 网格、记号与方向约定

将区域划分为不重叠控制体 $P$。当前为二维实现：

- $V_P$：控制体面积（程序接口名为 `cellVolume`）；
- $f$：控制体的一条边，$A_f$ 为其长度（接口名为 `faceArea`）；
- $\boldsymbol{x}_P$、$\boldsymbol{x}_f$：单元中心和面中心；
- $\boldsymbol{n}_f$：单位面法向；
- $\phi_P$ 或 $u_P$：保存在单元中心的数值未知量。

每条内面只存储一次，指定 owner 单元 $P$ 和 neighbor 单元 $N$，法向 $\boldsymbol{n}_f$ 从 $P$ 指向 $N$。边界面的法向从 owner 指向区域外。因而对 owner，正的法向通量代表流出；对 neighbor，完全相同的通量代表流入。这是所有守恒性质的共同基础。

两点通量近似使用沿面法向的投影距离：

$$
d_{PN} = (\boldsymbol{x}_N-\boldsymbol{x}_P)\cdot\boldsymbol{n}_f,\qquad
d_{Pb} = (\boldsymbol{x}_f-\boldsymbol{x}_P)\cdot\boldsymbol{n}_f.
$$

程序要求它们为正。若不为正，说明网格中的 owner/neighbor 关系、面法向或几何位置违反了上述约定。正交网格上，$d_{PN}$ 就是两个单元中心的间距。

## 2. 泊松方程：含义、边界与可解性

### 2.1 连续方程

项目中的稳态泊松（扩散）方程是

$$
-\nabla\cdot(\epsilon\nabla\phi)=\rho\qquad \text{in }\Omega.
$$

$\phi$ 为未知标量，$\epsilon\ge0$ 为扩散、导热或介电系数，$\rho$ 为单位面积源项。该式有多种物理解释：

- 稳态导热：$\phi=T$ 是温度，$\epsilon=k$ 是导热系数，$\rho$ 是热源；
- 静电学：$\phi$ 是电势，$\epsilon$ 是介电系数，$\rho$ 是电荷源（比例常数可并入 $\rho$）；
- 稳态物质扩散：$\phi$ 是浓度，$\epsilon=D$ 是扩散系数，$\rho$ 是产生或消耗率。

定义沿外法向的扩散通量密度

$$
q=-\epsilon\nabla\phi\cdot\boldsymbol{n}.
$$

负号就是 Fourier/Fick 定律：热量或物质从高 $\phi$ 向低 $\phi$ 扩散。对任一控制体 $P$ 积分，借助高斯定理得到

$$
\int_{\partial P}q\,\mathrm{d}s=\int_P\rho\,\mathrm{d}V.
$$

即“$\text{边界净流出量}=\text{内部产生量}$”。有限体积法离散的正是这条积分守恒律。

### 2.2 边界条件与符号

支持两类以区域外法向为准的边界条件。

1. Dirichlet 条件：直接指定场值

   $$
   \phi=\phi_b.
   $$

   例如指定墙温、电极电势或入口浓度。

2. Neumann 条件：直接指定外法向扩散通量

   $$
   -\epsilon\nabla\phi\cdot\boldsymbol{n}=q_b.
   $$

   程序的 Neumann `value` 就是 $q_b$。$q_b>0$ 表示量从区域流出，$q_b=0$ 表示绝热/无通量边界。

纯 Neumann 泊松问题只能确定 $\phi$ 的梯度；给解加任意常数仍是解，离散矩阵也因此
有常数零空间。此时构造 `PoissonSolver` 必须显式选择规范条件：`PinCellGauge` 固定
一个参考单元的电势，`ZeroMeanGauge` 则通过拉格朗日乘子施加体积加权零均值约束。
求解前还必须满足离散相容性

$$
\sum_P \rho_P V_P-\sum_{f\subset\partial\Omega}q_f A_f=0,
$$

也就是装配后 $\sum_P b_P=0$。求解器会按配置的绝对、相对容差检查该条件，不满足时
返回 `SolverStatus::IncompatibleRhs`。

## 3. 泊松方程的有限体积装配

### 3.1 每个控制体的一行方程

将积分边界切成若干面，用面中心通量近似每段积分：

$$
\sum_{f\subset\partial P}q_f A_f\approx\rho_P V_P.
$$

$\rho_P$ 是在单元中心采样的源项。这就是全局线性系统的一行：每个单元对应一个未知量 $\phi_P$，每条面决定这一行与自身及邻居的耦合。

### 3.2 内面：两点中心差分

对 owner $P$ 与 neighbor $N$ 之间的内面，法向梯度近似为

$$
\nabla\phi\cdot\boldsymbol{n}_f\approx\frac{\phi_N-\phi_P}{d_{PN}}.
$$

代入通量定义并乘面长：

$$
q_f A_f\approx
-\epsilon\frac{\phi_N-\phi_P}{d_{PN}}A_f
=c_f(\phi_P-\phi_N),
\qquad
c_f=\frac{\epsilon A_f}{d_{PN}}.
$$

因此该面给 $P$ 的方程添加 $+c_f\phi_P-c_f\phi_N$，给 $N$ 的方程添加 $+c_f\phi_N-c_f\phi_P$。写为矩阵局部块：

$$
\begin{bmatrix}
K_{PP}&K_{PN}\\
K_{NP}&K_{NN}
\end{bmatrix}
\mathrel{+}=
\begin{bmatrix}
c_f&-c_f\\
-c_f&c_f
\end{bmatrix}.
$$

这个式子同时解释三个重要事实：

- 内面通量在两个控制体之间严格抵消，所以守恒；
- 常系数扩散矩阵对称；
- 施加足够 Dirichlet 约束后，矩阵为正定，可使用 Cholesky 型求解器。

### 3.3 Dirichlet 面：已知边界值进入右端

边界面没有 neighbor。把面中心当作已知边界值位置，则

$$
\nabla\phi\cdot\boldsymbol{n}_f\approx\frac{\phi_b-\phi_P}{d_{Pb}}.
$$

记 $c_b=\epsilon A_f/d_{Pb}$，面流出量为

$$
q_fA_f\approx c_b(\phi_P-\phi_b).
$$

将未知量放在左端、已知量移到右端，装配规则是

$$
K_{PP}\mathrel{+}=c_b,\qquad b_P\mathrel{+}=c_b\phi_b.
$$

所以边界值不需要成为新的自由度，却会正确影响相邻单元的解。

### 3.4 Neumann 面：通量已经给定

若边界给定 $q_b$，面通量没有未知量。单元平衡为

$$
\sum_{\text{未知面}}q_fA_f+q_bA_f=\rho_PV_P.
$$

于是 Neumann 面不增加矩阵元素，只给右端项

$$
b_P\mathrel{+}=-q_bA_f.
$$

这说明符号的实际含义：正的 Neumann `value` 是向外流出，故会减少可由未知通量平衡的右端量。

### 3.5 两单元实例

`two_quads.msh` 是两个并排的 $1\times1$ 方形。取 $\epsilon=1$，全部边界为零 Dirichlet。共享内面的系数为 $1/1=1$；每条外边到单元中心的距离为 $1/2$，系数为 $1/(1/2)=2$。每个单元有三条外边和一条内边，故

$$
K=
\begin{bmatrix}
7&-1\\
-1&7
\end{bmatrix}.
$$

这正是 `AssemblesExpectedMatrixForTwoQuads` 用例检查的矩阵。若源项为 $\rho_0$、$\rho_1$，因两个面积均为 $1$，源项右端向量正好是 $(\rho_0,\rho_1)^\mathsf{T}$，此外再加边界贡献。

### 3.6 适用条件与空间精度

该近似称为两点通量近似（TPFA），本质是沿面法向的中心差分。对正交网格、常扩散系数和光滑解，它有二阶空间精度 $O(h^2)$：网格尺度减半时主导误差约为原来的四分之一。

强非正交网格上，单元中心连线不再与面法向平行，简单两点差分会漏掉切向梯度。一般需要非正交修正、多点通量近似（MPFA）或有限元法。当前实现没有这些修正，也不支持空间变系数、张量各向异性扩散或三维几何。

## 4. 显式扩散通量与离散散度

除了组装泊松矩阵，项目还提供可单独组合的面通量和散度算子。

### 4.1 diffusionFlux

给定单元场 $u$，`diffusionFlux` 在每个面存储

$$
\Gamma_f\cdot\boldsymbol{n}_f=-D\nabla u\cdot\boldsymbol{n}_f.
$$

内面和 Dirichlet 面使用与泊松装配相同的两点差分，但该算子不乘面长；Neumann 面直接使用给定通量。因而在线性场 $u=x$ 的正交网格上，数值通量精确为 $-D n_{f,x}$。

### 4.2 divergence 与离散高斯定理

`divergence` 将积分通量 $\Gamma_f\cdot\boldsymbol{n}_f A_f$ 加到 owner、从 neighbor 减去，最后除以单元面积：

$$
(\nabla\cdot\Gamma)_P\approx
\frac{1}{V_P}\sum_{f\subset\partial P}\Gamma_f\cdot\boldsymbol{n}_{P,f}A_f.
$$

全域求和时，每条内面精确抵消，因此满足离散高斯定理

$$
\sum_PV_P(\nabla\cdot\Gamma)_P
=\sum_{f\subset\partial\Omega}\Gamma_f\cdot\boldsymbol{n}_fA_f.
$$

这既是有限体积守恒性的代数表达，也是相关测试重点。

## 5. 对流方程与一阶迎风法

对流运输通常写为

$$
\frac{\partial u}{\partial t}+\nabla\cdot(vu)=s.
$$

当前项目提供其面通量算子。输入是面法向速度 $v_{n,f}=\boldsymbol{v}_f\cdot\boldsymbol{n}_f$；由于信息沿流线从上游传播，面状态必须取上游值：

$$
F_f=v_{n,f}u_{\mathrm{upwind}},\qquad
u_{\mathrm{upwind}}=
\begin{cases}
u_P,&v_n\ge0,\\
u_N,&v_n<0.
\end{cases}
$$

内面上，$v_n\ge0$ 表示 owner 流向 neighbor；$v_n<0$ 则反向。边界面上，$v_n\ge0$ 为出流，内部状态就是上游值；$v_n<0$ 为入流，必须给出 Dirichlet 外部状态。Neumann 扩散通量不含外部 $u$ 值，不能代替对流入流条件。

一阶迎风的优点是稳健，能避免高 $\mathrm{Pe}$ 数下中心插值常见的非物理振荡；代价是人为数值扩散，光滑区域通常只有一阶空间精度 $O(h)$。当前尚未把该算子与时间推进组合成完整的对流扩散方程求解器。

## 6. 瞬态扩散与后向欧拉

瞬态模块处理无体源扩散方程

$$
\frac{\partial u}{\partial t}-\nabla\cdot(D\nabla u)=0.
$$

对控制体积分，时间项为 $V_P\,\mathrm{d}u_P/\mathrm{d}t$，空间扩散项按第 3 节形成刚度矩阵 $K$。令 $M=\operatorname{diag}(V_P)$ 为对角质量矩阵，半离散方程为

$$
M\frac{\mathrm{d}u}{\mathrm{d}t}+Ku=b_{\mathrm{boundary}}.
$$

### 6.1 后向欧拉推导

令 $t^{n+1}=t^n+\Delta t$，并在新时刻评价扩散项：

$$
M\frac{u^{n+1}-u^n}{\Delta t}+Ku^{n+1}
=b_{\mathrm{boundary}}.
$$

整理为每一步所求的线性系统：

$$
\left(K+\frac{M}{\Delta t}\right)u^{n+1}
=\frac{M}{\Delta t}u^n+b_{\mathrm{boundary}}.
$$

这对应实际实现：矩阵的每个对角元加 $V_P/\Delta t$，旧时刻状态 $u^n$ 只进入右端向量。

后向欧拉是隐式方法。对扩散问题它无条件稳定：没有显式欧拉那种 $\Delta t$ 必须与 $h^2$ 同阶的稳定性限制；但仍有一阶时间误差 $O(\Delta t)$。零 Neumann 边界下，边界没有通量、内面通量成对抵消，故 $\sum_P V_Pu_P$ 每一步守恒；扩散同时会平滑空间起伏并降低单元间差异。

## 7. 稀疏线性系统与直接求解

泊松问题和每个后向欧拉时间步最终均为

$$
Ax=b.
$$

每个单元只和共享面的邻居耦合，所以 $A$ 是稀疏矩阵。稀疏存储避免保存大量零元，计算规模随网格邻接关系而非单元数平方增长。

求解器接口分三个阶段：

1. `analyzePattern(A)`：只分析非零元位置，确定重排与消元结构；
2. `factorize(A)`：对具体矩阵数值进行分解；
3. `solve(b,x)`：用已分解的因子求解。

三个阶段统一返回 `SolverResult`。这不是数值算法的一部分：求解器只在公共 API 的输入
检查以及后端报告失败的出口附加结构化 diagnostic，矩阵分析、分解和回代的核心计算中
不产生事件。方程层传播该结果，simulation 才负责把诊断交给外部 sink。

对带足够 Dirichlet 约束的对称正定扩散矩阵，`CholmodSolver` 使用 CHOLMOD 的 Cholesky 型分解，结构为 $A=LL^\mathsf{T}$。`UmfpackSolver` 使用 UMFPACK 的稀疏 LU 型分解，可处理一般非对称矩阵，适合未来含对流的系统。当前方程求解器会在首次调用时完成矩阵装配、模式分析与数值分解，之后只重建右端向量并复用分解；一旦网格、系数、边界或时间步改变，应重置或新建求解器。

## 8. 制造解、残差与收敛阶

测试不仅检查“线性方程能解”，还使用已知解析解验证离散是否正确。泊松制造解为

$$
\phi(x,y)=\sin(\pi x)\sin(\pi y),\qquad
-\Delta\phi=2\pi^2\sin(\pi x)\sin(\pi y),
$$

并在单位正方形边界施加零 Dirichlet 条件。瞬态扩散对应解析解为

$$
u(x,y,t)=e^{-2\pi^2Dt}\sin(\pi x)\sin(\pi y).
$$

误差使用面积加权离散 L² 范数：

$$
\lVert e\rVert_{L^2,h}=
\sqrt{\sum_P\left(u_P-u(x_P)\right)^2V_P}.
$$

相邻两层网格的观测空间阶为

$$
p=\frac{\log(e_h/e_{h/2})}{\log 2}.
$$

$p\approx2$ 说明二阶空间行为；固定空间网格、时间步逐半而终止时间相同时，$p\approx1$ 验证后向欧拉的一阶时间行为。代数残差 $\lVert A\boldsymbol{x}-\boldsymbol{b}\rVert/\lVert\boldsymbol{b}\rVert$ 检验线性求解器是否把离散系统解得足够准确；残差小不等于连续方程的离散误差小，因此测试分别检查残差和制造解误差。

## 9. 自适应等离子体推进的端到端算例与判据

当前的 `AdaptiveStepPlasmaSimulation` 已可用于在任意满足当前二维有限体积假设的
Gmsh 网格上推进多物种漂移扩散—泊松—反应系统。工程内置的较大网格文件名为
`data/meshfiles/poisson_64x64.msh`；它包含 $64\times64=4096$ 个控制体。求解结果
保存在调用方传入的 `SpeciesCellFields` 中，可通过 `simulation.density()` 读取；目前
可通过 `FieldOutputOptions` 接入 VTKHDF 后处理输出；它会写出同步的物种密度、电势和
电场等场数据，详见 `doc/output.md`。

### 9.1 64×64 均匀电离回归算例

`AdaptiveStepPlasmaSimulation64x64Test.SolvesUniformElectronImpactIonizationAcrossMultipleSteps`
在该网格上运行如下无量纲的五步算例。初始时电子、正离子密度均为

$$
n_e^0=n_i^0=1,
$$

四个电势边界均为零，物种 Dirichlet 边界也取 $1$。反应网络只含成对电离：

$$
\mathrm{e}+\mathrm{Ar}\longrightarrow2\mathrm{e}+\mathrm{Ar}^+.
$$

测试中的反应率模型为

$$
R=k n_e n_N,
\qquad
k=0.5,
\qquad
n_N=4,
$$

故初始反应率和两个带电物种的源项均为

$$
R=S_e=S_i=2.
$$

该算例将扩散系数取为严格正但极小的 $10^{-16}$，只是为了满足 SG 输运算子的参数
要求，同时使 $\Delta t=0.01$ 的多步算例中扩散误差远小于测试容差；它不是某种特定
气体的参数标定。`max_dt` 固定为 $0.01$，终止时间为 $0.05$，故 `run()` 必须完成
五个自适应步。电子、离子严格成对生成，因此每个电静力更新时都满足

$$
\rho=-n_e+n_i=0,
\qquad
\phi=0,
\qquad
\boldsymbol{E}=0.
$$

由于每一步的电离率随电子密度重算，均匀离散参考解满足

$$
n^{k+1}=\left(1+\Delta t k n_N\right)n^k=1.02n^k.
$$

因此五步后的预期结果为

$$
n_e^5=n_i^5=1.02^5\approx1.1040808032.
$$

对应测试逐单元和全域积分检查上述结果，因此它验证了网格读取、泊松求解、反应率、
化学计量源项、多步自适应时间推进和状态访问这一整条链路确实可在 4096 单元网格上运行。

### 9.2 64×64 厘米网格上的 400 V 平行板电离算例

`AdaptiveStepPlasmaSimulation64x64Test.ParallelPlate400VDrivesOppositeDriftAndIonizationInCentimeterMesh`
为网格坐标赋予明确的厘米单位。因此内置网格的范围 $[0,1]\times[0,1]$ 表示
$1\,\mathrm{cm}\times1\,\mathrm{cm}$ 的二维截面。该二维有限体积模型的积分量按单位外法向深度理解。

左右两条边界是平行板电极：

$$
\phi(0,y)=0\,\mathrm{V},
\qquad
\phi(1\,\mathrm{cm},y)=400\,\mathrm{V}.
$$

上下边界取绝缘条件：

$$
\frac{\partial\phi}{\partial n}=0.
$$

初始时电子与氩离子均匀且等密度：

$$
n_e^0=n_i^0=10^6\,\mathrm{cm^{-3}}.
$$

电荷、介电常数、迁移率、扩散系数、密度、反应系数和时间首先以 mp-units quantity
声明；只在进入求解器的边界通过 `numerical_value_in(...)` 转换为厘米制裸数值。场本身
继续存储 `double`，并通过 `FieldMetadata` 记录存储单位。电荷和介电常数为：

$$
q_e=-1.602176634\times10^{-19}\,\mathrm{C},
\qquad
q_i=+1.602176634\times10^{-19}\,\mathrm{C},
$$

$$
\varepsilon_0=8.8541878128\times10^{-14}\,\mathrm{F/cm}.
$$

因而初始净电荷为零，泊松方程的解析解是

$$
\phi(x,y)=400\frac{x}{1\,\mathrm{cm}},
\qquad
E_x=-400\,\mathrm{V/cm}.
$$

测试首先逐控制体检查该线性电势，并在左右边界面检查外法向电场。输运参数为

$$
\mu_e=10^3\,\mathrm{cm^2/(V\,s)},
\quad
D_e=10^2\,\mathrm{cm^2/s},
$$

$$
\mu_i=1.5\,\mathrm{cm^2/(V\,s)},
\quad
D_i=4\times10^{-2}\,\mathrm{cm^2/s}.
$$

化学模型仍是给定系数的局部成对电离：

$$
R_{\mathrm{ion}}=k_{\mathrm{ion}}n_en_N,
\qquad
n_N=2.5\times10^{19}\,\mathrm{cm^{-3}},
\qquad
k_{\mathrm{ion}}=10^{-13}\,\mathrm{cm^3/s}.
$$

物种边界均取初始密度的 Dirichlet 条件，以保持原 400 V 回归算例的边界定义不变。
Scharfetter--Gummel 物种输运现在也接受 Neumann：其值指定外向扩散通量，总通量还会
保留边界漂移出流；因此齐次 Neumann 在零法向漂移时才等同于严格 zero-flux。仿真运行到

$$
t_{\mathrm{end}}=2\times10^{-7}\,\mathrm{s},
$$

并将 $\Delta t_{\max}=10^{-6}\,\mathrm{s}$ 设得远大于实际漂移稳定性上限。因此测试要求
不止一个时间步，且末步小于该上限，证明步长并非由人为的 `max_dt` 固定。它还检查每个
控制体的反应率和两个成对源项为正，所有密度有限且非负；全域电子、离子数量均高于初始值，
并且电子质心在离子质心右侧。

最后一项不能简化为“离子质心必在初始位置左侧”。电离源与电子密度成正比；电子向右漂移后，
右侧会生成更多离子，因此离子质心可能也右移。电子与离子的相对质心位置才是该耦合算例中
稳健的输运方向判据。

该测试说明数值变量在厘米制下可量纲自洽地联用，但不把固定的
$k_{\mathrm{ion}}$ 宣称为真实气体放电的碰撞模型；后者仍需要 $E/N$ 依赖的反应系数、
气体状态和更完整的电极/鞘层边界模型。

### 9.3 “求解成功”不等于“物理模型已验证”

对一次数值仿真，建议按以下层次判定结果。

1. **代数求解成功。** 每次 `advance()` 返回的 `SolverResult` 必须满足
   `success()`，这表明泊松线性后端完成了模式分析、分解和求解，未报告奇异、非正定
   或求解失败。当前 CHOLMOD 与 UMFPACK 后端尚未填充 `residual_norm` 和
   `relative_residual`，因此它们的默认零值不能单独当作残差合格证据；若需要运行时
   残差，应由后端补充计算，或在调用端显式计算 $\lVert A\boldsymbol{x}-\boldsymbol{b}\rVert$。

2. **时间推进成功。** 仿真应到达目标时间，即 `finished()` 为真、
   $t=t_{\mathrm{end}}$，每个已接受步满足 $\Delta t>0$。自适应步长还应满足输运稳定、
   正性、`max_dt`、增长限制和最后剩余时间的共同约束。

3. **状态保持数值可接受。** 所有密度、电势和电场分量应有限；带输运物种密度不应为
   负。对存在解析守恒关系的问题，还应检查总量、净电荷或边界通量平衡。例如上述成对
   电离算例要求 $\int_\Omega n_e\,\mathrm{d}V=\int_\Omega n_i\,\mathrm{d}V=1.02^5$ 且
   $\rho=0$。

4. **离散模型正确且足够精确。** 这需要比一次成功运行更强的证据：制造解、网格加密
   收敛、时间步加密收敛，以及与实验或可信基准的比较。项目已有泊松与瞬态扩散的制造
   解/收敛测试；64×64 电离用例属于端到端回归与功能验证，不构成非均匀放电的物理验证。
