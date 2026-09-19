# particle_band_asr —— SPH-ASR 粒子带自适应加密模块

基于 Yang, Kong & Liu 2021 (PRE 104:055308) 的 SPH-ASR 思想，实现于本仓库的
Eulerian 多相 Godunov（Kapila 五方程）框架之上。2D，一阶精度，NO shifting。

## 论文方程 ↔ 代码对照

| 论文 | 含义 | 代码 |
|---|---|---|
| Eq. (3) | 对称核梯度 ½(∇W(h_i)+∇W(h_j)) | `ParticleBandInnerRelation` 缓存的 dW_ij（审计：bitwise 对称） |
| Eq. (5) | 双曲核（抑制拉伸不稳定） | `KernelHyperbolic`（ASR case 默认；WendlandC2 用于 SPHAdaptation 标定） |
| Eq. (6)-(9) | 变光滑长度 h_r=1.5V^{1/d}，夹在 [0.5h_r, 2h_r] | `UpdateSmoothingLengthByBand`（h 夹在 [0.5h_r, min(1.5h_r, h_ref)]） |
| Eq. (34)-(38) | 带间距比 C_r=2^{1/d}、带宽 ΔS=5Δs（2D） | `ParticleBandAdaptation`（ds_max_factor/band_coef/band_width_factor） |
| Eq. (37) 约束 | 邻居只跨 ±1 带 | 构建时 band 交互限制 + cross-band 审计 |
| Eq. (39)-(40) | 参考质量 m_r=ρ_r·Δs^d、γ=m/m_r | `ParticleSplittingByBand` / `ParticleMergingByBand` 判据 |
| Eq. (41)-(45) | γ_s=1.5 分裂、λ=0.6 偏移 | 同上（分裂方向取最近邻连线法向） |
| Eq. (50)-(53) | γ_m=0.7 合并 | `ParticleMergingByBand`（守恒合并 + EOS 恢复） |

本模块**不**包含论文的：弱可压 EOS (Eq.12)、人工黏性 (Eq.16)、粒子位移
技术、CSF 表面张力。动量/质量/能量输运由仓库既有的 Eulerian 多相积分器
（HLLC + 反射墙 ghost）承担，α 用 upwind 非守恒形式。

## 组件与适用边界

| 组件 | 作用 | 适用边界（实测，见 case findings.md） |
|---|---|---|
| `ParticleBandAdaptation` + `BandedLattice` | 分级带格子 + 初始布置 | 通用 |
| `UpdateParticleBands` | 特征跟踪分带：界面（等体积气水差指标 + 空bin跳过 + 滞后带 + Dijkstra 回退）；**激波**（`shock_band=on`：粒子级三重门种子——α 同纯相对 + div u<0 压缩 + 对压跳≥1% 局部压力——与界面种子合并走统一多种子 Dijkstra，按到最近特征的图距离分带） | 平面界面精确；复杂拓扑走 Dijkstra；激波跟踪为论文之外的方法论扩展（论文只跟相界面）：无 bin、坐标无关、3D 直推、多激波/无 argmax 竞争 |
| `UpdateSmoothingLengthByBand` | 局部 h 演化 | 通用 |
| `ParticleSplittingByBand` / `ParticleMergingByBand` | 分裂/合并 | 事件守恒 1e-13 级；x 墙向越界未处理（分裂 y 已 wrap） |
| `ComputeGradientCorrection` | **零阶一致核梯度修正** c_i=M_i/ΣV_j（内区+墙 stencil）：质量/动量/能量通量把 `dW_ij e_ij V_j` 换成 `(dW_ij e_ij − c_i) V_j`，使均匀态净通量恒零。只去常数偏置、保留真实梯度，**动态可用**（区别于背景压力修正）。α 差分项免疫 M_i 不修正。变量零初始化，均匀格子 c_i≈0 → 非 ASR 精确无操作 | 通用（分级格子带边界 + 镜像墙 cutoff 截断的 M_i≠0）；代价是守恒从逐位降为容差（~1e-6 量级） |
| `MultiphaseBackgroundPressureCorrection` | **静态 well-balancing**：抵消分级带边界一阶矩 M_i≠0 的伪通量（三阶段时序 + 界面 upwind 端点） | **仅静止/近均匀态**（静止保持、分裂/合并事件、平移界面）。动态激波管中必须 `steady_correction=off`：光滑流中逐对平均扣除会减掉物理散度本身（波系淬灭）。动态场改用 `ComputeGradientCorrection` |
| `ShepardDensityFilter` | 密度滤波 | 可选、默认 off。气水强瞬变（EOS α-spike）与滤波正反馈会爆炸；平静双气引入 2.1e-4 相对压力偏置 |
| `MultiphaseAcousticTimeStepSizeLocalH` | 局部 h 声速步长 | 通用 |

## 验证状态（tests/2d_examples/test_2d_eulerian_multiphase_shock_tube_asr）

- 静态分级格子：邻居对称/完备性审计、跨带限制、wall 镜像覆盖 —— PASS
- 静止分级：伪速度 4.8e-11 ≈ uniform 基线 —— PASS
- 分裂/合并：事件守恒 ~1e-13、19× 事后审计 —— PASS
- 平移界面：双气解析跟踪（0.07 带宽、守恒 exact）；气水 parity 门 —— PASS
- 双套激波管验收（三方 USR-fine/USR-coarse/ASR）：双气五判据全过
  （L2(u) ASR≤coarse）；气水守恒/稳健/形态/成本过、精度部分
  （激波宽化 + smear 伪速度，known limitation）—— PASS（带记录）
- 激波跟踪带（gaswater/twogas `shock_band=on`）：激波带连续跟随运动激波
  （gaswater x_shock→0.939 vs 精确 0.9439），激波区加密到 band 0，抹宽
  0.0715 ≈ 均匀细网格 0.070、粒子数 57%；quiescent 均匀场零误报 —— PASS

原 case（test_2d_eulerian_multiphase_shock_tube）零回归：4000 粒子、
u*=0.5531/p* 平台与项目前基线一致。
