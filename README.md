# Eulerian-SPH

一个轻量独立的 2D Eulerian-SPH (光滑粒子流体动力学) 求解器，从 [SPHinXsys](https://github.com/Xiangyu-Hu/SPHinXsys) 中提炼通用接口，仅依赖 Eigen3，无其他重型依赖。

## 功能特性

- **两种求解器**：弱可压缩 Eulerian-SPH（低速不可压近似）和全可压缩 Eulerian-SPH（含 HLLC Riemann 求解器）
- **多种边界条件**：周期边界、远场边界（Riemann 特征波分解）、壁面边界、Ghost 边界
- **Laguerre-Gauss 核函数**：支持查表加速
- **自由表面检测**：基于符号距离场的自由面与表面法向更新
- **OpenMP 并行**：邻域搜索、修正矩阵、主循环均支持可选并行加速
- **结构化输出**：流场写 VTP（ParaView 兼容），观测量写 DAT 时间序列

## 示例算例

| 算例 | 物理问题 | 求解器 | 边界条件 |
|------|----------|--------|----------|
| `taylor_green` | Taylor-Green 涡（验证可压缩 SPH 精度） | 可压缩 Eulerian | 双向周期 |
| `flow_around_cylinder_lg` | 圆柱绕流（卡门涡街） | 弱可压缩 Eulerian | 远场 + 壁面 |
| `supersonic_flow_new_bc` | 超声速绕流（激波捕捉） | 可压缩 Eulerian | Ghost + 壁面 |

## 目录结构

```
eulerian-SPH/
├── src/
│   ├── core/           # 粒子状态、网格、邻域搜索、仿真主循环
│   ├── dynamics/       # 弱可压/可压缩求解器、表面指示、核修正
│   ├── kernel/         # Laguerre-Gauss 核函数及查表实现
│   ├── boundary/       # 周期、远场、Ghost、壁面边界
│   ├── math/           # 数学工具
│   ├── io/             # VTP 流场输出、DAT 观测量记录
│   └── eulerian_sph.h  # 统一头文件入口
├── cases/
│   ├── taylor_green/
│   ├── flow_around_cylinder_lg/
│   └── supersonic_flow_new_bc/
└── CMakeLists.txt
```

## 依赖

| 依赖 | 版本要求 | 说明 |
|------|----------|------|
| C++ | C++17 | 必须 |
| Eigen3 | ≥ 3.3 | 必须，线性代数 |
| OpenMP | 可选 | 并行加速 |
| CMake | ≥ 3.20 | 构建系统 |

## 构建

```powershell
# 默认构建（Release，开启 OpenMP）
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release

# 不启用 OpenMP
cmake -S . -B build -DEULERIAN_SPH_ENABLE_OPENMP=OFF
cmake --build build --config Release

# 只构建库，不构建算例
cmake -S . -B build -DEULERIAN_SPH_BUILD_CASES=OFF
cmake --build build --config Release
```

## 运行

```powershell
# Taylor-Green 涡（周期可压缩）
./build/cases/taylor_green/Release/taylor_green_case.exe

# 圆柱绕流（弱可压缩 + 远场边界）
./build/cases/flow_around_cylinder_lg/Release/flow_around_cylinder_lg_case.exe

# 超声速绕流（可压缩 + Ghost 边界）
./build/cases/supersonic_flow_new_bc/Release/supersonic_flow_new_bc_case.exe
```

输出文件将写入各算例的 `cases/<case_name>/output/` 目录，包括：
- `*.vtp`：ParaView 可读流场快照
- `*_TotalKineticEnergy.dat` / `*_MaximumSpeed.dat`：时间历程观测量

## 可视化

使用 [ParaView](https://www.paraview.org/) 打开 `output/*.vtp` 文件即可可视化粒子流场（密度、速度、压力等）。

## 关系：与 SPHinXsys 的对应关系

本项目从 SPHinXsys 的以下三个示例中提炼而来，仅保留 Eulerian 框架的核心链路：

- `test_2d_eulerian_taylor_green_vortex`
- `test_2d_flow_around_cylinder_LG`
- `test_2d_eulerian_supersonic_flow_new_BC`

主要简化：去除 Lagrangian-SPH 模块、粒子自适应、并行框架（TBB）依赖，保留弱可压/可压缩两路 Eulerian 积分、Riemann 边界、表面法向与核修正的完整数值链路。

## 当前状态

- Taylor-Green 涡：编译通过，数值链路已验证
- 圆柱绕流：弱可压链路完整，正在对齐卡门涡街长时程稳定性
- 超声速绕流：Ghost 边界 split advance 已接入，步长与长时程推进仍在调优
