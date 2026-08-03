# eulerian-SPH

基于 [SPHinXsys](https://github.com/Xiangyu-Hu/SPHinXsys) 框架的欧拉 SPH 求解器，覆盖 2D 与 3D 弱可压缩与可压缩流动验证算例。

## 项目简介

`eulerian-SPH` 是 SPHinXsys 上游（`my_SPHinXsys`）的欧拉 SPH 专用分支。与一般 Lagrangian SPH 不同，欧拉 SPH 在固定背景网格上更新粒子场，借鉴有限体积法的 Riemann 解（HLLC + MUSCL 重构）处理对流项，从而：

- 在大变形、强间断（如激波）下保持稳定；
- 直接输出密度、压力、速度等守恒量，便于与网格参考解对照；
- 同时支持弱可压缩（low-Mach）与可压缩（high-Mach）两种流体模型。

本仓库聚焦于欧拉 SPH 相关算例与代码，包含：

- **基础算例**：Taylor–Green 涡（粘性可压缩）、圆柱绕流（Re=100 弱可压缩 / 可压缩 / 超声速 Ma=2）；
- **工业级算例**：3D 欧拉通道、跨声速压气机转子 Rotor 67；
- **统一回归测试**：基于时间平均 / DTW / 集合平均，可在没有 baseline 时自动跳过，避免单点失败阻断整个 case。

代码层面，所有算例都使用 **Laguerre–Gauss 核函数 + 线性梯度修正矩阵** 这套组合，并配合 Laguerre–Gauss 核的 `NoKernelCorrection` 选项控制壁面附近的核修正幅度，以抑制 contact 边界带来的奇偶振荡。

## 计算结果展示

### 2D 算例

**Lax 激波管（`test_2d_eulerian_shock_tube_LG`）**  
密度 / 压力 / 速度沿 x 方向与精确解（Riemann 解）对比，t=0.2：

![Lax shock tube](docs/eulerian_shock_tube.png)

可看到 SPH 解在激波、接触间断、稀疏波三段都能跟踪精确解，壁面附近无系统性跳变。

**超声速圆柱绕流（`test_2d_eulerian_supersonic_flow_new_BC`，Ma=2）**  
弓形激波 + 尾迹反射激波结构清晰可辨：

![Supersonic cylinder flow](docs/2d_eulerian_supersonic_flow_new_BC.png)

### 3D 算例

**3D 欧拉通道（`test_3d_eulerian_channel`）**  
弱可压缩 + 周期性 + 壁面边界，最终时刻速度大小云图，剖面呈抛物线分布：

![3D Eulerian channel](docs/3d_eulerian_channel.png)

**3D 弱可压缩圆柱绕流（`test_3d_eulerian_flow_around_cylinder_LG`，Re=100）**  
绕流后 Mach 数云图，可见卡门涡街：

![3D flow around cylinder (incompressible)](docs/3d_eulerian_flow_around_cylinder.png)

**3D 可压缩圆柱绕流（`test_3d_eulerian_compressible_flow_around_cylinder_LG`，Re=100）**  
可压缩模型下的速度大小云图，能看到前驻点 + 后侧回流区：

![3D flow around cylinder (compressible)](docs/3d_eulerian_compressible_flow_around_cylinder.png)

> 完整结果输出在 `tests/<case>/output/` 下，VTP 格式可用 [ParaView](https://www.paraview.org) 打开。

## 算例说明

### 2D 算例

| 算例 | 物理问题 | 流体模型 | 边界条件 |
|------|----------|----------|----------|
| `test_2d_eulerian_taylor_green_LG` | Taylor–Green 涡，Re=100 | 可压缩 + 粘性 | 全周期 |
| `test_2d_eulerian_flow_around_cylinder_LG` | 圆柱绕流，Re=100 | 弱可压缩 + 粘性 | 非反射远场 |
| `test_2d_eulerian_supersonic_flow_new_BC` | 超声速圆柱绕流，Ma=2 | 可压缩（无粘）| 幽灵粒子（反射壁 + 远场）|
| `test_2d_eulerian_shock_tube_LG` | Lax 激波管 | 可压缩 + Riemann | 反射壁 |

### 3D 算例

| 算例 | 物理问题 | 流体模型 | 边界条件 |
|------|----------|----------|----------|
| `test_3d_eulerian_channel` | 欧拉通道流 | 弱可压缩 | 周期性 + 壁面 |
| `test_3d_eulerian_flow_around_cylinder_LG` | 3D 圆柱绕流，Re=100 | 弱可压缩 + 粘性 | 非反射远场 |
| `test_3d_eulerian_compressible_flow_around_cylinder_LG` | 3D 可压缩圆柱绕流 | 可压缩 + 粘性 | 远场 |
| `test_3d_eulerian_rotor67` | 跨声速压气机转子 Rotor 67 | 可压缩 | 转子边界条件 |

## 目录结构

```
eulerian-SPH/
├── CMakeLists.txt
├── docs/                       # 算例结果展示图（本 README 引用）
├── src/
│   ├── shared/                 # SPHinXsys 共享核心库
│   │   ├── adaptations/        # 自适应比
│   │   ├── bodies/             # SPH body
│   │   ├── body_relations/     # inner / contact / complex 关系
│   │   ├── common/             # 基础数据结构、容器、变量
│   │   ├── geometries/         # 几何形状（level set、复杂几何）
│   │   ├── include/            # sphinxsys.h 聚合头
│   │   ├── io_system/          # VTK / VTP / XML / 参数化
│   │   ├── kernels/            # Wendland / Cubic / Laguerre–Gauss 等
│   │   ├── materials/          # 弱可压缩流体、可压缩流体、固体…
│   │   ├── mesh_dynamics/      # Level Set、邻居方法
│   │   ├── meshes/             # CellLinkedList、SparseMeshField
│   │   ├── particle_dynamics/  # 流体 / 固体 / 通用 / 连续介质
│   │   ├── particle_generator/ # 粒子生成器
│   │   ├── particle_neighborhood/
│   │   ├── particles/          # 粒子基类与操作
│   │   ├── regression_test/    # 回归测试（TA / DTW / EnsembleAverage）
│   │   ├── simbody_sphinxsys/  # Simbody 封装（保留 FSI 所需 SpatialVec 等）
│   │   ├── sphinxsys_system/   # SPH 系统管理
│   │   └── tools/              # tinyxml2、XML 解析
│   ├── for_2D_build/           # 2D 构建专属
│   ├── for_3D_build/           # 3D 构建专属（SDF / STL / 极分解等）
│   └── 3rd_party/              # 内嵌 tinyxml2
└── tests/
    ├── 2d_examples/
    │   ├── test_2d_eulerian_taylor_green_LG/
    │   ├── test_2d_eulerian_flow_around_cylinder_LG/
    │   ├── test_2d_eulerian_supersonic_flow_new_BC/
    │   └── test_2d_eulerian_shock_tube_LG/
    └── 3d_examples/
        ├── test_3d_eulerian_channel/
        ├── test_3d_eulerian_flow_around_cylinder_LG/
        ├── test_3d_eulerian_compressible_flow_around_cylinder_LG/
        └── test_3d_eulerian_rotor67/
```

> **说明**：相对上游 `my_SPHinXsys`，本项目显式剔除了 Cell-linked List 加速后端（`shared_ck`），CMake 在 `src/shared/` 下做了 glob 防护，避免未来上游同步再次引入；Simbody 部分保留 `SimTK::SpatialVec` 等用于 FSI 与 3D 几何。

## 项目同步策略

本项目从上游 `my_SPHinXsys` 周期同步欧拉 SPH 相关代码：

- 上游 `eulerian_fluid_dynamics`、`eulerian_ghost_boundary`、`eulerian_open_boundary`、`domain_bounding`、`continuum_dynamics` 等同步到 `src/shared`；
- 每次同步后只保留 ESPH 入口路径需要的文件，明确剔除 `shared_ck`；
- 同步流程与冲突解决详见 `task_plan.md`。

---

## 安装指南

### 依赖

| 依赖 | 版本要求 | 说明 |
|------|----------|------|
| CMake | ≥ 3.20 | 构建系统 |
| GCC / Clang | C++20（GCC 10+，Clang 12+）| 编译器 |
| [Eigen3](https://eigen.tuxfamily.org) | ≥ 3.4 | 向量/矩阵（仅头文件） |
| [TBB](https://github.com/oneapi-src/oneTBB) | ≥ 2021 | 并行任务调度 |
| [Boost](https://www.boost.org) | ≥ 1.74 | 多边形几何与命令行选项（`geometry` + `program_options`，必须） |
| [spdlog](https://github.com/gabime/spdlog) | ≥ 1.11 | 日志系统（`io_environment`、`sphinxsys_system`，必须） |
| [fmt](https://github.com/fmtlib/fmt) | ≥ 9.0 | spdlog 外部格式化后端（spdlog 以 `SPDLOG_FMT_EXTERNAL=ON` 构建时必须） |
| [Simbody](https://github.com/simbody/simbody) | ≥ 3.7 | FSI 耦合与 3D 几何（`TriangleMeshShape`、`SimTK::SpatialVec` 等） |
| [GTest](https://github.com/google/googletest) | ≥ 1.10 | 部分测试用例（可选，CMake `FetchContent` 自动下载） |

### Docker 快速开始（推荐）

如果希望在一台新的 Linux 电脑上快速部署本项目，同时避免宿主机依赖冲突，推荐使用 Docker 作为"依赖环境容器"。本项目代码保存在宿主机本地目录，Docker 只负责提供隔离好的编译和运行依赖；后续修改本地代码后，可以直接在容器内重新编译运行，无需重新安装宿主机依赖。

#### 1. 安装 Docker（Linux）

如果新机器还没有 Docker，请先按照 Docker 官方文档安装 Docker Engine。

**Ubuntu / Debian：**

```bash
sudo apt-get update
sudo apt-get install -y ca-certificates curl
sudo install -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "${UBUNTU_CODENAME:-$VERSION_CODENAME}") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

sudo apt-get update
sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
sudo systemctl enable --now docker
sudo docker run hello-world
```

如果是 Debian，请将上面仓库地址中的 `ubuntu` 替换为 `debian`。

**CentOS / RHEL 系：**

```bash
sudo dnf -y install dnf-plugins-core
sudo dnf config-manager --add-repo https://download.docker.com/linux/centos/docker-ce.repo
sudo dnf install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
sudo systemctl enable --now docker
sudo docker run hello-world
```

如果希望当前用户不加 `sudo` 直接运行 Docker，可参考 [Docker 官方后置步骤](https://docs.docker.com/engine/install/linux-postinstall/)。

#### 2. 删除旧目录并重新克隆仓库

```bash
mkdir -p /home/DataBank/SPH_solver
cd /home/DataBank/SPH_solver
rm -rf /home/DataBank/SPH_solver/eulerian-SPH
git clone https://github.com/KIYOYOZU/eulerian-SPH eulerian-SPH
cd /home/DataBank/SPH_solver/eulerian-SPH
```

#### 3. 构建 Docker 镜像

```bash
docker build -t eulerian-sph:latest .
```

镜像构建完成后，容器内具备以下依赖：`cmake`、`g++ / build-essential`、`Eigen3`、`TBB`、`Boost::geometry`、`Boost::program_options`、`Ninja`。

#### 4. 启动开发容器

把宿主机仓库整体挂载进容器，构建产物和输出会直接落在本地仓库中，本地改代码后能立即重新编译。

```bash
cd /home/DataBank/SPH_solver/eulerian-SPH

docker run --rm -it \
  -v "$PWD":/workspace \
  -e ROOT_DIR=/workspace \
  -e BUILD_DIR=/workspace/build \
  --entrypoint bash \
  eulerian-sph:latest
```

容器内项目根目录对应为 `/workspace`。

#### 5. 在容器内编译并运行某个 test

以 Taylor-Green 涡为例：

```bash
cd /workspace

cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release

cmake --build build --target test_2d_eulerian_taylor_green_LG --parallel 8

cd tests/2d_examples/test_2d_eulerian_taylor_green_LG
../../../build/tests/2d_examples/test_2d_eulerian_taylor_green_LG/bin/Release/test_2d_eulerian_taylor_green_LG
```

> **注意**：MSVC 多配置生成器会在 `bin/` 下再建 `Release/` 子目录，完整路径为  
> `build/tests/2d_examples/<test_name>/bin/Release/<test_name>.exe`

如果之后在宿主机修改了 `src/`、`tests/` 或 `CMakeLists.txt`，只需回到容器内重新执行：

```bash
cd /workspace
cmake --build build --target <target_name> --parallel 8
```

如果改动了 CMake 配置或新增源文件，推荐重新执行：

```bash
cd /workspace
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --target <target_name> --parallel 8
```

#### 6. 粒子松弛与重载运行

圆柱绕流算例如果要先做粒子松弛，再加载松弛结果正式运行，可在容器内或通过 `docker run` 传参执行。

```bash
/workspace/build/tests/2d_examples/test_2d_eulerian_flow_around_cylinder_LG/bin/Release/test_2d_eulerian_flow_around_cylinder_LG --relax=true

/workspace/build/tests/2d_examples/test_2d_eulerian_flow_around_cylinder_LG/bin/Release/test_2d_eulerian_flow_around_cylinder_LG --reload=true
```

#### 7. 常用环境变量

| 变量 | 默认值 | 说明 |
|------|--------|------|
| `ROOT_DIR` | `/opt/eulerian-SPH` | 项目根目录，可改成挂载进去的源码目录 |
| `BUILD_DIR` | `$ROOT_DIR/build` | Docker 内使用的构建目录 |
| `BUILD_TYPE` | `Release` | CMake 构建类型 |
| `BUILD_JOBS` | `nproc` | 并行编译线程数 |
| `FORCE_CONFIGURE` | `0` | 设为 `1` 时强制重新执行 CMake 配置 |
| `CONFIGURE_ONLY` | `0` | 设为 `1` 时只配置不编译、不运行 |
| `SKIP_BUILD` | `0` | 设为 `1` 时跳过编译，仅运行已有二进制 |
| `RUN_ARGS` | 空 | 用环境变量传递简单运行参数 |

---

### 系统包管理器安装（有 root 权限）

```bash
# Ubuntu / Debian
sudo apt install cmake g++ libeigen3-dev libtbb-dev libboost-dev libboost-program-options-dev libspdlog-dev libfmt-dev libsimbody-dev libgtest-dev

# CentOS / RHEL
sudo dnf install cmake gcc-c++ eigen3-devel tbb-devel boost-devel boost-program-options spdlog-devel fmt-devel simbody-devel gtest-devel
```

### 超算 / 无 root 环境（手动编译安装）

超算节点通常无法使用包管理器，需将依赖编译到用户目录。以下以 `$HOME/deps` 为安装根目录为例。

**Eigen3**（纯头文件，解压即用）：

```bash
wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
tar -xzf eigen-3.4.0.tar.gz
cmake -S eigen-3.4.0 -B eigen-build \
      -DCMAKE_INSTALL_PREFIX=$HOME/deps/eigen3
cmake --build eigen-build --target install
```

**TBB**：

```bash
wget https://github.com/oneapi-src/oneTBB/archive/refs/tags/v2021.10.0.tar.gz
tar -xzf v2021.10.0.tar.gz
cmake -S oneTBB-2021.10.0 -B tbb-build \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=$HOME/deps/tbb \
      -DTBB_TEST=OFF
cmake --build tbb-build -j$(nproc)
cmake --build tbb-build --target install
```

**Boost**（仅需头文件，解压即用）：

```bash
wget https://boostorg.jfrog.io/artifactory/main/release/1.82.0/source/boost_1_82_0.tar.gz
tar -xzf boost_1_82_0.tar.gz
# 无需编译，直接指定解压目录即可
```

**fmt**（spdlog 外部格式化后端，必须先于 spdlog 安装）：

```bash
wget https://github.com/fmtlib/fmt/archive/refs/tags/10.2.1.tar.gz
tar -xzf 10.2.1.tar.gz
cmake -S fmt-10.2.1 -B fmt-build \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=$HOME/deps/fmt \
      -DFMT_TEST=OFF -DFMT_DOC=OFF
cmake --build fmt-build -j$(nproc)
cmake --build fmt-build --target install
```

**spdlog**（日志系统，需启用外部 fmt）：

```bash
wget https://github.com/gabime/spdlog/archive/refs/tags/v1.13.0.tar.gz
tar -xzf spdlog-1.13.0.tar.gz
cmake -S spdlog-1.13.0 -B spdlog-build \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=$HOME/deps/spdlog \
      -DSPDLOG_FMT_EXTERNAL=ON \
      -Dfmt_DIR=$HOME/deps/fmt/lib/cmake/fmt \
      -DSPDLOG_BUILD_EXAMPLE=OFF -DSPDLOG_BUILD_TESTS=OFF
cmake --build spdlog-build -j$(nproc)
cmake --build spdlog-build --target install
```

**Simbody**（FSI 与 3D 几何必需）：

```bash
wget https://github.com/simbody/simbody/archive/refs/tags/Simbody-3.7.tar.gz
tar -xzf Simbody-3.7.tar.gz
cmake -S simbody-Simbody-3.7 -B simbody-build \
      -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_INSTALL_PREFIX=$HOME/deps/simbody \
      -DBUILD_TESTING=OFF
cmake --build simbody-build -j$(nproc)
cmake --build simbody-build --target install
```

---

### Linux 编译

#### 第一步：配置

`-S` 指定源码根目录，`-B` 指定构建目录（编译中间文件和可执行文件均在此，不污染源码树）。

**有 root 权限（依赖由包管理器安装，cmake 自动查找）：**

```bash
git clone https://github.com/KIYOYOZU/eulerian-SPH.git
cd eulerian-SPH

cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release
```

**超算 / 无 root 环境（需显式指定依赖安装路径）：**

```bash
ROOT_DIR=$HOME/eulerian-SPH
BUILD_DIR=$ROOT_DIR/build
DEPS=$HOME/deps

cmake -S "$ROOT_DIR" -B "$BUILD_DIR" \
  -Wno-dev \
  -DCMAKE_BUILD_TYPE=Release \
  -DEigen3_DIR=$DEPS/eigen3/share/eigen3/cmake \
  -DTBB_DIR=$DEPS/tbb/lib64/cmake/TBB \
  -Dfmt_DIR=$DEPS/fmt/lib/cmake/fmt \
  -Dspdlog_DIR=$DEPS/spdlog/lib/cmake/spdlog \
  -DSimbody_DIR=$DEPS/simbody/lib/cmake/simbody \
  -DBOOST_ROOT=$DEPS/boost_1_82_0
```

#### 第二步：编译

```bash
# 2D 算例
cmake --build "$BUILD_DIR" --target test_2d_eulerian_taylor_green_LG              --config Release -j$(nproc)
cmake --build "$BUILD_DIR" --target test_2d_eulerian_flow_around_cylinder_LG      --config Release -j$(nproc)
cmake --build "$BUILD_DIR" --target test_2d_eulerian_supersonic_flow_new_BC       --config Release -j$(nproc)

# 3D 算例
cmake --build "$BUILD_DIR" --target test_3d_eulerian_channel                      --config Release -j$(nproc)
cmake --build "$BUILD_DIR" --target test_3d_eulerian_compressible_flow_around_cylinder_LG --config Release -j$(nproc)
cmake --build "$BUILD_DIR" --target test_3d_eulerian_flow_around_cylinder_LG      --config Release -j$(nproc)
cmake --build "$BUILD_DIR" --target test_3d_eulerian_rotor67                      --config Release -j$(nproc)

# 或一次编译全部算例
cmake --build "$BUILD_DIR" --config Release -j$(nproc)
```

#### 编译产物位置

```
build/tests/2d_examples/<test_name>/bin/Release/<test_name>
build/tests/3d_examples/<test_name>/bin/Release/<test_name>
```

---

### Windows 编译（vcpkg）

Windows 下推荐用 [vcpkg](https://github.com/microsoft/vcpkg) 管理依赖，需先安装 Visual Studio 2022 或 Visual Studio 2026（含 C++ 桌面开发组件）。

**第一步：用 vcpkg 安装依赖**

```powershell
cd D:\path\to\vcpkg
.\vcpkg install eigen3 tbb boost spdlog fmt simbody gtest --triplet x64-windows
```

**第二步：配置**

```powershell
cmake -S . -B build `
  -DCMAKE_TOOLCHAIN_FILE="$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake" `
  -DVCPKG_TARGET_TRIPLET=x64-windows `
  -DTBB_DIR="$VCPKG_ROOT/installed/x64-windows/share/tbb" `
  -DSimbody_DIR="$VCPKG_ROOT/installed/x64-windows/share/simbody" `
  -DBOOST_ROOT="$VCPKG_ROOT/installed/x64-windows"
```

> **注意**：如果使用 Visual Studio 18 2026 Build Tools，请改用 `-G "Visual Studio 18 2026" -A x64` 参数。

**第三步：编译**

```powershell
cmake --build build --target test_2d_eulerian_taylor_green_LG         --config Release
cmake --build build --target test_2d_eulerian_flow_around_cylinder_LG --config Release
cmake --build build --target test_2d_eulerian_supersonic_flow_new_BC  --config Release

cmake --build build --target test_3d_eulerian_channel                      --config Release
cmake --build build --target test_3d_eulerian_flow_around_cylinder_LG      --config Release
cmake --build build --target test_3d_eulerian_rotor67                      --config Release

# 或编译全部
cmake --build build --config Release
```

**第四步：运行**

```powershell
cd tests\2d_examples\test_2d_eulerian_taylor_green_LG
..\..\..\..\build\tests\2d_examples\test_2d_eulerian_taylor_green_LG\bin\Release\test_2d_eulerian_taylor_green_LG.exe
```

> **注意**：MSVC 多配置生成器会在 `bin\` 下再建 `Release\` 子目录，完整路径为  
> `build/tests/2d_examples/<test_name>/bin\Release\<test_name>.exe`

---

### 运行

程序运行时会在**当前工作目录**下写入 `output/`、`restart/`、`reload/` 等文件夹，因此必须先 `cd` 到对应 test 源目录再执行。

**超算环境还需提前设置动态库路径（TBB、Simbody 为动态库时）：**

```bash
export LD_LIBRARY_PATH=$HOME/deps/tbb/lib64:$HOME/deps/simbody/lib64:$LD_LIBRARY_PATH
```

```bash
# 2D Taylor-Green 涡
cd tests/2d_examples/test_2d_eulerian_taylor_green_LG
../../../../build/tests/2d_examples/test_2d_eulerian_taylor_green_LG/bin/Release/test_2d_eulerian_taylor_green_LG

# 2D 圆柱绕流（弱可压缩，Re=100）
cd tests/2d_examples/test_2d_eulerian_flow_around_cylinder_LG
../../../../build/tests/2d_examples/test_2d_eulerian_flow_around_cylinder_LG/bin/Release/test_2d_eulerian_flow_around_cylinder_LG

# 2D 超声速圆柱绕流（Ma=2）
cd tests/2d_examples/test_2d_eulerian_supersonic_flow_new_BC
../../../../build/tests/2d_examples/test_2d_eulerian_supersonic_flow_new_BC/bin/Release/test_2d_eulerian_supersonic_flow_new_BC

# 3D 欧拉通道
cd tests/3d_examples/test_3d_eulerian_channel
../../../../build/tests/3d_examples/test_3d_eulerian_channel/bin/Release/test_3d_eulerian_channel

# 3D 圆柱绕流
cd tests/3d_examples/test_3d_eulerian_flow_around_cylinder_LG
../../../../build/tests/3d_examples/test_3d_eulerian_flow_around_cylinder_LG/bin/Release/test_3d_eulerian_flow_around_cylinder_LG
```

结果写入 `output/`（VTP 格式），可用 [ParaView](https://www.paraview.org) 打开。