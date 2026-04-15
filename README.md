# eulerian-SPH

基于 [SPHinXsys](https://github.com/Xiangyu-Hu/SPHinXsys) 框架的轻量级 2D 欧拉 SPH 求解器，包含弱可压缩与可压缩流动验证算例。

---

## Docker 快速开始

如果希望在一台新的 Linux 电脑上快速部署本项目，同时避免宿主机依赖冲突，推荐使用 Docker 作为“依赖环境容器”。本项目代码保存在宿主机本地目录，Docker 只负责提供隔离好的编译和运行依赖；后续修改本地代码后，可以直接在容器内重新编译运行，无需重新安装宿主机依赖。

### 1. 删除旧目录并重新克隆仓库

以下命令示例假定项目统一放在 `/home/DataBank/SPH_solver/eulerian-SPH`。如果该目录已经存在旧版本仓库，先删除后重新克隆。下面使用实际仓库地址 `https://github.com/KIYOYOZU/eulerian-SPH`。

```bash
mkdir -p /home/DataBank/SPH_solver
cd /home/DataBank/SPH_solver
rm -rf /home/DataBank/SPH_solver/eulerian-SPH
git clone https://github.com/KIYOYOZU/eulerian-SPH eulerian-SPH
cd /home/DataBank/SPH_solver/eulerian-SPH
```

### 2. 构建 Docker 镜像

在仓库根目录执行：

```bash
docker build -t eulerian-sph:latest .
```

镜像构建完成后，容器内会具备以下依赖环境：

- `cmake`
- `g++ / build-essential`
- `Eigen3`
- `TBB`
- `Boost::geometry`
- `Boost::program_options`
- `Ninja`

### 3. 启动开发容器

推荐将宿主机本地仓库整体挂载到容器内，这样容器里的编译目录和运行输出会直接落在本地仓库中，本地改代码后也能立即重新编译。这里 Docker 只提供依赖环境，不保存开发中的源码状态。

```bash
cd /home/DataBank/SPH_solver/eulerian-SPH

docker run --rm -it \
  -v "$PWD":/workspace \
  -e ROOT_DIR=/workspace \
  -e BUILD_DIR=/workspace/build \
  --entrypoint bash \
  eulerian-sph:latest
```

进入容器后，项目根目录对应为：

```bash
/workspace
```

这样容器内构建产生的文件会直接落在宿主机本地仓库：

- `/home/DataBank/SPH_solver/eulerian-SPH/build`
- `/home/DataBank/SPH_solver/eulerian-SPH/cases/2d_eulerian_flow_around_cylinder_LG/output`

### 4. 在容器内编译并运行某个 case

以圆柱绕流算例为例，在容器内执行：

```bash
cd /workspace

cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release

cmake --build build --target 2d_eulerian_flow_around_cylinder_LG --parallel 8

cd /workspace/cases/2d_eulerian_flow_around_cylinder_LG

/workspace/build/cases/2d_eulerian_flow_around_cylinder_LG/bin/2d_eulerian_flow_around_cylinder_LG
```

上面这些命令虽然是在 Docker 容器里执行，但由于 `/home/DataBank/SPH_solver/eulerian-SPH` 已经挂载到容器内的 `/workspace`，因此结果会直接保存在宿主机本地：

- 编译结果：`/home/DataBank/SPH_solver/eulerian-SPH/build`
- 运行输出：`/home/DataBank/SPH_solver/eulerian-SPH/cases/2d_eulerian_flow_around_cylinder_LG/output`
- 其他运行目录：`reload/`、`restart/`、`logs/` 也会保存在对应本地 case 目录下

如果之后在宿主机修改了 `src/`、`cases/` 或 `CMakeLists.txt`，只需要回到容器内重新执行：

```bash
cd /workspace
cmake --build build --target 2d_eulerian_flow_around_cylinder_LG --parallel 8
```

如果改动了 CMake 配置或新增源文件，推荐重新执行：

```bash
cd /workspace
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build --target 2d_eulerian_flow_around_cylinder_LG --parallel 8
```

### 5. 粒子松弛与重载运行

圆柱绕流算例如果要先做粒子松弛，再加载松弛结果正式运行，可在容器内或通过 `docker run` 传参执行。

在已经进入容器的前提下，粒子松弛运行：

```bash
/workspace/build/cases/2d_eulerian_flow_around_cylinder_LG/bin/2d_eulerian_flow_around_cylinder_LG --relax=true
```

加载松弛结果继续运行：

```bash
/workspace/build/cases/2d_eulerian_flow_around_cylinder_LG/bin/2d_eulerian_flow_around_cylinder_LG --reload=true
```

### 6. 常用环境变量

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

## 依赖

| 依赖 | 版本要求 | 说明 |
|------|----------|------|
| CMake | ≥ 3.20 | 构建系统 |
| GCC / Clang | C++20（GCC 10+，Clang 12+）| 编译器 |
| [Eigen3](https://eigen.tuxfamily.org) | ≥ 3.4 | 向量/矩阵运算（仅头文件） |
| [TBB](https://github.com/oneapi-src/oneTBB) | ≥ 2021 | 并行任务调度 |
| [Boost](https://www.boost.org) | ≥ 1.74 | 多边形几何域构造与命令行选项解析（`geometry` + `program_options`，必须） |

### 系统包管理器安装（有 root 权限）

```bash
# Ubuntu / Debian
sudo apt install cmake g++ libeigen3-dev libtbb-dev libboost-dev libboost-program-options-dev

# CentOS / RHEL
sudo dnf install cmake gcc-c++ eigen3-devel tbb-devel boost-devel boost-program-options
```

### 超算 / 无 root 环境（手动编译安装）

超算节点通常无法使用包管理器，需将依赖编译到用户目录。以下以 `$HOME/deps` 为安装根目录为例：

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

---

## Linux 编译

### 第一步：配置（cmake -S ... -B ...）

`-S` 指定源码根目录（含根 `CMakeLists.txt`），`-B` 指定构建目录（编译中间文件和可执行文件均在此，不污染源码树）。

**有 root 权限（依赖由包管理器安装，cmake 自动查找）：**
```bash
git clone https://github.com/<your-name>/eulerian-SPH.git
cd eulerian-SPH

cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release
```

**超算 / 无 root 环境（需显式指定依赖安装路径）：**
```bash
# 根据实际安装位置修改以下变量
ROOT_DIR=$HOME/eulerian-SPH   # 本仓库根目录
BUILD_DIR=$ROOT_DIR/build     # 构建输出目录
DEPS=$HOME/deps               # 依赖库安装根目录

# *_DIR 变量指向各依赖的 cmake 配置文件所在目录（含 *Config.cmake）
# BOOST_ROOT 指向 Boost 解压根目录（含 include/boost/）
cmake -S "$ROOT_DIR" -B "$BUILD_DIR" \
  -Wno-dev \
  -DCMAKE_BUILD_TYPE=Release \
  -DEigen3_DIR=$DEPS/eigen3/share/eigen3/cmake \
  -DTBB_DIR=$DEPS/tbb/lib64/cmake/TBB \
  -DBOOST_ROOT=$DEPS/boost_1_82_0
```

### 第二步：编译（cmake --build --target）

`--target` 指定只编译某一个算例（不加则编译全部），`--config Release` 在多配置生成器下指定 Release 模式，`-j` 并行编译加速。

```bash
# 编译单个算例（推荐，节省时间）
cmake --build "$BUILD_DIR" --target 2d_eulerian_taylor_green_LG          --config Release -j$(nproc)
cmake --build "$BUILD_DIR" --target 2d_eulerian_flow_around_cylinder_LG  --config Release -j$(nproc)
cmake --build "$BUILD_DIR" --target 2d_eulerian_supersonic_flow_around_cylinder --config Release -j$(nproc)

# 或一次编译全部算例
cmake --build "$BUILD_DIR" --config Release -j$(nproc)
```

### 编译产物位置

```
build/cases/<case_name>/bin/<case_name>

# 例如：
build/cases/2d_eulerian_taylor_green_LG/bin/2d_eulerian_taylor_green_LG
build/cases/2d_eulerian_flow_around_cylinder_LG/bin/2d_eulerian_flow_around_cylinder_LG
build/cases/2d_eulerian_supersonic_flow_around_cylinder/bin/2d_eulerian_supersonic_flow_around_cylinder
```

---

## Windows 编译（vcpkg）

Windows 下推荐用 [vcpkg](https://github.com/microsoft/vcpkg) 管理依赖，需先安装 Visual Studio 2022（含 C++ 桌面开发组件）。

**第一步：用 vcpkg 安装依赖**

```powershell
# 安装三个依赖（x64-windows）
cd D:\path\to\vcpkg   # 替换为你的 vcpkg 根目录
.\vcpkg install eigen3 tbb boost --triplet x64-windows
```

**第二步：配置**

```powershell
# VCPKG_ROOT 替换为你的 vcpkg 安装路径，例如 D:\code\vcpkg
cmake -S . -B build `
  -DCMAKE_TOOLCHAIN_FILE="$VCPKG_ROOT/scripts/buildsystems/vcpkg.cmake" `
  -DVCPKG_TARGET_TRIPLET=x64-windows `
  -DTBB_DIR="$VCPKG_ROOT/installed/x64-windows/share/tbb" `
  -DBOOST_ROOT="$VCPKG_ROOT/installed/x64-windows"
```

**第三步：编译**

```powershell
# 编译单个算例
cmake --build build --target 2d_eulerian_taylor_green_LG         --config Release
cmake --build build --target 2d_eulerian_flow_around_cylinder_LG --config Release
cmake --build build --target 2d_eulerian_supersonic_flow_around_cylinder --config Release

# 或编译全部
cmake --build build --config Release
```

**第四步：运行**

```powershell
# 切换到 case 源目录后执行，output 将写入当前目录
cd cases\2d_eulerian_taylor_green_LG
..\..\build\cases\2d_eulerian_taylor_green_LG\bin\Release\2d_eulerian_taylor_green_LG.exe
```

> **注意**：MSVC 多配置生成器会在 `bin\` 下再建 `Release\` 子目录，完整路径为  
> `build/cases/<case_name>/bin/Release/<case_name>.exe`

---

## 运行

程序运行时会在**当前工作目录**下写入 `output/`、`restart/`、`reload/` 等文件夹，因此必须先 `cd` 到对应 case 源目录再执行。

**超算环境还需提前设置动态库路径（TBB 为动态库时）：**
```bash
export LD_LIBRARY_PATH=$HOME/deps/tbb/lib64:$LD_LIBRARY_PATH
```

```bash
# Taylor-Green 涡
cd cases/2d_eulerian_taylor_green_LG            # 切换到 case 源目录，output 将写入此处
../../build/cases/2d_eulerian_taylor_green_LG/bin/2d_eulerian_taylor_green_LG

# 圆柱绕流（弱可压缩，Re=100）
cd cases/2d_eulerian_flow_around_cylinder_LG
../../build/cases/2d_eulerian_flow_around_cylinder_LG/bin/2d_eulerian_flow_around_cylinder_LG

# 超声速圆柱绕流（Ma=2）
cd cases/2d_eulerian_supersonic_flow_around_cylinder
../../build/cases/2d_eulerian_supersonic_flow_around_cylinder/bin/2d_eulerian_supersonic_flow_around_cylinder
```

结果写入 `output/`（VTP 格式），可用 [ParaView](https://www.paraview.org) 打开。

## 算例说明

| 算例 | 物理问题 | 流体模型 | 边界条件 |
|------|----------|----------|----------|
| `2d_eulerian_taylor_green_LG` | Taylor-Green 涡，Re=100 | 可压缩 + 粘性 | 全周期 |
| `2d_eulerian_flow_around_cylinder_LG` | 圆柱绕流，Re=100 | 弱可压缩 + 粘性 | 非反射远场 |
| `2d_eulerian_supersonic_flow_around_cylinder` | 超声速圆柱绕流，Ma=2 | 可压缩（无粘）| 幽灵粒子（反射壁 + 远场）|

所有算例均使用 Laguerre-Gauss 核函数与线性梯度修正矩阵。

---

## 目录结构

```
eulerian-SPH/
├── CMakeLists.txt
├── src/                  # SPHinXsys 核心库
│   ├── shared/
│   ├── for_2D_build/
│   └── 3rd_party/        # 内嵌 tinyxml2
└── cases/
    ├── 2d_eulerian_taylor_green_LG/
    ├── 2d_eulerian_flow_around_cylinder_LG/
    └── 2d_eulerian_supersonic_flow_around_cylinder/
```
