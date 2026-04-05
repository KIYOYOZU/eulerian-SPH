# eulerian-SPH

基于 [SPHinXsys](https://github.com/Xiangyu-Hu/SPHinXsys) 框架的轻量级 2D 欧拉 SPH 求解器，包含弱可压缩与可压缩流动验证算例。

---

## 依赖

| 依赖 | 版本要求 | 说明 |
|------|----------|------|
| CMake | ≥ 3.20 | 构建系统 |
| GCC / Clang | C++20（GCC 10+，Clang 12+）| 编译器 |
| [Eigen3](https://eigen.tuxfamily.org) | ≥ 3.4 | 向量/矩阵运算（仅头文件） |
| [TBB](https://github.com/oneapi-src/oneTBB) | ≥ 2021 | 并行任务调度 |
| [Boost](https://www.boost.org) | ≥ 1.74（仅头文件）| 多边形几何域构造 |

### 系统包管理器安装（有 root 权限）

```bash
# Ubuntu / Debian
sudo apt install cmake g++ libeigen3-dev libtbb-dev libboost-dev

# CentOS / RHEL
sudo dnf install cmake gcc-c++ eigen3-devel tbb-devel boost-devel
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

## 编译

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

cmake -S "$ROOT_DIR" -B "$BUILD_DIR" \
  -Wno-dev \
  -DCMAKE_BUILD_TYPE=Release \
  -DEigen3_DIR=$DEPS/eigen3/share/eigen3/cmake \  # Eigen3 的 cmake 配置文件目录
  -DTBB_DIR=$DEPS/tbb/lib64/cmake/TBB \           # TBB 的 cmake 配置文件目录
  -DBOOST_ROOT=$DEPS/boost_1_82_0                 # Boost 根目录（含 include/）
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

### 粒子松弛（仅圆柱算例）

对于含固体边界的算例，初次运行前可先执行粒子松弛，生成贴体粒子分布以提升精度：

```bash
cd cases/2d_eulerian_flow_around_cylinder_LG

# 第一步：松弛运行，生成贴体粒子坐标，结果写入 reload/
../../build/cases/2d_eulerian_flow_around_cylinder_LG/bin/2d_eulerian_flow_around_cylinder_LG --r=true

# 第二步：正式仿真，从 reload/ 加载松弛后的粒子分布
../../build/cases/2d_eulerian_flow_around_cylinder_LG/bin/2d_eulerian_flow_around_cylinder_LG --reload=true
```

---

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
