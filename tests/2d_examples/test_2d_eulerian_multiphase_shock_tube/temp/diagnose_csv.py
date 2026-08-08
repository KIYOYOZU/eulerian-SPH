# 诊断脚本：读 output/ 里最新一帧 particles_*.csv，报告 rho/p/u/alpha 的极值及位置。
# 用法: python temp/diagnose_csv.py [可选: 指定 csv 路径]
import glob
import os
import sys

import numpy as np

case_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
out_dir = os.path.join(case_dir, "output")

if len(sys.argv) > 1:
    csv_path = sys.argv[1]
else:
    files = glob.glob(os.path.join(out_dir, "particles_*.csv"))
    if not files:
        print("no csv found in", out_dir)
        sys.exit(1)

    def step_of(f):
        return int(os.path.basename(f)[len("particles_"):-len(".csv")])

    csv_path = max(files, key=step_of)

time_path = csv_path.replace("particles_", "time_").replace(".csv", ".txt")
t = float(open(time_path).read().strip()) if os.path.exists(time_path) else -1.0

data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
x, y, rho, p, u, alpha = data.T

print(f"file: {os.path.basename(csv_path)}  t={t:.6g}  N={len(x)}")
for name, arr in [("rho", rho), ("p", p), ("u", u), ("alpha", alpha)]:
    imin, imax = int(np.argmin(arr)), int(np.argmax(arr))
    print(f"  {name:5s}: min={arr[imin]:12.5g} @x={x[imin]:.4f} | "
          f"max={arr[imax]:12.5g} @x={x[imax]:.4f}")

# 物理参考范围（两气体 Sod）：rho∈[0.125,1], p∈[0.1,0.425], u∈[0,~0.45]
print("  bounded(rho<=1.1):", bool((rho <= 1.1).all()),
      " | p>0:", bool((p > 0).all()),
      " | u>=-0.1:", bool((u >= -0.1).all()),
      " | alpha in [0,1]:", bool(((alpha >= -1e-12) & (alpha <= 1 + 1e-12)).all()))
