# 打印指定 CSV 在界面窗口 [x0, x1] 内的粒子剖面（按 x 排序）。
# 用法: python temp/profile_interface.py <csv> [x0=0.66] [x1=0.74]
import sys

import numpy as np

csv_path = sys.argv[1]
x0 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.66
x1 = float(sys.argv[3]) if len(sys.argv) > 3 else 0.74

time_path = csv_path.replace("particles_", "time_").replace(".csv", ".txt")
try:
    t = float(open(time_path).read().strip())
except Exception:
    t = -1.0

data = np.loadtxt(csv_path, delimiter=",", skiprows=1)
x, y, rho, p, u, alpha = data.T
idx = np.where((x >= x0) & (x <= x1))[0]
idx = idx[np.argsort(x[idx])]
print(f"t={t:.6g}  window [{x0},{x1}]  N={len(idx)}")
for j in idx:
    print(f"  x={x[j]:.4f} rho={rho[j]:12.5g} p={p[j]:12.5g} u={u[j]:12.5g} alpha={alpha[j]:8.5g}")
