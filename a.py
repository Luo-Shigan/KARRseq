import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import numpy as np

# 定义主曲线的端点和控制点
# 外贝塞尔曲线（二次）
P0 = np.array([0, 0])    # 起点
C0 = np.array([1.5, -0.5])    # 控制点
P1 = np.array([3, 0])    # 终点

# 内贝塞尔曲线（二次）
P2 = np.array([0, -1])   # 起点
C1 = np.array([1.5, -0.5])   # 控制点
P3 = np.array([3, -1])   # 终点

# 计算闭合曲线的控制点（模拟圆弧连接）
def calculate_arc_controls(start, end, curvature=0.5):
    """计算近似圆弧的二次贝塞尔控制点"""
    mid = (start + end) / 2
    direction = end - start
    normal = np.array([-direction[1], direction[0]]) * curvature
    return mid + normal

# 闭合贝塞尔曲线的控制点（调整curvature参数控制曲率）
C2 = calculate_arc_controls(P1, P3, curvature=0.6)  # 上闭合曲线控制点
C3 = calculate_arc_controls(P2, P0, curvature=0.6) # 下闭合曲线控制点

# 构建路径
path = Path(
    vertices=[
        P0, C0, P1,        # 外曲线（二次贝塞尔）
        C2, P3,            # 上闭合曲线（二次贝塞尔）
        C1, P2,            # 内曲线（二次贝塞尔）
        C3, P0,            # 下闭合曲线（二次贝塞尔）
        P0                 # 闭合路径
    ],
    codes=[
        Path.MOVETO,
        Path.CURVE3,
        Path.CURVE3,
        Path.CURVE3,
        Path.CURVE3,
        Path.CURVE3,
        Path.CURVE3,
        Path.CURVE3,
        Path.CURVE3,
        Path.CLOSEPOLY
    ]
)

# 绘制图形
fig, ax = plt.subplots(figsize=(8, 6))
patch = patches.PathPatch(path, facecolor='lightblue', edgecolor='darkblue', lw=2)
ax.add_patch(patch)

# 绘制控制点和辅助线
for point, label in zip([P0, P1, P2, P3], ['P0', 'P1', 'P2', 'P3']):
    ax.plot(*point, 'ro')
    ax.text(*point, label, fontsize=12)
for cp, label in zip([C0, C1, C2, C3], ['C0', 'C1', 'C2', 'C3']):
    ax.plot(*cp, 'go')
    ax.text(*cp, label, fontsize=12)

ax.set_xlim(-2, 5)
ax.set_ylim(-4, 3)
ax.set_aspect('equal')
plt.grid(True)
plt.title('Geometric Wide Chord with Bezier Curves')
plt.show()