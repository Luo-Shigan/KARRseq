import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from scipy.sparse import random
def calculate_arc_controls(start, end, curvature=0.5):
    """计算近似圆弧的二次贝塞尔控制点"""
    mid = (start + end) / 2
    direction = end - start
    normal = np.array([-direction[1], direction[0]]) * curvature
    return mid + normal

def bezier_control_points(center, p0, p3):
    """
    计算拟合圆弧的三次贝塞尔曲线控制点
    参数:
        center: (x, y) 圆心
        p0: (x, y) 起点
        p3: (x, y) 终点
    返回:
        (p1, p2): 两个控制点
    """
    center = np.array(center)
    p0 = np.array(p0)
    p3 = np.array(p3)

    # 半径
    r = np.linalg.norm(p0 - center)
    
    # 起点和终点对应的角度
    def angle(v):
        return np.arctan2(v[1], v[0])

    theta0 = angle(p0 - center)
    theta1 = angle(p3 - center)
    
    # 确保角度差在 -pi 到 pi 之间
    delta_theta = theta1 - theta0
    while delta_theta <= -np.pi:
        delta_theta += 2 * np.pi
    while delta_theta > np.pi:
        delta_theta -= 2 * np.pi

    # 计算 k 值
    k = (4/3) * np.tan(delta_theta / 4)

    # 起点的切线方向（逆时针旋转90°）
    t0 = np.array([-np.sin(theta0), np.cos(theta0)])
    t1 = np.array([-np.sin(theta1), np.cos(theta1)])

    # 控制点
    p1 = p0 + k * r * t0
    p2 = p3 - k * r * t1

    return p1.tolist(), p2.tolist()
def chordDiagram_fixed(data, ax, colors=None, gap=2 * (np.pi / 180), 
                        linewidth_range=(0.5, 5), alpha_range=(0.3, 0.9)):
    """
    data: (n x n) numpy array，表示互作矩阵
    ax: matplotlib axes
    colors: 节点颜色
    gap: 相邻节点间隔，单位是弧度
    """

    n = data.shape[0]
    # total = np.triu(data).sum() 
    total = np.sum(data)  # 计算总互作强度
    total_strength = data.sum(axis=1)  # 每个节点的总互作强度
    print(f"total: {total}")
    # 每个节点占据的弧长（弧度）
    angles = total_strength / total * (2 * np.pi - n * gap)
    starts = np.zeros(n)
    
    # 计算每个节点的起始位置（弧度）
    for i in range(1, n):
        starts[i] = starts[i-1] + angles[i-1] + gap

    # angles = [np.degrees(angle) for angle in angles]
    # starts = [np.degrees(start) for start in starts]
    # print(f"starts{starts}")
    # print(f"angles{angles}")
    # print(f"{sum(angles)}")
    center = (0, 0)
    radius = 0.85
    nodePos = []

    # 计算节点的位置，并绘制每个节点的弧
    for i in range(n):
        theta1 = np.degrees(starts[i])  # 起始角度
        theta2 = np.degrees(starts[i] + angles[i])  # 结束角度
        # print(f"theta1: {theta1}, theta2: {theta2}")
        x = np.cos(np.radians((theta1 + theta2) / 2)) * radius * 1.2
        y = np.sin(np.radians((theta1 + theta2) / 2)) * radius * 1.2
        angle_deg = (theta1 + theta2) / 2
        nodePos.append((x, y, angle_deg))
        # 绘制每个节点的弧
        ax.add_patch(
            patches.Arc(center, 2 * radius, 2 * radius, 
                        theta1=theta1, theta2=theta2,
                        color=colors[i] if colors else 'grey', lw=10)
        )

    # 画弦
    radius = 0.81
    angles_cp = angles.copy()

    for i in range(n):
        # 当前节点弧的起点
        start_angle_i = starts[i] + angles[i]
        for j in range(i, n):  # 只画上三角（避免重复）
            strength = data[i, j]
            # print(f"data[i,j]:{strength}")
            if strength == 0:
                continue
            
            # 计算这个互作占节点i和节点j各自弧度的比例
            prop_i = strength / total_strength[i]
            prop_j = strength / total_strength[j]

            # 小段弧度
            angle_span_i = angles_cp[i] * prop_i
            angle_span_j = angles_cp[j] * prop_j

            # 起止弧度
            theta_i_start = start_angle_i
            theta_i_end = theta_i_start - angle_span_i


            theta_j_start = starts[j] + angles[j]
            theta_j_end = theta_j_start - angle_span_j
            angles[j] -= angle_span_j 


            # 外贝塞尔曲线
            P0 = np.array([np.cos(theta_i_start), np.sin(theta_i_start)]) * radius
            C0 = np.array([np.cos((theta_i_start+theta_j_end)/2), np.sin((theta_i_start+theta_j_end)/2)]) * radius * 0.2
            P1 = np.array([np.cos(theta_j_end), np.sin(theta_j_end)]) * radius
            # 内贝塞尔曲线
            P2 = np.array([np.cos(theta_i_end), np.sin(theta_i_end)]) * radius
            C1 = np.array([np.cos((theta_i_end+theta_j_start)/2), np.sin((theta_i_end+theta_j_start)/2)]) * radius * 0.2
            P3 = np.array([np.cos(theta_j_start), np.sin(theta_j_start)]) * radius
            
            # 闭合贝塞尔曲线的控制点（调整curvature参数控制曲率）
            C2 = calculate_arc_controls(P1, P3, curvature=-0.235)  # 上闭合曲线控制点
            C3 = calculate_arc_controls(P2, P0, curvature=-0.235) # 下闭合曲线控制点
            C4,C5 = bezier_control_points(center, P1, P3) # 上闭合曲线控制点
            C6,C7 = bezier_control_points(center, P2, P0) # 下闭合曲线控制点
            # 构建路径
            path = Path(
                vertices=[
                    P0, C0, P1,        # 外曲线（二次贝塞尔）
                    C4,C5,P3,            # 上闭合曲线（三次贝塞尔）
                    C1, P2,            # 内曲线（二次贝塞尔）
                    C6,C7, P0,            # 下闭合曲线（三次贝塞尔）
                    P0                 # 闭合路径
                ],
                codes=[
                    Path.MOVETO,Path.CURVE3,Path.CURVE3,
                    Path.CURVE4,Path.CURVE4,Path.CURVE4,
                    Path.CURVE3,Path.CURVE4,
                    Path.CURVE4,Path.CURVE4,Path.CURVE4,
                    Path.CLOSEPOLY
                ]
            )
            cmap = plt.get_cmap('tab20')
            if i == j:
                patch = patches.PathPatch(path, facecolor=colors[i], edgecolor='none', alpha = 0.5)
            else:
                patch = patches.PathPatch(path, facecolor=cmap(colorH), edgecolor='none', alpha = 0.5)
            ax.add_patch(patch)
            start_angle_i -= angle_span_i
            colorH += 1


    return nodePos
if __name__ == "__main__":
    """
    输入为稀疏矩阵
    """
    # 9个节点的互作矩阵
    nodes = ["A", "B", "C", "D", "E", "F", "G", "H", "I"]
    # data = np.random.poisson(5, size=(9, 9))
    # data = np.triu(data, 1)  # 上三角
    # data = data + data.T     # 对称
    # print(data)
    # seed = 40
    # np.random.seed(seed)
    data = random(9, 9, density=0.2, format='csr', random_state=42)  # 稀疏矩阵
    data = data.multiply(np.random.poisson(5, size=(9, 9)))  # 乘以 Poisson 分布的值
    data = np.triu(data.toarray(), 1)  # 上三角矩阵
    data = data + data.T  # 对称
    diag_values = np.random.poisson(2, size=9)
    np.fill_diagonal(data, diag_values)
    print(data)

    colors = ["#c23f76", "#156697", "#f2cc21", "#9994c2", "#34a198",
            "#5b5b6d", "#9e1f63", "#a0522d", "#f99533"]

    # 画图
    fig = plt.figure(figsize=(3, 3))  # 调整画布尺寸
    ax = plt.axes([0, 0, 1, 1])
    nodePos = chordDiagram_fixed(data, ax, colors=colors)
    ax.axis('off')
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    # 标注节点
    for i in range(len(nodes)):
        ax.text(nodePos[i][0], nodePos[i][1], nodes[i], 
                rotation=nodePos[i][2], ha='center', va='center', fontsize=14)

    plt.savefig("chord_diagram.png", dpi=300, bbox_inches='tight')