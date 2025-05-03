import numpy as np
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
from scipy.sparse import random
import pandas as pd
from typing import List, Tuple
def calculate_arc_controls(start, end, curvature=0.5):
    """计算近似圆弧的二次贝塞尔控制点"""
    mid = (start + end) / 2
    direction = end - start
    normal = np.array([-direction[1], direction[0]]) * curvature
    return mid + normal

def bezier_control_points(center, p0, p3):
    """
    计算拟合圆弧的三次贝塞尔曲线控制点(<180度)
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

def bezierPlot(theta_i_start,theta_i_end,theta_j_start,theta_j_end,radius,center):
    # 外贝塞尔曲线
    P0 = np.array([np.cos(theta_i_start), np.sin(theta_i_start)]) * radius
    C0 = np.array([np.cos((theta_i_start+theta_j_end)/2), np.sin((theta_i_start+theta_j_end)/2)]) * radius * 0.2
    P1 = np.array([np.cos(theta_j_end), np.sin(theta_j_end)]) * radius
    # 内贝塞尔曲线
    P2 = np.array([np.cos(theta_i_end), np.sin(theta_i_end)]) * radius
    C1 = np.array([np.cos((theta_i_end+theta_j_start)/2), np.sin((theta_i_end+theta_j_start)/2)]) * radius * 0.2
    P3 = np.array([np.cos(theta_j_start), np.sin(theta_j_start)]) * radius
    
    # 闭合贝塞尔曲线的控制点（调整curvature参数控制曲率）
    C2 = calculate_arc_controls(P1, P3, curvature=-0.235)  # 上闭合曲线控制点（二阶贝塞尔）
    C3 = calculate_arc_controls(P2, P0, curvature=-0.235) # 下闭合曲线控制点（二阶贝塞尔）

    C4,C5 = bezier_control_points(center, P1, P3) # 上闭合曲线控制点（三阶贝塞尔）
    C6,C7 = bezier_control_points(center, P2, P0) # 下闭合曲线控制点（三阶贝塞尔）
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
    return path

def generate_closed_bezier_path(
    thetas: List[float],
    curve_types: List[int],
    radius:float,
    center:Tuple[float,float]
) -> Path:
    """
    根据任意数量的点和贝塞尔阶数生成闭合路径。

    参数：
        thetas: 弧度列表代表点
        curve_types: 与每个点对对应的整数列表 (2或3)，表示是二次还是三次贝塞尔
        calculate_arc_controls(p0, p1): 计算二次贝塞尔控制点
        bezier_control_points(p0, p1): 计算三次贝塞尔控制点对（默认中心为中点）
    
    返回：
        Path 对象，可用于 PathPatch 绘图
    """
    points = [np.array([np.cos(theta), np.sin(theta)]) * radius for theta in thetas]
    assert len(points) >= 2, "在圆上只需要2个点就能构成闭合路径"
    assert len(points) == len(curve_types), "闭合路径点数与曲线类型数必须一致"

    vertices = []
    codes = []

    n = len(points)
    for i in range(n):
        p0 = np.array(points[i])
        p1 = np.array(points[(i + 1) % n]) # 闭合路径，当i=n-1,最后一个点为最开始的点
        curve_type = curve_types[i]

        if i == 0:
            vertices.append(p0)
            codes.append(Path.MOVETO)

        if curve_type == 2:
            ctrl = np.array([np.cos((thetas[i]+thetas[(i + 1) % n])/2), np.sin((thetas[i]+thetas[(i + 1) % n])/2)]) * radius * 0.2
            vertices.extend([ctrl, p1])
            codes.extend([Path.CURVE3, Path.CURVE3])
        elif curve_type == 3:
            c1, c2 = bezier_control_points(center, p0, p1)
            vertices.extend([c1, c2, p1])
            codes.extend([Path.CURVE4, Path.CURVE4, Path.CURVE4])
        else:
            raise ValueError("curve_types 只能为 2 或 3")

    # 闭合
    vertices.append(points[0])
    codes.append(Path.CLOSEPOLY)

    return Path(vertices, codes)





    

def chordDiagram_fixed(data, ax, colors=None, gap=2 * (np.pi / 180)):
    """
    data: (n x n) numpy array, 表示互作矩阵
    ax: matplotlib axes
    colors: 节点颜色
    gap: 相邻节点间隔，单位是弧度
    """

    n = data.shape[0]
    # total = np.triu(data).sum() 
    total = np.sum(data)  # 计算总互作强度
    total_strength = data.sum(axis=1)  # 每个节点的总互作强度
    print(f"total: {total}")
    # 每个节点总互作强度占据的弧度
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
    colorH = 0
    cmap = plt.get_cmap('tab20')
    for i in range(n):
        # 当前节点弧的起点
        start_angle_i = starts[i] + angles[i]

        for j in range(i, n):  # 只画上三角（避免重复）
            strength = data[i, j]
            # print(f"data[i,j]:{strength}")
            if strength == 0:
                continue

            if i == j: # 自身连接主点只有两个
                prop = strength / total_strength[i]
                angle_span = angles_cp[i] * prop
                angles[i] -= angle_span

                theta_start = start_angle_i
                theta_end = theta_start - angle_span
                if angle_span <= np.pi: 
                    thetas = [theta_start,theta_end]
                    curve_types = [2,3]
                    path = generate_closed_bezier_path(thetas,curve_types,radius,center)

                else: # 当角度大于180度，三次贝塞尔曲线不能很好地拟合圆弧
                    theta_sep = theta_start - (angle_span - np.pi)
                    thetas = [theta_start,theta_end,theta_sep]
                    curve_types = [2,3,3]
                    path = generate_closed_bezier_path(thetas,curve_types,radius,center)                   

                
                patch = patches.PathPatch(path, facecolor=colors[i], edgecolor='none', alpha = 0.5)

                ax.add_patch(patch)

                start_angle_i -= angle_span #更新角度

            else:   # 不同节点连接主点有四个
                # 计算这个互作占节点i和节点j各自弧度的比例
                prop_i = strength / total_strength[i]
                prop_j = strength / total_strength[j]

                # 小段弧度
                angle_span_i = angles_cp[i] * prop_i
                angle_span_j = angles_cp[j] * prop_j
                
                theta_i_start = start_angle_i
                theta_i_end = theta_i_start - angle_span_i

                theta_j_start = starts[j] + angles[j]
                theta_j_end = theta_j_start - angle_span_j
                angles[j] -= angle_span_j 

                if angle_span_i <= np.pi and angle_span_j <= np.pi:
                    thetas = [theta_i_start, theta_j_end, theta_j_start, theta_i_end]
                    curve_types = [2,3,2,3]
                    path = generate_closed_bezier_path(thetas,curve_types,radius,center)
                elif angle_span_i > np.pi:
                    theta_i_sep = theta_i_start - (angle_span_i - np.pi)
                    thetas = [theta_i_start, theta_j_end, theta_j_start, theta_i_end,theta_i_sep]
                    curve_types = [2,3,2,3,3]
                    path = generate_closed_bezier_path(thetas,curve_types,radius,center)
                elif angle_span_j > np.pi:
                    theta_j_sep = theta_j_start - (angle_span_j - np.pi)
                    thetas = [theta_i_start, theta_j_end, theta_j_sep, theta_j_start, theta_i_end]
                    curve_types = [2,3,3,2,3]
                    path = generate_closed_bezier_path(thetas,curve_types,radius,center)
                else: #当 i !=j 的情况下，这种条件不可能成立
                    theta_i_sep = theta_i_start - (angle_span_i - np.pi)
                    theta_j_sep = theta_j_start - (angle_span_j - np.pi)
                    thetas = [theta_i_start, theta_j_end,theta_j_sep, theta_j_start, theta_i_end,theta_i_sep]
                    curve_types = [2,3,3,2,3,3]
                    path = generate_closed_bezier_path(thetas,curve_types,radius,center)
                
                patch = patches.PathPatch(path, facecolor=cmap(colorH), edgecolor='none', alpha = 0.5)
                
                ax.add_patch(patch)
                start_angle_i -= angle_span_i #更新角度
            colorH += 1 #更新颜色


    return nodePos
def parse_results1(data, nodes, f_results, dtype):
    """
    get Interaction matrix from annotated interaction file
    """
    if dtype == "inter":
        with open(f_results, "r") as f:
            for line in f:
                row = line.strip("\r\n").split("\t")
                t1, t2, count = row[0], row[1], int(row[2])
                
                i = nodes.index("inter-"+t1)
                j = nodes.index("inter-"+t2)
                
                if i == j:
                    data[i,j] += count
                else:
                    data[i,j] += count
                    data[j,i] += count
                
    elif dtype == "intra":
        with open(f_results, "r") as f:
            for line in f:
                row = line.strip("\r\n").split("\t")
                t1, t2, count = row[0], row[1], int(row[2])

                if t1 == "intergenic" or t2 == "intergenic":
                    continue
                
                i = nodes.index("intra-"+t1)
                j = nodes.index("intra-"+t2)

                if i == j:
                    data[i,j] += count
                else:
                    data[i,j] += count
                    data[j,i] += count
    
    return data
def parse_results(pair:str,data,nodes):
    """
    get statistic count from annotated interaction file
    """
    df = pd.read_table(pair,sep="\t")
    sta_df = df[['Target_Label','Control_Label']].apply(lambda row: tuple(sorted(row)), axis=1).value_counts()
    results = []
    for index,count in sta_df.items():
        l1,l2 = index
        results.append((l1,l2,count))
        i = nodes.index(l1)
        j = nodes.index(l2)
        if i == j:
            data[i,j] += count
        else:
            data[i,j] += count
            data[j,i] += count
    return data

        

def main():
    nodes = ["Exon (cis)", "3' UTR (cis)","5' UTR (cis)",
            "Exon (trans)","3' UTR (trans)", "5' UTR (trans)"]


    data = np.array(np.zeros((len(nodes), len(nodes))))

    data = parse_results("b.csv",data, nodes)
    colors = ["#c23f76", "#156697", "#f2cc21", 
            "#5b5b6d", "#9e1f63", "#a0522d"]
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
                rotation=nodePos[i][2], ha='center', va='center', fontsize=3)

    plt.savefig("chord_diagram.png", dpi=300, bbox_inches='tight')

def simulated():
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
                rotation=nodePos[i][2], ha='center', va='center', fontsize=2)

    plt.savefig("chord_diagram.png", dpi=300, bbox_inches='tight')
if __name__ == "__main__":
    # simulated()
    main()

