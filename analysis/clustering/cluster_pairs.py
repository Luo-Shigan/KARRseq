#!/usr/bin/env python3

import matplotlib.pyplot as plt
import sys
import gzip
import argparse

import numpy as np
from sklearn.manifold import TSNE
from sklearn.cluster import SpectralClustering
from sklearn.cluster import DBSCAN
from sklearn.neighbors import kneighbors_graph

import networkx as nx
from networkx.algorithms.community import greedy_modularity_communities
"""
    1. 输入处理:

    - 读取gzip压缩的配对文件 (f_pairs)，提取同一转录本的配对。

    - 过滤掉跨度过小或距离过近的配对 (p_long, p_span) 。

    2. 聚类分析:

    - 对每个转录本的数据构建图，节点代表配对，边基于重叠阈值。

    - 使用贪心模块度算法检测社区（聚类）。

    - 保留满足最小读数 (p_min_reads) 的聚类。

    3. 结果输出:

    - BED文件: 记录每个聚类的具体区间，用于基因组浏览器可视化。

    - 文本文件: 记录各聚类的中位坐标。

    - 可选可视化: 生成t-SNE降维图 (需指定-t参数)。
"""

def overlap_function(i1, i2):
    """
    功能: 计算两个区间对 (i1和i2) 在左右两侧的重叠长度。

    输入: 两个区间，每个区间包含左右两端的起止坐标。

    逻辑:

    - 定义内部函数calc_overlap计算单侧重叠。

    - 分别计算左右两侧的重叠，若任一侧无重叠则返回(0, 0)。

    输出: 左右两侧的重叠长度元组。
    """
    i1_s1, i1_e1, i1_s2, i1_e2 = i1[1], i1[2], i1[3], i1[4]
    i2_s1, i2_e1, i2_s2, i2_e2 = i2[1], i2[2], i2[3], i2[4]

    def calc_overlap(a_s, a_e, b_s, b_e):
        if (b_s > a_e) or (a_s > b_e):
            return 0
        else:
            o_s = max(a_s, b_s)
            o_e = min(a_e, b_e)
            return abs(o_e-o_s)

    left = calc_overlap(i1_s1, i1_e1, i2_s1, i2_e1)
    right = calc_overlap(i1_s2, i1_e2, i2_s2, i2_e2)

    if (left == 0) or (right == 0):
        return 0, 0
    else:
        #return min(left, right)
        return (left, right)

def load_sizes(f_sizes):
    """
    功能: 从文件加载染色体/转录本大小信息。

    输出: 返回字典，键为名称，值为大小。
    """
    sizes = {}
    with open(f_sizes, "r") as f:
        for line in f:
            row = line.strip("\r\n").split("\t")
            sizes[row[0]] = int(row[1])
    return sizes

def load_eligible_transcript(f_pairs):
    """
    功能: 统计每个转录本的配对数量。

    逻辑: 仅处理同一转录本内的配对 (row[1] == row[3]）。

    输出: 字典记录各转录本的配对数，用于后续过滤。
    """
    eligible = {}
    with gzip.open(f_pairs, "rt") as f:
        for line in f:
            row = line.strip("\r\n").split("\t")
            if row[1] != row[3]:
                continue
            try:
                eligible[row[1]] += 1
            except KeyError:
                eligible[row[1]] = 1
    return eligible

# currently not functional
def community_detection():
    import community
    partition = community.best_partition(G)

    labels = np.array(list(partition.values()))

    for num in list(set(labels)):

        condition = (labels == num)
        c = np.where(condition)[0]

        if len(c) < 3:
            continue
        out = open("clusters_%s.bed" % num, "w")
        for j in list(c):
            out.write(info[j] + "\n")
        out.close()

def spectral_clustering():
    cluster = SpectralClustering(n_clusters=p_clusters,
                                 affinity="nearest_neighbors",
                                 random_state=386).fit(simimat)
    labels = cluster.labels_

def spectral_clustering_eigengap():
    connectivity = kneighbors_graph(data, n_neighbors=10)
    affinity_matrix = 0.5 * (connectivity + connectivity.T)
    cluster = SpectralClustering(n_clusters=p_clusters,
                                 affinity="precomputed", random_state=386,
                                 assign_labels="kmeans").fit(affinity_matrix)

    k = predict_k(affinity_matrix)
    p_clusters = k
    labels = cluster.labels_

def spectral_clustering_precomputed():
    cluster = SpectralClustering(n_clusters=20, affinity="precomputed",
                                 random_state=386)
    cluster.fit(simimat)
    labels = cluster.labels_

def dbscan():
    cluster = DBSCAN(eps=5, min_samples=3)
    labels = cluster.labels_

def hdbscan():
    cluster = hdbscan.HDBSCAN(min_cluster_size=10, metric="precomputed")
    labels = cluster.fit_predict(distmat)

def clustering():
    print(max(labels))

def call_clusters(data, counter, o_bed, o_output,
                  p_overlap, p_min_reads, p_method="greedy"):
    """
    核心函数: 对输入的配对数据进行聚类。

    步骤:

    1. 构建图结构: 基于重叠阈值 (p_overlap) 在节点间添加边。

    2. 社区发现: 使用贪心模块度算法 (greedy_modularity_communities) 进行聚类。

    输出结果:

    - 生成BED格式文件记录聚类区间。

    - 输出聚类的中位坐标到文本文件。

    关键参数:

    - p_overlap: 边连接的重叠阈值。

    - p_min_reads: 聚类的最小读数要求。
    """
    txt = [] # 保存聚类的中位坐标信息
    
    G = nx.Graph() # 无向图，用于存储节点和边
    m = len(data) # 数据中的配对数量

    fullmat = np.zeros((m, m)) # 存储每对配对的最小片段长度的矩阵
    simimat = np.zeros((m, m)) # 存储没对配对之间的相似性（重叠量）

    ### 构建图与计算相似性
    for i in range(len(data)):
        G.add_node(i) # 添加每个配对为节点
    # 通过双重循环，计算没对配对之间的相似性（重叠量），调用overlap_function来计算它们的重叠
    for i in range(len(data)):
        for j in range(i+1):
            overlap = overlap_function(data[i], data[j])
            simimat[i,j] = min(overlap)
            simimat[j,i] = min(overlap)

            fullmat[i,j] = min(data[i][5], data[i][6], data[j][5], data[j][6])
            fullmat[j,i] = min(data[i][5], data[i][6], data[j][5], data[j][6])

            # 如果两队配对的重量大于阈值p_overlap，则在图中添加一条边连接这两队配对
            #if overlap > p_overlap:
            if overlap[0] > p_overlap and overlap[1] > p_overlap:
                G.add_edge(i, j)

    #DEBUG print G.number_of_nodes()
    #DEBUG print G.number_of_edges()

    ### 聚类
    # 如果选择的方法是 "biconnected"，使用 networkx 库的 biconnected_components 来寻找图中的双连通分量，作为聚类结果。
    if p_method == "biconnected":
        clusters = nx.biconnected_components(G)
    # 如果选择的是 "greedy" 方法，调用 greedy_modularity_communities 函数基于贪心模块度算法来聚类
    # greedy_modularity_communities 会尝试最大化图的模块度，找到密集的子图作为聚类。如果发生除以零的错误（通常表示图中没有连通的节点），函数会返回 counter。
    elif p_method == "greedy":

        try:
            clusters = greedy_modularity_communities(G)
        except ZeroDivisionError:
            return counter
    ### 给每个聚类打标签，并进行输出
    # 初始化 labels 数组，用于标记每个节点属于哪个聚类。
    labels = np.zeros(m, dtype=int)
    # 对每个聚类 c，给它编号（num）。
    for num, c in enumerate(clusters):
        # 将聚类 c 转换成列表形式，方便后续处理。
        c = list(c)
        # 如果当前聚类的大小小于最小读取数 p_min_reads，则跳过这个聚类
        if len(c) < p_min_reads:
            continue
        # 给当前聚类中的所有节点（配对）打上编号 num，表示它们属于这个聚类。
        labels[c] = num

        # 提取当前聚类中所有配对的相关信息：
        # c1：转录本 ID
        # s1 和 e1：片段 1 的起始和结束位置
        # s2 和 e2：片段 2 的起始和结束位置
        c1 = [ data[j][0] for j in c ]
        s1 = [ data[j][1] for j in c ]
        e1 = [ data[j][2] for j in c ]
        s2 = [ data[j][3] for j in c ]
        e2 = [ data[j][4] for j in c ]

        # 计算聚类中每个片段的中位坐标（50th percentile）。
        s1_m = np.percentile(s1, 50)
        e1_m = np.percentile(e1, 50)
        s2_m = np.percentile(s2, 50)
        e2_m = np.percentile(e2, 50)

        ### 输出结果
        # ==============
        # write median
        # ==============
        #print c1[0], int((s1_m+e1_m)//2), int((s2_m+e2_m)//2)

        # ==============
        # write bed
        # ==============
        # 对每个聚类生成一个 BED 格式的记录，描述聚类区域。
        # 使用 counter 和聚类的中位坐标生成唯一的 cid（聚类 ID），并写入 o_bed 文件

        chrom = c1[0]
        cid = "%s|%s_%s|%s:%s:%s:%s:%s" % (counter, chrom, num,
                                           chrom, int(s1_m), int(e1_m),
                                           int(s2_m), int(e2_m))
        for j in c:
            start1, start2, length1, length2, strand = (data[j][1], data[j][3],
                                                        data[j][5], data[j][6],
                                                        data[j][8])

            bed_coords = [chrom, start1, start2+length2, cid, 1000,
                          strand, start1, start2+length2, 0, 2,
                          "%s,%s" % (length1, length2),
                          "%s,%s" % (0, start2-start1)]
            
            if o_bed is not None:
                o_bed.write("\t".join(map(str, bed_coords)) + "\n")

        # =============
        # write txt
        # =============
        # 生成一个包含中位坐标的文本行，并写入到 o_output 文件中。
        out_coords = [c1[0], int(s1_m), int(e1_m), c1[0], int(s2_m), int(e2_m)]
        txt.append(out_coords)
        if o_output is not None:
            o_output.write("\t".join(map(str, out_coords)) + "\n")
        # 更新计数器并返回
        counter += 1
    # 最终返回更新后的计数器 counter 和保存聚类中位坐标的列表 txt。
    return counter, txt

def visualize(simimat, fullmat, f_tsne=None, ax=None):
    """
    功能: 使用t-SNE降维并可视化聚类结果 (需指定-t参数触发)。
    """
    distmat = (fullmat - simimat)

    # compute embedding
    model = TSNE(random_state=386, metric="precomputed")
    embeddings = model.fit_transform(distmat)

    # visualization
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(embeddings[:,0], embeddings[:,1], s=4, c=labels, cmap="Spectral")

    for i in range(max(labels)+1):
        condition = (labels == i)

        ax.annotate(i, embeddings[condition,:].mean(axis=0),
                    horizontalalignment="center",
                    verticalalignment="center", size=10,
                    weight="bold",
                    color="black")

    if f_tsne is not None:
        plt.savefig(f_tsne, dpi=300)

    return fig, ax

def cluster_chimeric_pairs(p_iid, p_genome,
                           p_long, p_span, p_overlap, p_min_reads,
                           p_min_interactions, p_max_interactions,
                           f_sizes, f_tsne, f_pairs, f_output, f_bed):
    """
    主流程函数: 处理输入文件，按转录本分组并调用聚类。

    逻辑:

    - 过滤不符合长度条件 (p_long, p_span) 的配对。

    - 对每个转录本的数据调用call_clusters。

    过滤条件:

    - p_long/p_span: 排除过近或跨度过小的配对。

    - p_min/max_interactions: 根据配对数筛选转录本。
    """
    results = []
    
    eligible = load_eligible_transcript(f_pairs)

    if f_output is not None:
        o_bed = open(f_bed, "w")
        o_output = open(f_output, "w")
    else:
        o_bed = None
        o_output = None
    
    current_iid = None # 当前处理转录本ID
    counter = 0 # 聚类编号计数器
    data = [] # 当前转录本所有嵌合对数据
    with gzip.open(f_pairs, "rt") as f:

        for line in f:

            row = line.strip("\r\n").split("\t")

            # 发生相互作用的两个端不是同一个转录本跳过
            if row[1] != row[3]:
                continue
            # all表示不过滤，如果不是特定的目标转录本跳过
            if row[1] != "all" and row[1] != p_iid:
                continue
            
            if current_iid == None:
                current_iid = row[1]

            # 如果当前转录本的配对数在合格返回，调用 call_cluster对上一组转录本的数据聚类
            elif current_iid != row[1]:
                if (eligible[current_iid] >= p_min_interactions
                    and eligible[current_iid] <= p_max_interactions):
                    counter, txt = call_clusters(data, counter, o_bed, o_output,
                                                 p_overlap, p_min_reads)
                    results.append(txt)
                # 清空data,开始收集新的转录本
                data = []
                current_iid = row[1]


            # append these information
            # 记录符合条件的配对数据
            # c1 、c2是转录本ID；s1、s2是起始位置；l1、l2是片段长度；iid是配对ID；strand1、strand2是链方向
            c1, s1, l1, c2, s2, l2 = (row[1], int(row[2]), int(row[9]),
                                      row[3], int(row[4]), int(row[10]))
            iid, strand1, strand2 = (row[0], row[5], row[6])

            # 如果两个片段太近，跳过
            if abs(s1-(s2+l2)) < p_long:
                continue
            
            # 如果跨度太小，跳过
            if abs(s2-(s1+l1)) < p_span:
                continue
            
            # 如果两个片段在同一条链上，统一处理链方向
            if strand1 == strand2:
                if strand1 == "-":
                    strand = "+"
                else:
                    strand = "-"
            # 将整理好的配对信息加入data，为后续的按转录本聚类做准备
            data.append((c1, s1, s1+l1, s2, s2+l2, l1, l2, iid, strand))

    # ===========
    # last call
    # ===========
    # 检查最后一个 current_iid 是否是合格的转录本（配对数在[min, max]之间）。

    # 如果合格，调用 call_clusters 进行聚类，把结果追加到 results。

    if (eligible[current_iid] >= p_min_interactions
        and eligible[current_iid] <= p_max_interactions):
        counter, txt = call_clusters(data, counter, o_bed, o_output,
                                     p_overlap, p_min_reads)
        results.append(txt)
        
    if f_output is not None:
        o_bed.close()
        o_output.close()

    return results

def main(p_iid, p_genome,
         p_long, p_span, p_overlap, p_min_reads,
         p_min_interactions, p_max_interactions,
         f_sizes, f_tsne, f_pairs, f_output, f_bed):

    cluster_chimeric_pairs(p_iid, p_genome,
                           p_long, p_span, p_overlap, p_min_reads,
                           p_min_interactions, p_max_interactions,
                           f_sizes, f_tsne, f_pairs, f_output, f_bed)

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-c", default="all", help="all|id")
    parser.add_argument("-b", default=None, help="bed")
    parser.add_argument("-o", default=None, help="clusters")
    parser.add_argument("-t", default=None, help="tsne")
    parser.add_argument("-s", default=None, help="sizes file")
    parser.add_argument("-g", default=None, help="genome")
    parser.add_argument("-l", default=200, type=int, help="long range threshold")
    parser.add_argument("-p", default=10, type=int, help="span threshold")
    parser.add_argument("-ol", default=10, type=int, help="overlapping threshold")
    parser.add_argument("-m", default=3, type=int, help="minimum no. of reads per cluster")
    parser.add_argument("-i", default=200, type=int, help="minimum no. of interactions per transcript")
    parser.add_argument("-j", default=20000, type=int, help="max no. of interactions per transcript")
    parser.add_argument("-pairs", default=None, help="paired txt file")
    args = parser.parse_args()

    p_iid = args.c
    p_genome = args.g
    p_long = args.l
    p_span = args.p
    p_overlap = args.ol
    p_min_reads = args.m
    p_min_interactions = args.i
    p_max_interactions = args.j

    f_sizes = args.s
    f_tsne = args.t
    f_pairs = args.pairs
    f_output = args.o
    f_bed = args.b

    main(p_iid=p_iid, p_genome=p_genome,
         p_long=p_long, p_span=p_span, p_overlap=p_overlap,
         p_min_reads=p_min_reads,
         p_min_interactions=p_min_interactions,
         p_max_interactions=p_max_interactions,
         f_sizes=f_sizes, f_tsne=f_tsne, f_pairs=f_pairs, f_output=f_output, f_bed=f_bed)
    
