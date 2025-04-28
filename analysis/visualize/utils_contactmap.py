#!/usr/bin/env python3

import matplotlib

import sys
import gzip
import argparse
import pypairix
import numpy as np
import scipy.stats as spstats

import pandas as pd
import pyBigWig
from pyfaidx import Fasta
import matplotlib.pyplot as plt
from .utils_suboptimal_structure import get_matrix, get_bppm



def parse_structure(f_structure, p_res, p_display):
    """
    ​功能: 解析RNA二级结构文件, 生成接触矩阵\n
    ​参数:
        f_structure: 结构文件路径
        p_res: 分辨率(每个bin包含的碱基数)
        p_display: 显示模式("full", "lower"或"upper")\n
    ​返回值: 表示RNA二级结构接触的矩阵 
    """
    struct = ""
    with open(f_structure, "r") as f:
        for line in f:
            struct += line.strip("\r\n")
    n_bins = len(struct)//p_res
    data = np.zeros((n_bins+1, n_bins+1))
    #print(struct.count("("), struct.count(")"), struct.count("."))
    ds = []
    for j, fold in enumerate(struct):
        if fold == "(":
            ds.append(j)
        elif fold == ")":
            i = ds.pop()
            bin_i = i//p_res
            bin_j = j//p_res
            if p_display in ("full", "lower"):
                data[bin_i][bin_j] += 1
            if p_display in ("full", "upper"):
                data[bin_j][bin_i] += 1
    return data

def get_size(f_size):
    """
    ​功能: 从染色体大小文件中读取染色体名称和大小\n
    ​参数: 
        f_size: 染色体大小文件路径\n
    ​返回值: 字典{染色体名: 大小}
    """
    sizes = {}
    with open(f_size, "r") as f:
        for line in f:
            row = line.strip("\r\n").split("\t")
            sizes[row[0]] = int(row[1])
    return sizes

def get_annotation(f_annotation):
    """
    ​功能: 从注释文件中获取基因注释信息\n
    ​参数: 
        f_annotation: 注释文件路径\n
    ​返回值: 字典{基因ID: (起始位置, 结束位置)}
    """
    annotations = {}
    with open(f_annotation, "r") as f:
        for line in f:
            row = line.strip("\r\n").split("\t")
            if int(row[2]) == int(row[3]):  # entire CDS is UTR
                annotations[row[0]] = (int(row[1]), int(row[4]))
            else:
                annotations[row[0]] = (int(row[2]), int(row[3]))
    return annotations

def get_splice(f_bed12):
    """
    功能: 从BED12格式文件中获取剪接位点信息\n
    ​参数: 
        f_bed12: BED12文件路径\n
    ​返回值: 字典{转录本ID: 剪接位点数组}
    """
    splices = {}
    with open(f_bed12, "r") as f:
        for line in f:
            row = line.strip("\r\n").split("\t")
            iid = row[3]
            strand = row[5]
            exon_blocks = np.array(list(map(int, row[10].split(",")[:-1])))
            if strand == "+":
                splice_sites = exon_blocks.cumsum()
            elif strand == "-":
                splice_sites = exon_blocks[::-1].cumsum()

            splices[iid] = splice_sites[:-1]
    return splices

def get_longest_by_gene(f_index, gene):
    """
    功能: 获取基因的最长转录本信息\n
    ​参数:
        f_index: 索引文件路径
        gene: 基因名称\n
    ​返回值: (转录本ID, 起始位置, 结束位置)
    """
    index = {}
    with open(f_index, "r") as f:
        for line in f:
            row = line.strip("\r\n").split("\t")
            try:
                index[row[1]].append(row)
            except KeyError:
                index[row[1]] = [row]

    try:
        query = index[gene][0]  # longest
        rid, start, end = (query[0], 0, int(query[2])-1)
    except KeyError:
        rid, start, end = (None, None, None)

    return rid, start, end

def get_by_rid(f_index, rid):
    """
    功能: 根据转录本ID获取其坐标信息\n
    ​参数:
        f_index: 索引文件路径
        rid: 转录本ID\n
    ​返回值: (起始位置, 结束位置)
    """
    index = {}
    with open(f_index, "r") as f:
        for line in f:
            row = line.strip("\r\n").split("\t")
            index[row[0]] = row

    try:
        query = index[rid]
        start, end = (0, int(query[2])-1)
    except KeyError:
        start, end = (None, None)

    return start, end

def load_loops(chrom, p_start, p_end, f_loop_path,
               p_resolution, p_point=False):
    """
    ​功能: 加载环状结构(loops)数据\n
    ​参数:
        chrom: 染色体名称
        p_start: 起始位置
        p_end: 结束位置
        f_loop_path: loops文件路径
        p_resolution: 分辨率
        p_point: 是否作为点处理\n
    ​返回值: loops列表
    """
    loop = []
    with open(f_loop_path, "r") as f:
        for line in f:
            if line.startswith("track"):
                continue
            row = line.strip("\r\n").split("\t")
            #lchrom, lstart1, lend1, lstart2, lend2 = (row[0],
            #                                          int(row[1]), int(row[1])+10,
            #                                          int(row[2]), int(row[2])+10)

            if len(row) > 6:
                (lchrom1, lstart1, lend1,
                 lchrom2, lstart2, lend2, category) = (row[0], int(row[1]), int(row[2]),
                                                       row[3], int(row[4]), int(row[5]),
                                                       row[6])
            else:
                (lchrom1, lstart1, lend1,
                 lchrom2, lstart2, lend2) = (row[0], int(row[1]), int(row[2]),
                                             row[3], int(row[4]), int(row[5]))

            # TODO intra: check whether they are in the same chromosome
            # TODO inter: no checks
            #if lchrom != chrom:
            #    continue

            #if lstart1 >= p_start and lend2 <= p_end:
            if len(row) > 6:
                if p_point:
                    p_pres = p_resolution * 2
                    loop.append((lstart1-p_pres, lstart1+p_pres,
                                 lstart2-p_pres, lstart2+p_pres, category))
                else:
                    loop.append((lstart1, lend1, lstart2, lend2, category))
            else:
                if p_point:
                    p_pres = p_resolution * 2
                    loop.append((lstart1-p_pres, lstart1+p_pres,
                                 lstart2-p_pres, lstart2+p_pres))
                else:
                    loop.append((lstart1, lend1, lstart2, lend2))
    return loop

def tabulate(data, f_pair, loci, p_start1, p_start2,
             p_end1, p_end2,
             p_chrom1, p_chrom2, p_resolution,
             p_downsample=1.0):
    """
    ## ​功能: 从pairix格式文件中统计接触频率\n
    ## 参数:
    - data: 预初始化的零矩阵，用于存储结果
    - f_pair: pairix 格式数据文件路径
    - loci: 查询的基因组区域列表
    - p_start1/2, p_end1/2: 区域起止位置
    - p_chrom1/2: 染色体名称
    - p_resolution: 分析分辨率 (bin 大小) 
    - p_downsample: 下采样比例 (默认1.0, 即不降采样）
    ## ​返回值: 
    填充后的接触矩阵
    """
 # 设置随机种子以确保可重复的下采样
    np.random.seed(386)

    # 遍历每个查询区域（loci可能包含多个区域）
    for no, locus in enumerate(loci):

        # 打开pairix文件并获取指定区域的迭代器
        pairs = pypairix.open(f_pair)
        iterator = pairs.querys2D(locus)  # 获取该区域的所有接触对

        # 遍历每个接触对
        for row in iterator:

            # 下采样：随机跳过部分数据（默认100%保留）
            if np.random.random() > p_downsample:
                continue

            # 解析接触对的两个坐标信息
            chrom1, start1, chrom2, start2 = (row[1], int(row[2]),
                                              row[3], int(row[4]))

            # 分支1：处理两个独立定义区域的接触（例如不同染色体）
            if p_start1 and p_end1 and p_start2 and p_end2:
                # 计算bin索引（基于各自区域的起始位置）
                binth1 = (start1 - p_start1) // p_resolution
                binth2 = (start2 - p_start2) // p_resolution
                # 非对称填充（仅单方向）
                data[binth1][binth2] += 1

            # 分支2：处理单一定义区域内的接触（例如同一染色体内）
            elif p_start1 and p_end1:
                # 计算bin索引（基于同一区域的起始位置）
                binth1 = (start1 - p_start1) // p_resolution
                binth2 = (start2 - p_start1) // p_resolution
                # 对称填充（交互是双向的）
                data[binth1][binth2] += 1
                data[binth2][binth1] += 1

            # 分支3：处理跨染色体且需要翻转坐标的情况
            elif p_chrom1 and p_chrom2:  # flip
                # 直接按分辨率计算bin（假设全局坐标）
                binth1 = start1 // p_resolution
                binth2 = start2 // p_resolution
                # 根据遍历的locus序号决定填充方向
                if no == 0:
                    data[binth1][binth2] += 1
                elif no == 1:
                    data[binth2][binth1] += 1

            # 分支4：处理单染色体内的全局统计
            elif p_chrom1:
                # 直接按分辨率计算bin
                binth1 = start1 // p_resolution
                binth2 = start2 // p_resolution
                # 对称填充
                data[binth1][binth2] += 1
                data[binth2][binth1] += 1

    return data

def normalize_by_diag(X):
    """
    ​功能: 按对角线距离对接触矩阵进行归一化\n
    ​参数: 
        X: 原始接触矩阵\n
    ​返回值: 归一化后的矩阵
    """
    pseudocount = 1e-10
    N = np.zeros(X.shape)
    for i in range(X.shape[0]):
        pos_indices = np.where(np.eye(X.shape[0], k=i)==1)
        neg_indices = np.where(np.eye(X.shape[0], k=-i)==1)
        N[pos_indices] = np.log2((X[pos_indices] + pseudocount) / (np.diagonal(X, i).mean() + pseudocount))
        N[neg_indices] = np.log2((X[neg_indices] + pseudocount) / (np.diagonal(X, -i).mean() + pseudocount))
    return N

def normalize_diff(X0, X1):
    """
    ​功能: 计算两个接触矩阵的差异并归一化\n
    ​参数:
        X0: 第一个矩阵
        X1: 第二个矩阵\n
    ​返回值: 归一化后的差异矩阵
    """
    pseudocount = 1e-10
    Xm = (X0 + X1)/2
    Xd = X0 - X1
    print(Xm.shape)
    print(Xd.shape)
    N = np.zeros(Xm.shape)
    for i in range(Xm.shape[0]):
        pos_indices = np.where(np.eye(Xm.shape[0], k=i)==1)
        neg_indices = np.where(np.eye(Xm.shape[0], k=-i)==1)
        N[pos_indices] = (Xd[pos_indices] + pseudocount) / (np.diagonal(Xm, i).mean() + pseudocount)
        N[neg_indices] = (Xd[neg_indices] + pseudocount) / (np.diagonal(Xm, -i).mean() + pseudocount)
    return N

def get_ref(p_genome):
    """
    ​功能: 根据基因组版本获取染色体大小文件和fasta文件路径\n
    ​参数: 
        p_genome: 基因组版本(如"mm10", "hg19")\n
    ​返回值: (染色体大小文件路径, fasta文件路径)
    """
    proj_dir = ""
    dict_sizes = {"mm10": proj_dir + "/ChIP_seq_2/Data/index/Mus_musculus/NCBI/GRCm38.p6/mm10refseq_wo_version.fa.fai",
                  "hg19": proj_dir + "/ChIP_seq_2/Data/index/Homo_sapiens/NCBI/GRCh37.p13/hrefseq_wo_version.fa.fai",
                  "GRCm39": proj_dir + "/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa.fai",
                  "GRCh38": proj_dir + "/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/GRCh38.primary_assembly.genome.fa.fai",
                  "GRCm39T": proj_dir + "/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/transcripts.fai"}
    dict_fasta = {"mm10": proj_dir + "/ChIP_seq_2/Data/index/Mus_musculus/NCBI/GRCm38.p6/mm10refseq_wo_version.fa",
                  "hg19": proj_dir + "/ChIP_seq_2/Data/index/Homo_sapiens/NCBI/GRCh37.p13/hrefseq_wo_version.fa",
                  "GRCm39": proj_dir + "/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa",
                  "GRCh38": proj_dir + "/ChIP_seq_2/Data/index/Homo_sapiens/GENCODE/GRCh38/GRCh38.primary_assembly.genome.fa",
                  "GRCm39T": proj_dir + "/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/transcripts.fa" }
    try:
        f_size = dict_sizes[p_genome]
        f_fasta = dict_fasta[p_genome]
    except KeyError:
        f_size = None
        f_fasta = None
    return f_size, f_fasta

def get_loci(p_chrom1=None, p_chrom2=None,
             p_start1=None, p_start2=None,
             p_end1=None, p_end2=None, p_gene1=None, p_gene2=None,
             p_resolution=100, p_genome=None, f_size=None, f_index=None, **kwargs):
    """
    ​功能: 根据输入参数生成查询区域和bin信息\n
    ​参数: 
        - p_chrom1 和 p_chrom2: 染色体名称
        - p_start1 和 p_start2: 起始位置
        - p_end1 和 p_end2: 结束位置
        - p_gene1 和 p_gene2: 基因名称
        - p_resolution: 分辨率, 默认100
        - p_genome: 基因组标识
        - f_size 和 f_index: 文件相关参数
        - **kwargs: 其他任意关键字参数 \n
    ​返回值: 包含区域信息和bin信息的元组
    """
    if p_genome:
        f_size, f_fasta = get_ref(p_genome)

    sizes = get_size(f_size)

    bins1 = None
    bins2 = None
    # 仅提供染色体名称时，生成对应区间；为分析区域内相互作用初始化默认值
    if p_chrom1 and (p_start1 is None) and (p_end1 is None):
        p_start1, p_end1 = 0, sizes[p_chrom1]-1
    if p_chrom2 and (p_start2 is None) and (p_end2 is None):
        p_start2, p_end2 = 0, sizes[p_chrom2]-1
    # 分析两个特定区间的相互作用
    if (p_start1 is not None) and (p_end1 is not None) and (p_start2 is not None) and (p_end2 is not None):
        bins1 = (p_end1 - p_start1)//p_resolution
        bins2 = (p_end2 - p_start2)//p_resolution
        locus = "%s:%s-%s|%s:%s-%s" % (p_chrom1, p_start1, p_end1, p_chrom2, p_start2, p_end2)
        loci = [locus]
        m_bins = bins1+1
        n_bins = bins2+1
        #DEBUG print("A")
    # 分析区域内相互作用
    elif (p_start1 is not None) and (p_end1 is not None):
        p_start1, p_end1 = (int(p_start1), int(p_end1))
        bins1 = (p_end1 - p_start1)//p_resolution
        bins2 = None
        locus = "%s:%s-%s|%s:%s-%s" % (p_chrom1, p_start1, p_end1, p_chrom1, p_start1, p_end1)
        loci = [locus]
        m_bins = bins1+1
        n_bins = bins1+1
        #DEBUG print("B")
    # 分析两个染色体之间的相互作用
    elif (p_chrom1 is not None) and (p_chrom2 is not None):
        bins1 = sizes[p_chrom1]//p_resolution
        bins2 = sizes[p_chrom2]//p_resolution
        locus1 = "%s|%s" % (p_chrom1, p_chrom2)
        locus2 = "%s|%s" % (p_chrom2, p_chrom1)
        loci = [locus1, locus2]
        m_bins = bins1+1
        n_bins = bins2+1
        #DEBUG print("C")
    print(f"Debug: After get_ref(), f_size = {f_size}") 
    return (loci, p_chrom1, p_start1, p_end1,
            p_chrom2, p_start2, p_end2,
            m_bins, n_bins, bins1, bins2)

def get_contactmap(f_pairs, p_chrom1=None, p_chrom2=None,
                   p_start1=None, p_start2=None,
                   p_end1=None, p_end2=None, p_gene1=None, p_gene2=None,
                   p_resolution=100, p_multi=False, p_diff=False,
                   p_genome=None, f_size=None, f_index=None, p_downsample=1.0,
                   p_point=False, **kwargs):
    """
    get_contactmap 是一个用于从交互数据生成接触矩阵 (contact matrix) 的核心函数, 主要用于构建转录本交互的热图数据。

    ## 函数功能
    1. 获取基因组区域信息：通过调用 get_loci 确定分析区域
    2. ​构建接触矩阵：根据输入数据创建 N*N 的交互矩阵
    3. 支持多种数据处理模式：
        - 单样本累积模式
        - 多样本独立模式
        - 差异比较模式
    4. ​数据下采样：支持对原始数据进行降采样处理
    ## 参数解析
    ### 主要输入参数
    - **f_pairs**: Hi-C 交互数据文件路径（单个文件或多个文件）
    ### 基因组坐标参数:
    - p_chrom1, p_chrom2: 染色体名称
    - p_start1, p_end1: 区域起始和结束位置
    - p_gene1, p_gene2: 基因名称（可替代坐标）
    ### 分析参数:
    - p_resolution: 分析分辨率 (默认100bp) 
    - p_multi: 是否处理多个独立样本
    - p_diff: 是否准备差异分析数据
    - p_downsample: 下采样比例 (0.0-1.0) 
    """
    (loci, p_chrom1, p_start1, p_end1,
     p_chrom2, p_start2, p_end2,
     m_bins, n_bins, bins1, bins2) = get_loci(p_chrom1=p_chrom1, p_chrom2=p_chrom2,
                                              p_start1=p_start1, p_start2=p_start2,
                                              p_end1=p_end1, p_end2=p_end2, p_gene1=p_gene1,
                                              p_gene2=p_gene2,
                                              p_resolution=p_resolution, p_genome=p_genome,
                                              f_size=f_size, f_index=f_index)

    raw = []
    if p_multi:

        for f_pair in f_pairs:
            for f_split in f_pair.split(","):
                datum = np.zeros((m_bins, n_bins))
                datum = tabulate(datum, f_split, loci,
                                 p_start1, p_start2,
                                 p_end1, p_end2,
                                 p_chrom1, p_chrom2, p_resolution,
                                 p_downsample)
            raw.append(datum)

    elif p_diff:

        if len(f_pairs) != 2:
            assert False

        for f_pair in f_pairs:
            for f_split in f_pair.split(","):
                datum = np.zeros((m_bins, n_bins))
                datum = tabulate(datum, f_split, loci,
                                 p_start1, p_start2,
                                 p_end1, p_end2,
                                 p_chrom1, p_chrom2, p_resolution,
                                 p_downsample)
            raw.append(datum)

    else:  # accumulate

        datum = np.zeros((m_bins, n_bins))
        for f_pair in f_pairs:
            for f_split in f_pair.split(","):
                datum = tabulate(datum, f_split, loci,
                                 p_start1, p_start2,
                                 p_end1, p_end2,
                                 p_chrom1, p_chrom2, p_resolution,
                                 p_downsample)

        raw.append(datum)

    return raw

def normalize_contactmap(raw, p_mode="none"):
    """
    功能: 对接触矩阵进行归一化\n
    ​参数:
        raw: 原始接触矩阵列表
        p_mode: 归一化模式\n
    ​返回值: 归一化后的矩阵列表
    """
    norm = []
    for datum in raw:
        if p_mode == "self-relative":
            depth = datum.sum()//2
            datum = datum/float(depth)
            print("relative using self depth")
        elif p_mode == "KR":
            import scipy.sparse as sps
            from .jupyter_KR import knightRuizAlg, removeZeroDiagonalCSR, addZeroes

            #test = np.random.randint(0,20+1, size=(5,5))
            #test = (test + test.T)/2
            #test = sps.csr_matrix(test)
            #print(test.toarray())
            #DEBUG print(datum.sum(axis=0))
            sdatum = sps.csr_matrix(datum) #+1e-5)

            # remove sparse or diagonals with zeros
            percentOfSparseToRemove = 0
            mtxAndRemoved = removeZeroDiagonalCSR(sdatum, percentOfSparseToRemove)
            initialSize = sdatum.shape[0]
            sdatum, removed = mtxAndRemoved[0], mtxAndRemoved[1]
            newSize = sdatum.shape[0]

            #DEBUG print(sdatum.toarray().sum(axis=0))
            #DEBUG print(removed)

            CC = ((sdatum.sum())/(sdatum.shape[0]*2))
            CCother = sdatum.sum()/sdatum.size
            result = knightRuizAlg(sdatum)
            col = result[0]
            x = sps.diags(col.flatten(), 0, format='csr')
            mtx = x.dot(sdatum.dot(x))

            #print(mtx.toarray())
            CCmtx = CC * mtx
            CCothermtx = CCother * mtx

            # restore zeros
            datum = addZeroes(CCmtx.toarray(), removed)
            #datum = CCmtx.toarray()
            #datum = CCothermtx.toarray()

            print("balanced")
        elif p_mode == "VC":
            rowsum = datum.sum(axis=0)
            datum = datum / np.sqrt(np.outer(rowsum, rowsum))
            datum = np.nan_to_num(datum)
            print("VC")
        elif p_mode == "normalize":
            #depth = datum.sum()/2
            #datum = datum/float(depth)
            datum = normalize_by_diag(datum)
        norm.append(datum)
    return norm

def compute_differential(norm, p_diff="norm_diff"):
    """
    ​功能: 计算两个接触矩阵的差异\n
    ​参数:
        norm: 归一化后的矩阵列表
        p_diff: 差异计算方法\n
    ​返回值: 差异矩阵
    """
    pseudocount = 1e-10
    if p_diff == "log2FC":
        results = np.log2((norm[0]+pseudocount)/(norm[1]+pseudocount))
    elif p_diff == "diff":
        results = (norm[0]+pseudocount)-(norm[1]+pseudocount)
    elif p_diff == "pct_diff":
        results = ((norm[0]+pseudocount)-(norm[1]+pseudocount))/(norm[1]+pseudocount)
    elif p_diff == "norm_diff":
        results = normalize_diff(norm[0], norm[1])
    return results

def set_n_max(norm, p_max, p_mode):
    """
    ​功能: 设置颜色映射的最大值\n
    ​参数:
        norm: 归一化矩阵
        p_max: 用户指定的最大值
        p_mode: 归一化模式\n
    ​返回值: 计算得到的最大值
    """
    if p_max is not None:
        n_max = p_max
    else:
        if p_mode == "relative":
            n_max = max([datum.mean()*3 for datumn in norm])
        elif p_mode == "VC":
            n_max = 0.3
        else:
            n_max = 1
    return n_max

def set_p_mode(p_relative, p_self_relative, p_balance, p_vc, p_normalize, **kwargs):
    """
    功能: 根据参数设置归一化模式\n
    ​参数: 各种归一化选项\n
    ​返回值: 归一化模式字符串
    """
    p_mode = kwargs["p_mode"]
    if p_relative:
        p_mode = "relative"
    if p_self_relative:
        p_mode = "self-relative"
    if p_balance:
        p_mode = "KR"
    if p_vc:
        p_mode = "VC"
    if p_normalize:
        p_mode = "normalize"
    return p_mode

def set_dims(f_pairs, p_col, p_multi, **kwargs):
    """
    功能: 设置子图的行列数\n
    ​参数:
        f_pairs: 输入文件列表
        p_col: 列数
        p_multi: 是否多图模式\n
    ​返回值: (行数, 列数)
    """
    if p_multi:
        n_rows = ((len(f_pairs)-1)//p_col)+1
        n_cols = ((len(f_pairs)-1)%p_col)+1
    else:
        n_rows = 1
        n_cols = 1
    return n_rows, n_cols

def add_structure_overlay(ax,
                          f_dot, p_resolution, p_display, p_structure, f_fasta,
                          p_chrom1, p_start1, p_end1, p_dot_size, p_col, **kwargs):
    """
    功能: 
    提供dot文件时 : 在接触图上添加RNA二级结构叠加,#256baf(蓝色)表示外部结构数据,red表示预测的结构数据\n
    不提供dot : 蓝色表示预测结果数据 
    ​参数: 各种结构相关参数\n
    ​返回值: 更新后的axes对象
    """
    datum_structs = []
    color_structs = ["#256baf", "red"]
    #color_structs = ["purple", "red"]

    if f_dot:
        datum_structs.append(parse_structure(f_dot, p_resolution, p_display))
    if p_structure:
        o_fasta = Fasta(f_fasta)
        seq = str(o_fasta[p_chrom1][p_start1:p_end1]).upper()
        datum_structs.append(get_bppm(seq, p_res=p_resolution, p_display=p_display))

    for datum_struct, color in zip(datum_structs, color_structs):
        X = []
        Y = []
        size = []
        for i in range(datum_struct.shape[0]-1):
            for j in range(datum_struct.shape[1]-1):
                x = datum_struct[i,j]
                X.append(i)
                Y.append(j)
                size.append(x)
        size = np.array(size)
        size = size * p_dot_size

        # filter out super small points in the overlay
        size[size<0.001] = 0.0
        
        for i in range(ax.shape[0]*ax.shape[1]):
            r, c = i//p_col, i%p_col
            ax[r][c].scatter(X, Y, marker="o", s=size, c=color, alpha=1.0)
        
    return ax

def add_loop_overlay(ax,
                     p_chrom1, p_start1, p_end1,
                     f_loop, p_resolution, p_point, **kwargs):
    """
    ​功能: 在接触图上添加环状结构叠加\n
    ​参数: 各种环状结构相关参数\n
    ​返回值: 更新后的axes对象
    """
    dict_color = {"loops": "black", "lstripe": "blue", "rstripe": "green"}
    if f_loop is not None:
        p_boxwidth = 1.0

        loop = load_loops(p_chrom1, p_start1, p_end1, f_loop,
                          p_resolution, p_point=p_point)

        for l in loop:
            if len(l) == 5:
                lstart1, lend1, lstart2, lend2, category = l
                color = dict_color.get(category, "orange")
            else:
                lstart1, lend1, lstart2, lend2 = l
                color = "black" #"blue"

            loop_start1 = ((lstart1-p_start1)//p_resolution)-0.5
            loop_end1 = ((lend1-p_start1)//p_resolution)-0.5
            loop_start2 = ((lstart2-p_start1)//p_resolution)-0.5
            loop_end2 = ((lend2-p_start1)//p_resolution)-0.5

            def define_boxes(ax):
                ax.hlines(loop_start1-p_boxwidth, loop_start2-p_boxwidth,
                          loop_end2+p_boxwidth, color=color, alpha=1.0, linewidth=1.0)
                ax.hlines(loop_end1+p_boxwidth, loop_start2-p_boxwidth,
                          loop_end2+p_boxwidth, color=color, alpha=1.0, linewidth=1.0)

                ax.vlines(loop_start2-p_boxwidth, loop_start1-p_boxwidth,
                          loop_end1+p_boxwidth, color=color, alpha=1.0, linewidth=1.0)
                ax.vlines(loop_end2+p_boxwidth, loop_start1-p_boxwidth,
                          loop_end1+p_boxwidth, color=color, alpha=1.0, linewidth=1.0)
                return ax
            
            for i in range(ax.shape[0]*ax.shape[1]):
                r, c = i//p_col, i%p_col
                ax[r][c] = define_boxes(ax[r][c])
            
    return ax

def label_axes(ax, sizes, p_multi, p_label_res,
               p_chrom1, p_start1, p_end1,
               p_chrom2, p_start2, p_end2,
               bins1, bins2):
    """
    功能: 设置坐标轴标签\n
    ​参数: 各种区域和标签相关参数\n
    ​返回值: 更新后的axes对象
    """
    xloc = ax[0][0].get_xticks()
    yloc = ax[0][0].get_yticks()
    xlabel = ax[0][0].get_xticklabels()
    ylabel = ax[0][0].get_yticklabels()

    if p_start1 and p_end1 and p_start2 and p_end2:
        xnewloc = np.linspace(0, bins2, len(xloc)-1)
        ynewloc = np.linspace(0, bins1, len(yloc)-1)
        xnewlabel = [ "%.3f" % (i/p_label_res) for i in np.linspace(p_start2, p_end2, len(xloc)-1) ]
        ynewlabel = [ "%.3f" % (i/p_label_res) for i in np.linspace(p_start1, p_end1, len(yloc)-1) ]
    elif p_start1 and p_end1:
        xnewloc = np.linspace(0, bins1, len(xloc)-1)
        ynewloc = np.linspace(0, bins1, len(yloc)-1)
        xnewlabel = [ "%.1f" % (i/p_label_res) for i in np.linspace(p_start1, p_end1, len(xloc)-1) ]
        ynewlabel = [ "%.1f" % (i/p_label_res) for i in np.linspace(p_start1, p_end1, len(yloc)-1) ]
    elif p_chrom1 and p_chrom2:
        xnewloc = np.linspace(0, bins2, len(xloc)-1)
        ynewloc = np.linspace(0, bins1, len(yloc)-1)
        xnewlabel = [ "%.1f" % (i/p_label_res) for i in np.linspace(0, sizes[p_chrom2], len(xloc)-1) ]
        ynewlabel = [ "%.1f" % (i/p_label_res) for i in np.linspace(0, sizes[p_chrom1], len(yloc)-1) ]
    elif p_chrom1:
        xnewloc = np.linspace(0, bins1, len(xloc)-1)
        ynewloc = np.linspace(0, bins1, len(yloc)-1)
        xnewlabel = [ "%.1f" % (i/p_label_res) for i in np.linspace(0, sizes[p_chrom1], len(xloc)-1) ]
        ynewlabel = [ "%.1f" % (i/p_label_res) for i in np.linspace(0, sizes[p_chrom1], len(yloc)-1) ]

    plt.setp(ax, xticks=xnewloc, yticks=ynewloc,
             xticklabels=xnewlabel, yticklabels=ynewlabel)

    return ax

def create_axes(n_rows, n_cols, p_ax, p_figsize, p_dpi, **kwargs):
    """
    ​功能: 创建图形和axes对象\n
    ​参数:
        n_rows: 行数
        n_cols: 列数
        p_ax: 现有axes
        p_figsize: 图形大小
        p_dpi: 分辨率\n
    ​返回值: (fig, ax)元组
    """
    if p_ax is None:
        fig, ax = plt.subplots(n_rows, n_cols,
                               figsize=((p_figsize[0]*n_cols)+2, p_figsize[1]*n_rows),
                               dpi=p_dpi, squeeze=0)
    else:
        fig, ax = p_ax
        ax = np.array([[ax]])  # squeeze
    return fig, ax


def plot_contactmap(f_pairs, p_chrom1=None, p_chrom2=None,
                    p_start1=None, p_start2=None,
                    p_end1=None, p_end2=None, p_gene1=None, p_gene2=None,
                    p_resolution=100, p_max=None,
                    p_balance=False, p_relative=False, p_self_relative=False,
                    p_normalize=False, p_mode="none",
                    p_vc=False, p_multi=False, p_col=3,
                    p_genome=None, p_structure=False,
                    p_display="full",
                    f_size=None, f_out=None,
                    f_loop=None, f_dot=None, f_fasta=None,
                    f_index=None,
                    p_downsample=1.0, p_figsize=(4, 4), p_dot_size=5.0, p_point=False,
                    p_label_res=1000, p_dpi=100, p_ax=None, p_verbose=True):
    """
    # 功能概述
    plot_contactmap() 是一个综合性的绘图函数，主要功能包括：
        - 加载和预处理: 从输入的基因组交互数据 (如Hi-C数据)中提取指定区域的接触频率。
        - 归一化处理: 支持多种归一化方法 (如KR平衡、VC归一化等)。
        - 可视化: 生成热图形式的接触图, 并可叠加RNA二级结构、染色质环 (loops) 等注释信息。
        - 输出: 支持保存为图片文件或返回Matplotlib对象供进一步调整。
    # 核心参数说明
    1. 输入数据参数
        - f_pairs: 必选参数，输入文件的路径(pari.gz;支持多个文件(重复测序样本))。
        - p_chrom1/p_charm2: 染色体名称（如"chr1"），用于指定查询区域。
        - p_start1/p_end1: 起始和终止位置（坐标范围）。
        - p_gene1/p_gene2: 基因名称 (替代直接坐标,需配合f_index使用)。
    2. 数据处理参数
        - p_resolution: 分辨率 (bin大小, 默认100bp)。
        - p_downsample: 下采样比例 (0.0-1.0)，用于减少数据量。
        - p_mode: 归一化模式，支持:
            "none"（原始计数）
            "KR" (Knight-Ruiz平衡) 
            "VC"（方差校正）
            "self-relative"（相对自身深度归一化）
            "normalize"（按对角线距离归一化）
    3. 可视化参数
        - p_max: 颜色映射的最大值（自动计算或手动指定）。
        - p_figsize: 图像尺寸（默认(4, 4)）。
        - p_dpi: 图像分辨率 (默认100)。
        - f_out: 输出文件路径（如"output.png"）。
    4. 叠加注释参数
        - f_loop: 染色质环 (loops) 的BED文件路径, 用于在热图上标记环状结构。
        - f_dot: RNA二级结构文件 (点括号格式)，用于叠加碱基配对信息。
        - p_structure: 若为True, 从f_fasta序列预测RNA二级结构并叠加。
        - p_display: 显示模式 ("full"、"upper"或"lower"三角矩阵）。
    # ​工作流程\n
        1. 区域解析：
            调用get_loci()根据输入参数生成查询区域（如"chr1:1000000-2000000"）。
            若通过基因名称查询, 使用f_index文件解析坐标。
            ​数据加载与统计：
            调用get_contactmap()从输入文件中提取指定区域的交互频率，生成原始接触矩阵。
        2. ​归一化处理：
            调用normalize_contactmap()按指定方法 (如KR平衡) 归一化数据。
            ​绘图设置：
            根据数据量自动调整子图布局 (n_rows/n_cols)。
            创建Matplotlib的fig和ax对象。
        3. ​绘制热图：
            使用imshow()绘制接触矩阵，颜色映射为"Reds"（默认）或"coolwarm"（归一化差异时）。
            添加颜色条 (colorbar)。
        4. ​叠加注释：
            若有f_loop, 调用add_loop_overlay()用矩形框标记环状结构。
            若有f_dot或p_structure, 调用add_structure_overlay()用散点标记碱基配对。
        5. ​坐标轴调整：
            调用label_axes()将bin编号转换为实际基因组坐标 (如1.0Mb)。
        6. ​输出结果：
            保存图像 (若f_out指定) 或返回fig和ax对象。
    """
    if p_genome:
        f_size, f_fasta = get_ref(p_genome)

    dict_kwargs = {"f_pairs": f_pairs, "p_chrom1": p_chrom1, "p_chrom2": p_chrom2,
                   "p_start1": p_start1, "p_start2": p_start2,
                   "p_end1": p_end1, "p_end2": p_end2,
                   "p_gene1": p_gene1, "p_gene2": p_gene2,
                   "p_resolution": p_resolution, "p_max": p_max,
                   "p_balance": p_balance, "p_relative": p_relative, "p_self_relative": p_self_relative,
                   "p_normalize": p_normalize, "p_mode": p_mode,
                   "p_vc": p_vc, "p_multi": p_multi, "p_col": p_col, "p_genome": p_genome, "p_structure": p_structure,
                   "p_display": p_display, "f_size": f_size, "f_out": f_out, "f_loop": f_loop,
                   "f_dot": f_dot, "f_fasta": f_fasta, "f_index": f_index,
                   "p_downsample": p_downsample, "p_figsize": p_figsize,
                   "p_dot_size": p_dot_size, "p_point": p_point, "p_label_res": p_label_res,
                   "p_dpi": p_dpi, "p_ax": p_ax, "p_verbose": p_verbose}

    (loci, p_chrom1, p_start1, p_end1,
     p_chrom2, p_start2, p_end2,
     m_bins, n_bins, bins1, bins2) = get_loci(**dict_kwargs)

    sizes = get_size(f_size)    
    raw = get_contactmap(**dict_kwargs)
    p_mode = set_p_mode(**dict_kwargs)
    norm = normalize_contactmap(raw, p_mode=p_mode)
    n_max = set_n_max(norm, p_max, p_mode)
    n_rows, n_cols = set_dims(**dict_kwargs)
    fig, ax = create_axes(n_rows, n_cols, **dict_kwargs)
    
    if p_verbose:
        print(ax.shape)
        print(loci)        
        for datum in raw:
            print("Dim: %s | Total Counts: %d" % (datum.shape, datum.sum()//2))
        print(n_rows, n_cols)

    # ================
    # Plot Contactmap
    # ================
    for i in range(ax.shape[0]*ax.shape[1]):
        r, c = i//p_col, i%p_col
        indices = (r*p_col)+c
        if p_normalize:
            im = ax[r][c].imshow(norm[indices], cmap="coolwarm",
                                 vmin=-abs(float(n_max)), vmax=abs(float(n_max)))
        else:
            im = ax[r][c].imshow(norm[indices], cmap = "Reds", vmin=0, vmax=float(n_max))

        for tick in ax[r][c].get_xticklabels():
            tick.set_rotation(45)

    # ============================
    # Adjust fig and add colorbar
    # ============================
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax)
        
    # =======================
    # Label axes and ticks
    # =======================
    ax = label_axes(ax, sizes, p_multi, p_label_res,
                    p_chrom1, p_start1, p_end1,
                    p_chrom2, p_start2, p_end2,
                    bins1, bins2)
    
    # ======================================
    # Set up structure overlay onto heatmap
    # ======================================
    ax = add_structure_overlay(ax, **dict_kwargs)

    # ============
    # Plot Loops
    # ============
    ax = add_loop_overlay(ax, **dict_kwargs)

    if f_out:
        fig.savefig(f_out)

    return fig, ax


def plot_diffmap(f_pairs, p_chrom1=None, p_chrom2=None,
                 p_start1=None, p_start2=None,
                 p_end1=None, p_end2=None, p_gene1=None, p_gene2=None,
                 p_resolution=100, p_max=None,
                 p_balance=False, p_relative=False, p_self_relative=False,
                 p_normalize=False, p_mode="none",
                 p_vc=False, p_multi=False, p_col=3,
                 p_diff="norm_diff", p_genome=None, p_structure=False,
                 p_display="full",
                 f_size=None, f_out=None,
                 f_loop=None, f_dot=None, f_fasta=None,
                 f_index=None,
                 p_downsample=1.0, p_figsize=(4, 4), p_dot_size=5.0, p_point=False,
                 p_label_res=1000, p_dpi=100, p_ax=None, p_verbose=True):
    """
    plot_diffmap 是一个用于可视化基因组交互差异 (differential genomic interactions) 的函数。
    1. 函数功能
        获取基因组区域信息
        加载接触矩阵 (contact map) 数据
        标准化接触矩阵
        计算差异交互模式
        绘制差异热图
        添加各种注释和覆盖层（如结构信息、染色质环等）
    2. 参数详解
    主要输入参数
        **f_pairs**: 输入文件路径，包含交互对数据
    ​染色体/区域参数:
        - p_chrom1, p_chrom2: 染色体名称
        - p_start1, p_start2, p_end1, p_end2: 区域起止位置
        - p_gene1, p_gene2: 基因名称（可替代坐标）
    ​分辨率参数:
        p_resolution: 分析分辨率 (默认100bp) 
    ​显示控制参数:
        - p_max: 颜色标尺的最大值
        - p_balance: 是否平衡矩阵
        - p_relative: 是否显示相对值
        - p_normalize: 是否标准化
        - p_mode: 标准化模式
    ​输出控制参数:
        - f_out: 输出文件路径
        - p_figsize: 图形大小
        - p_dpi: 输出分辨率
    """
    if p_genome:
        f_size, f_fasta = get_ref(p_genome)
    
    dict_kwargs = {"f_pairs": f_pairs, "p_chrom1": p_chrom1, "p_chrom2": p_chrom2,
                   "p_start1": p_start1, "p_start2": p_start2,
                   "p_end1": p_end1, "p_end2": p_end2,
                   "p_gene1": p_gene1, "p_gene2": p_gene2,
                   "p_resolution": p_resolution, "p_max": p_max,
                   "p_balance": p_balance, "p_relative": p_relative, "p_self_relative": p_self_relative,
                   "p_normalize": p_normalize, "p_mode": p_mode,
                   "p_vc": p_vc, "p_multi": p_multi, "p_col": p_col, "p_diff": p_diff, "p_genome": p_genome, "p_structure": p_structure,
                   "p_display": p_display, "f_size": f_size, "f_out": f_out, "f_loop": f_loop,
                   "f_dot": f_dot, "f_fasta": f_fasta, "f_index": f_index,
                   "p_downsample": p_downsample, "p_figsize": p_figsize,
                   "p_dot_size": p_dot_size, "p_point": p_point, "p_label_res": p_label_res,
                   "p_dpi": p_dpi, "p_ax": p_ax, "p_verbose": p_verbose}

    (loci, p_chrom1, p_start1, p_end1,
     p_chrom2, p_start2, p_end2,
     m_bins, n_bins, bins1, bins2) = get_loci(**dict_kwargs)
    
    sizes = get_size(f_size)    
    raw = get_contactmap(**dict_kwargs)
    p_mode = set_p_mode(**dict_kwargs)
    norm = normalize_contactmap(raw, p_mode=p_mode)
    n_max = set_n_max(norm, p_max, p_mode)
    n_rows, n_cols = set_dims(**dict_kwargs)
    fig, ax = create_axes(n_rows, n_cols, **dict_kwargs)
    
    if p_verbose:
        print(loci)        
        for datum in raw:
            print("Dim: %s | Total Counts: %d" % (datum.shape, datum.sum()//2))
        print(n_rows, n_cols)

    # =================
    # Differential map
    # =================
    results = compute_differential(norm, p_diff=p_diff)
    im = ax[0][0].imshow(results, cmap="coolwarm", vmin=-abs(float(n_max)),
                         vmax=abs(float(n_max)))
    for tick in ax[0][0].get_xticklabels():
        tick.set_rotation(45)
    
    # ============================
    # Adjust fig and add colorbar
    # ============================
    fig.subplots_adjust(right=0.8)
    cbar_ax = fig.add_axes([0.85, 0.15, 0.02, 0.7])
    fig.colorbar(im, cax=cbar_ax)
        
    # =======================
    # Label axes and ticks
    # =======================
    ax = label_axes(ax, sizes, p_multi, p_label_res,
                    p_chrom1, p_start1, p_end1,
                    p_chrom2, p_start2, p_end2,
                    bins1, bins2)
    
    # ======================================
    # Set up structure overlay onto heatmap
    # ======================================
    ax = add_structure_overlay(ax, **dict_kwargs)

    # ============
    # Plot Loops
    # ============
    ax = add_loop_overlay(ax, **dict_kwargs)

    if f_out:
        fig.savefig(f_out)

    return fig, ax


def plot_contactmap_w_arcbands(f_pairs, f_cluster, p_chrom1=None, p_chrom2=None,
                               p_start1=None, p_start2=None,
                               p_end1=None, p_end2=None, p_gene1=None, p_gene2=None,
                               p_resolution=100, p_max=None,
                               p_balance=False, p_relative=False, p_self_relative=False,
                               p_normalize=False, p_mode="none",
                               p_vc=False, p_multi=False, p_col=3,
                               p_genome=None, p_structure=False,
                               p_display="full",
                               f_size=None, f_out=None,
                               f_loop=None, f_dot=None, f_fasta=None,
                               f_index=None,
                               p_downsample=1.0, p_figsize=(4, 4), p_dot_size=5.0, p_point=False,
                               p_label_res=1000, p_dpi=100, p_ax=None, p_verbose=True, p_order=None):
    """
    ## 主要功能
    - 绘制基因组接触热图
    - 在热图上方添加弧形条带注释
    - 支持多种数据标准化方式
    - 可添加染色质环和结构注释
    ## 参数解析
    ### 核心输入参数
    - **f_pairs**: Hi-C交互数据文件路径
    - **f_cluster**: 结构域/聚类信息文件路径（用于弧形条带）
    ### ​基因组坐标参数:
    - p_chrom1, p_chrom2: 染色体名称
    - p_start1, p_end1: 区域起始和结束位置
    - p_gene1, p_gene2: 基因名称（替代坐标）
    ### 可视化控制参数
    - **p_resolution**: 分析分辨率 (默认100bp) 
    - **p_max**: 颜色标尺最大值
    - **p_balance**: 是否进行矩阵平衡
    - **p_normalize**: 是否标准化数据
    - **p_mode**: 标准化模式
    - **p_figsize**: 图形尺寸
    - **p_dpi**: 输出分辨率
    """
    if p_genome:
        f_size, f_fasta = get_ref(p_genome)
    
    dict_kwargs = {"f_pairs": f_pairs, "f_cluster": f_cluster, "p_chrom1": p_chrom1, "p_chrom2": p_chrom2,
                   "p_start1": p_start1, "p_start2": p_start2,
                   "p_end1": p_end1, "p_end2": p_end2,
                   "p_gene1": p_gene1, "p_gene2": p_gene2,
                   "p_resolution": p_resolution, "p_max": p_max,
                   "p_balance": p_balance, "p_relative": p_relative, "p_self_relative": p_self_relative,
                   "p_normalize": p_normalize, "p_mode": p_mode,
                   "p_vc": p_vc, "p_multi": p_multi, "p_col": p_col, "p_genome": p_genome,
                   "p_structure": p_structure,
                   "p_display": p_display, "f_size": f_size, "f_out": f_out, "f_loop": f_loop,
                   "f_dot": f_dot, "f_fasta": f_fasta, "f_index": f_index,
                   "p_downsample": p_downsample, "p_figsize": p_figsize,
                   "p_dot_size": p_dot_size, "p_point": p_point, "p_label_res": p_label_res,
                   "p_dpi": p_dpi, "p_ax": p_ax, "p_verbose": p_verbose, "p_order": p_order}

    (loci, p_chrom1, p_start1, p_end1,
     p_chrom2, p_start2, p_end2,
     m_bins, n_bins, bins1, bins2) = get_loci(**dict_kwargs)
    
    sizes = get_size(f_size)    
    raw = get_contactmap(**dict_kwargs)
    p_mode = set_p_mode(**dict_kwargs)
    norm = normalize_contactmap(raw, p_mode=p_mode)
    n_max = set_n_max(norm, p_max, p_mode)
    n_rows, n_cols = set_dims(**dict_kwargs)

    import matplotlib.gridspec as grd
    fig = plt.figure(figsize=(10, 10))
    gs = grd.GridSpec(2, 2,
                      height_ratios=[2, 8],
                      width_ratios=[8, 2],
                      wspace=0.1)

    ax = np.array([[plt.subplot(gs[2])]])  # contactmap
    ax_arcband = plt.subplot(gs[0])  # arcband
    ax_colorbar = plt.subplot(gs[3])  # colorbar
    
    
    if p_verbose:
        print(ax.shape)
        print(loci)        
        for datum in raw:
            print("Dim: %s | Total Counts: %d" % (datum.shape, datum.sum()//2))
        print(n_rows, n_cols)

    # ================
    # Plot Contactmap
    # ================
    for i in range(ax.shape[0]*ax.shape[1]):
        r, c = i//p_col, i%p_col
        indices = (r*p_col)+c
        if p_normalize:
            im = ax[r][c].imshow(norm[indices], cmap="coolwarm",
                                 vmin=-abs(float(n_max)), vmax=abs(float(n_max)))
        else:
            im = ax[r][c].imshow(norm[indices], cmap = "Reds",
                                 vmin=0, vmax=float(n_max), aspect="auto")

        for tick in ax[r][c].get_xticklabels():
            tick.set_rotation(45)

    # ======================================
    # Add colorbar associated w contact map
    # ======================================
    fig.colorbar(im, cax=ax_colorbar)
    
    # =======================
    # Label axes and ticks
    # =======================
    ax = label_axes(ax, sizes, p_multi, p_label_res,
                    p_chrom1, p_start1, p_end1,
                    p_chrom2, p_start2, p_end2,
                    bins1, bins2)

    # =======================
    # Set up ArcBand
    # =======================
    from .utils_arcband import add_arcband_overlay
    ax_arcband = add_arcband_overlay(f_cluster, sizes,
                                     ax_arcband, p_chrom1, p_order)
    
    # ======================================
    # Set up structure overlay onto heatmap
    # ======================================
    ax = add_structure_overlay(ax, **dict_kwargs)

    # ============
    # Plot Loops
    # ============
    ax = add_loop_overlay(ax, **dict_kwargs)

    if f_out:
        fig.savefig(f_out)

    return fig, ax
