import re
import gffutils
import pandas as pd
from collections import defaultdict

def get_transcript_structure(transcript_id, db):
    """
    根据转录本ID获取其CDS和外显子结构
    返回: {
        'cds_ranges': [(start1, end1), ...],
        'exon_ranges': [(start1, end1), ...],
        'strand': '+/-'
    }
    """
    try:
        transcript = db[transcript_id]
    except:
        return None
    ### mRNA 有起始密码子和终止密码子
    start_codon = list(db.children(transcript_id, featuretype="start_codon", order_by="start"))
    stop_codon = list(db.children(transcript_id, featuretype="stop_codon", order_by="start"))
    if len(start_codon) !=0 and len(stop_codon) !=0:
        start_codonT = [s.start for s in start_codon] # 位于第一个外显子上
        stop_codonT = [s.end for s in stop_codon] # 位于最后一个外显子上

        exon_list = list(db.children(transcript_id, featuretype="exon", order_by="start"))
        exon_ranges = [(e.start,e.end) for e in exon_list]
        exon_lenghts = [e.end - e.start for e in exon_list]

        if transcript.strand == '+':
            UTR3 = start_codonT[0] - exon_ranges[0][0]
            UTR5 = sum(exon_lenghts) - (exon_ranges[len(exon_ranges)-1][1] - stop_codonT[0])
        else:
            UTR3 = stop_codonT[0] - exon_ranges[0][0]
            UTR5 = sum(exon_lenghts) - (exon_ranges[len(exon_ranges)-1][1] - start_codonT[0])
    else:
        UTR3=0
        UTR5=9999999999999999 # 一个非常大的数字

    # print(exon_ranges)
    # print(start_codonT)
    # print(stop_codonT)
    
    ### 其它RNA只有外显子

    
    # 获取外显子区域

    
    return {
        "UTR3": UTR3,
        "UTR5": UTR5
    }

def annotate_position_in_transcript(transcript_id, pos, db, label_type="t"):
    """
    根据转录本结构和位置返回对应的label
    :param transcript_id: 转录本ID
    :param pos: 位置 (1-based)
    :param db: gffutils数据库
    :param label_type: "t" (target) 或 "c" (control)
    :return: label字符串
    """

    structure = get_transcript_structure(transcript_id, db)
    if not structure:
        return f"Intergenic ({label_type})"
    # print(structure)
    pos = int(pos)
    UTR3 = structure["UTR3"]
    UTR5 = structure["UTR5"]
    if UTR3 < pos < UTR5:
        return f"Exon ({label_type})"
    elif pos < UTR3:
        return f"3' UTR ({label_type})"
    else:
        return f"5' UTR ({label_type})"

def main(input_file, output_file, db_path):
    # 加载数据库
    db = gffutils.FeatureDB(db_path)
    
    # 读取输入文件
    with open(input_file) as f_in, open(output_file, "w") as f_out:
        # 写入表头
        f_out.write("Read_ID\tTarget_Transcript\tTarget_Pos\tTarget_Label\tControl_Transcript\tControl_Pos\tControl_Label\n")
        
        for line in f_in:
            if line.startswith("#"):
                continue
            
            parts = line.strip().split()
            if len(parts) < 5:
                continue
                
            read_id = parts[0]
            l_transcript = parts[1]
            l_pos = parts[2]

            r_transcript = parts[3]
            r_pos = parts[4]

            # 分子内相互作用为cis，分子间相互作用为trans
            if l_transcript == r_transcript:
                l_label = annotate_position_in_transcript(l_transcript, l_pos, db, "cis")
                r_label = annotate_position_in_transcript(r_transcript, r_pos, db, "cis")
            else:
                l_label = annotate_position_in_transcript(l_transcript, l_pos, db, "trans")
                r_label = annotate_position_in_transcript(r_transcript, r_pos, db, "trans")
            
            # 输出结果
            f_out.write(f"{read_id}\t{l_transcript}\t{l_pos}\t{l_label}\t{r_transcript}\t{r_pos}\t{r_label}\n")

if __name__ == '__main__':
    gtf_file = "/ChIP_seq_2/Data/index/Mus_musculus/NCBI/GRCm38p6/GCF_000001635.26_GRCm38.p6_genomic.gtf"
    db_file = "NCBI.db"
    # db = gffutils.FeatureDB(db_file)
    # print(get_transcript_structure('NM_001001130.3',db))
    # print(db['XM_006495550.4'])
    # # 创建数据库（只需运行一次）
    # db = gffutils.create_db(
    #     gtf_file, 
    #     dbfn=db_file, 
    #     force=True,
    #     keep_order=True,
    #     merge_strategy="merge",
    #     sort_attribute_values=True,
    #     disable_infer_genes=True
    # )
    
    # 运行主函数
    main("a.csv", "b.csv", db_file)