shell.prefix("set -x; set -e;")
configfile: "config/KARRseq.yaml"

indir = config.get('indir', '../data')  # 如果没有传递 indir，则使用 'data' 作为默认值
outdir = config.get('outdir', '../output')  # 如果没有传递 outdir，则使用 'output' 作为默认值

import glob
from os.path import join
import re


SAMPLESPre = glob_wildcards(indir + "/{sample}.fastq.gz").sample
paired_samples = glob_wildcards(indir + "/{sample}_1.fastq.gz").sample
single_samples = [sample for sample in SAMPLESPre if re.match(r"^SRR\d+$", sample)]
SAMPLES = paired_samples + single_samples
print(SAMPLES)

# =================
# Targets
# =================


def request():
    ##对于能通过依赖关系寻找的中间文件不需要重复定义，否则执行次数会过多
    output = dict()
    output['star_single_end_ligation'] = expand(outdir + "/chimeric/inner/{genome}/{sample}/{sample}.dedup.pairs.gz",sample=SAMPLES,genome=['mm10_transcript'])
    return list(output.values())

rule targets_all:
    input:
        request()


# =====================
# Rules for Single-end (STAR)
# =====================


rule merge_fastq:
    input:
        read1 = indir + "/fq/{sample}_1.fastq.gz",
        read2 = indir + "/fq/{sample}_2.fastq.gz"
    output:
        read1 = outdir + "/merge/{sample}/{sample}_unmerge_1.fastq.gz",
        read2 = outdir + "/merge/{sample}/{sample}_unmerge_2.fastq.gz",
        read3 = outdir + "/merge/{sample}/{sample}_reject_1.fastq.gz",
        read4 = outdir + "/merge/{sample}/{sample}_reject_2.fastq.gz",
        merge = outdir + "/merge/{sample}/{sample}_merge.fastq.gz",
        report = outdir + "/merge/{sample}/{sample}_merge.report"
    conda:
        config['conda']['KARRseq']
    log:
        log = outdir + "/log/{sample}/merge_fastq.log"
    shell:
        """
	    SeqPrep \
            -f {input.read1} \
            -r {input.read2} \
            -1 {output.read1} \
            -2 {output.read2} \
            -3 {output.read3} \
            -4 {output.read4} \
            -s {output.merge} 2> {log.log}
        seqkit stats {output.merge} > {output.report} 2>> {log.log}
        """
def get_alignment_input(wildcards):
    """动态判断输入文件类型"""
    # 双端文件需要合并成单端；对于单端文件没必要
    paired = f"{outdir}/merge/{wildcards.sample}/{wildcards.sample}_merge.fastq.gz"
    single = f"{indir}/{wildcards.sample}.fastq.gz"
    
    # 检查文件实际存在情况
    if wildcards.sample in paired_samples:
        print(f"双端测序， 已合并: {[paired]}")
        return [paired]
    elif wildcards.sample in single_samples:
        print(f"单端测序：{[single]}")
        return [single]
    else:
        raise FileNotFoundError(
            f"Missing input files for sample {wildcards.sample}\n"
            f"Checked paths:\n- {paired}\n- {single}"
        )


rule align_with_STAR:
    input: 
        get_alignment_input
    output:
        aligned = outdir + "/bam/{genome}/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        chimeric = outdir + "/bam/{genome}/{sample}/{sample}_Chimeric.out.sam"
    log:
        log = outdir + "/log/{genome}/{sample}/STAR.log"
    params:
        prefix = outdir + "/bam/{genome}/{sample}/{sample}",
        index = lambda wildcards: config['index']['star'][wildcards.genome],
        STAR = config['procedure']['STAR']
    threads: 15
    shell:
        """
        {params.STAR} \
            --runMode alignReads \
            --runThreadN {threads} \
            --genomeDir {params.index} \
            --readFilesIn {input} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.prefix}_ \
            --outReadsUnmapped Fastq \
            --outFilterMultimapNmax 100 \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes All \
            --alignIntronMin 1 \
            --scoreGapNoncan -4 \
            --scoreGapATAC -4 \
            --chimSegmentMin 15 \
            --chimJunctionOverhangMin 15 \
            --limitOutSJcollapsed 10000000 \
            --limitIObufferSize 1500000000 > {log.log}
        samtools index {output.aligned} 2>> {log.log}
        """

rule star_single_end_chimeric_reads_to_pairs:
    input:
        aligned = outdir + "/bam/{genome}/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
        chimeric = outdir + "/bam/{genome}/{sample}/{sample}_Chimeric.out.sam"
    output: 
        outdir + "/chimeric/pairs/{genome}/{sample}/{sample}.txt.gz"
    log:
        outdir + "/log/{genome}/{sample}/star_single_end_chimeric_reads_to_pairs.log"
    params:
        mapq = 1,
        span = 0,
        tmp1 = outdir + "/chimeric/pairs/{genome}/{sample}/{sample}.tmp1",
        tmp2 = outdir + "/chimeric/pairs/{genome}/{sample}/{sample}.tmp2"
    threads: 4
    conda:
        config['conda']['KARRseq']
    shell:
        """
        samtools view {input.aligned} \
            | python src/get_STAR_reads.py Aligned {params.mapq} {params.span} > {params.tmp1}
        cat {input.chimeric} \
            | python src/get_STAR_reads.py Chimeric {params.mapq} {params.span} > {params.tmp2}
        cat {params.tmp1} {params.tmp2} \
            | /usr/bin/sort -k2,2 -k4,4 -k3,3n -k5,5n -k10,10n -k11,11n --parallel={threads} \
            | gzip -c > {output} 2>{log}
        rm {params.tmp1} {params.tmp2}
        """

rule star_single_end_remove_duplicates:
    input: 
        outdir + "/chimeric/pairs/{genome}/{sample}/{sample}.txt.gz"
    output:
        dedup = outdir + "/chimeric/pairs/{genome}/{sample}/{sample}.dedup.txt.gz",
        bed = outdir + "/chimeric/pairs/{genome}/{sample}/{sample}.dedup.bed",
        pairs = outdir + "/chimeric/pairs/{genome}/{sample}/{sample}.dedup.pairs.gz"
    log:
        outdir + "/log/{genome}/{sample}/star_single_end_remove_duplicates.log"
    params:
        mapq = 1,#与star_single_end_chimeric_reads_to_pairs规则中的mapq参数一致
        todedup = "dedup"
    conda:
        config['conda']['KARRseq']
    shell:
        """
        zcat {input} | python src/remove_duplicates.py {params.todedup} | gzip -c > {output.dedup}
        zcat {output.dedup} | python src/pairs_to_bed.py {params.mapq} > {output.bed}
        zcat {output.dedup} | cut -f1-7 | bgzip -c > {output.pairs} 2> {log}
        sleep 120
        pairix {output.pairs} > {log} 2>&1
        """

rule star_single_end_ligation:
    input: 
        outdir + "/chimeric/pairs/{genome}/{sample}/{sample}.dedup.txt.gz",
    output: 
        outdir + "/chimeric/inner/{genome}/{sample}/{sample}.dedup.pairs.gz"
    log:
        outdir + "/log/{genome}/{sample}/star_single_end_ligation.log"
    threads: 4
    conda:
        config['conda']['KARRseq']
    shell:
        """
        zcat {input} \
            | awk '{{print $1"\t"$2"\t"$3+$10"\t"$4"\t"$5"\t"$6"\t"$7}}' \
            | awk '{{if($2==$4){{if($3>$5){{print $1"\t"$4"\t"$5"\t"$2"\t"$3"\t"$7"\t"$6}}else{{print $0}}}}else{{print $0}}}}' \
            | sort -k2,2 -k4,4 -k3,3n -k5,5n --parallel={threads} \
            | bgzip -c > {output}
        sleep 60
        pairix {output}
        """
# 步骤	列顺序（以$1-$10为例）	说明
# ​原始输入	$1 $2 $3 $4 $5 $6 $7 $8 $9 $10	压缩文件中的原始列。
# ​第一个awk后	$1 $2 ($3+$10) $4 $5 $6 $7	第3列变为$3+$10，其余列保留。
# ​第二个awk后	可能交换$2/$4和$3/$5（若条件满足）	确保坐标按升序排列。
# ​最终输出	排序后的重组列	按染色体和坐标排序的BGZF文件。
# 如果作用对在同一个reference name（chr/或transcript name）;确保作用对第一个mapping起始位置小于第二个
# read id, 
# first reference name,first continous mapping segment end position
# second reference name,second continous mapping segment begining position, 
#first strand, second strand