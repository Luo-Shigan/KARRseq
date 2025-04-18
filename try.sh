#!/bin/bash
### conda: KARRseq
function merge_fq(){
    read1=$1
    read2=$2
    outRead1=$3
    outRead2=$4
    outRead3=$5
    outRead4=$6
    outMerge=$7
    SeqPrep \
        -f ${read1} \
        -r ${read2} \
        -1 ${outRead1} \
        -2 ${outRead2} \
        -3 ${outRead3} \
        -4 ${outRead4} \
        -s ${outMerge}

}
export -f merge_fq
read1=data/GSE166155/fq/SRR13628276_1.fastq.gz
read2=data/GSE166155/fq/SRR13628276_2.fastq.gz
outRead1=output/fq/merge/SRR13628276/SRR13628276_unmerge_P1.fastq.gz
outRead2=output/fq/merge/SRR13628276/SRR13628276_unmerge_P2.fastq.gz
outRead3=output/fq/merge/SRR13628276/SRR13628276_reject_P1.fastq.gz
outRead4=output/fq/merge/SRR13628276/SRR13628276_reject_P2.fastq.gz
outMerge=output/fq/merge/SRR13628276/SRR13628276_merge.fastq.gz
# merge_fq ${read1} ${read2} \
#     ${outRead1} ${outRead2} \
#     ${outRead3} ${outRead4} \
#     ${outMerge}
function index(){
    genome_index=$1
    genome_fa=$2
    genome_gtf=$3
    STAR=/opt/STAR-2.5.2a/bin/Linux_x86_64/STAR
    ${STAR} --runMode genomeGenerate \
    --runThreadN 25 \
    --genomeDir ${genome_index} \
    --genomeFastaFiles ${genome_fa} \
    --sjdbGTFfile ${genome_gtf}
}
export -f index
genome_index=/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/star2.5.2a
genome_fa=/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/GRCm39.primary_assembly.genome.fa
genome_gtf=/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/gencode.vM36.primary_assembly.annotation.gtf
# index ${genome_index} ${genome_fa} ${genome_gtf}

function align(){
    index=$1
    fq=$2
    outPrefix=$3
    STAR=/opt/STAR-2.5.2a/bin/Linux_x86_64/STAR
    ${STAR} \
    --runMode alignReads \
    --runThreadN 25 \
    --genomeDir ${index} \
    --readFilesIn ${fq} \
    --readFilesCommand zcat \
    --outFileNamePrefix ${outPrefix}_ \
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
    --limitIObufferSize 1500000000

    samtools index ${outPrefix}_Aligned.sortedByCoord.out.bam
}
export -f align
index=/ChIP_seq_2/Data/index/Mus_musculus/GENCODE/GRCm39/star2.5.2a
fq=/ChIP_seq_2/StemCells/RNAInteraction20250327/output/fq/merge/SRR13628276/SRR13628276_merge.fastq.gz
outPrefix=output/bam/mouse/SRR13628276/SRR13628276
# align ${index} ${fq} ${outPrefix}

function ssecrtp(){
    aligned=$1
    mapq=$2
    span=$3
    tmp1=$4
    tmp2=$5
    chimeric=$6
    output=$7
    samtools view ${aligned} \
        | python workflow/KARR-seq/src/get_STAR_reads.py Aligned ${mapq} ${span} > ${tmp1}
        cat ${chimeric} \
        | python workflow/KARR-seq/src/get_STAR_reads.py Chimeric ${mapq} ${span} > ${tmp2}
        cat ${tmp1} ${tmp2} \
        | /usr/bin/sort -k2,2 -k4,4 -k3,3n -k5,5n -k10,10n -k11,11n --parallel=4 \
        | gzip -c > ${output}
}
export -f ssecrtp
aligned=/ChIP_seq_2/StemCells/RNAInteraction20250327/output/bam/mouse/SRR13628276/SRR13628276_Aligned.sortedByCoord.out.bam
mapq=1
span=0
tmp1=output/bam/mouse/SRR13628276/SRR13628276.tmp1
tmp2=output/bam/mouse/SRR13628276/SRR13628276.tmp2
chimeric=/ChIP_seq_2/StemCells/RNAInteraction20250327/output/bam/mouse/SRR13628276/SRR13628276_Chimeric.out.sam
output=output/bam/mouse/SRR13628276/MAPQ1_SPAN0.txt.gz
# ssecrtp ${aligned} ${mapq} ${span} ${tmp1} ${tmp2} ${chimeric} ${output}
function sserd(){
    input=$1
    todedup=$2
    dedup=$3
    mapq=$4
    bed=$5
    pairs=$6
    zcat ${input} | python workflow/KARR-seq/src/remove_duplicates.py ${todedup} | gzip -c > ${dedup}
    zcat ${dedup} | python workflow/KARR-seq/src/pairs_to_bed.py ${mapq} > ${bed}
    zcat ${dedup} | cut -f1-7 | bgzip -c > ${pairs}
    sleep 120
    pairix ${pairs}
}
export -f sserd
input=output/bam/mouse/SRR13628276/MAPQ1_SPAN0.txt.gz
todedup=dedup
dedup=output/bam/mouse/SRR13628276/MAPQ1_SPAN0.dedup.txt.gz
mapq=1
bed=output/bam/mouse/SRR13628276/MAPQ1_SPAN0.dedup.bed
pairs=output/bam/mouse/SRR13628276/MAPQ1_SPAN0.dedup.pairs.gz
# sserd ${input} ${todedup} ${dedup} ${mapq} ${bed} ${pairs}
function ssel(){
    input=$1
    output=$2
    zcat ${input} \
        | awk '{{print $1"\t"$2"\t"$3+$10"\t"$4"\t"$5"\t"$6"\t"$7}}' \
        | awk '{{if($2==$4){{if($3>$5){{print $1"\t"$4"\t"$5"\t"$2"\t"$3"\t"$7"\t"$6}}else{{print $0}}}}else{{print $0}}}}' \
        | sort -k2,2 -k4,4 -k3,3n -k5,5n --parallel=4 \
        | bgzip -c > ${output}
    sleep 60
    pairix ${output}
}
export -f ssel
input=output/bam/mouse/SRR13628276/MAPQ1_SPAN0.dedup.txt.gz
output=output/bam/mouse/SRR13628276/inner/MAPQ1_SPAN0.dedup.pairs_inner.gz
ssel ${input} ${output}
