# Description

本流程是对[KARR-seq](https://github.com/ouyang-lab/KARR-seq.git)流程输入输出位置做了些修改，使得其适应个性化分析任务

## 使用方法

举例：

```shell
 snakemake -s workflow/KARR-seq/KARRseq.smk --use-conda --cores 30 --directory workflow/KARR-seq/ --config indir="../../data/GSE166155" outdir="../../output" 
```

- --directory指定snakfile所在目录，使得其可以加载同目录下"config/KARRseq.yaml"
- --config indir="../../data/GSE166155" outdir="../../output" 分别指定输入和输出位置，相对于snakfile所在目录
