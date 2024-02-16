#!/bin/bash
# 更新包列表可能不再需要，因为所有文件都已经在本地

# 参考基因组和注释文件的路径（这些文件应该已经被手动下载和解压到指定位置）
REF_GENOME="/home/ubuntu/workspace/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
ANNOTATION_GTF="/home/ubuntu/workspace/Homo_sapiens.GRCh38.104.gtf"

# HISAT2索引的路径（假设已经解压到指定位置）
HISAT2_INDEX="~/workspace/grch38/genome"

# 遍历指定的SRA文件目录
find ~/workspace/DATA -type f -name "*.sra" | while read SRA_FILE; do
    # 获取SRA文件的目录路径
    SRA_DIR=$(dirname "$SRA_FILE")
    # 转换为fastq格式
    echo "Converting $SRA_FILE to fastq format..."
    fastq-dump --split-files -O $SRA_DIR $SRA_FILE

    # 获取不带扩展名的文件名
    BASENAME=$(basename "$SRA_FILE" .sra)

    # # 质量控制
    # echo "Performing quality control with FastQC..."
    # fastqc -o $SRA_DIR ${SRA_DIR}/${BASENAME}_1.fastq

    # 修剪fastq文件，去除adapter
    echo "Trimming fastq files with fastp..."
    fastp -i ${SRA_DIR}/${BASENAME}_1.fastq -o ${SRA_DIR}/${BASENAME}_1_trimmed.fq -h ${SRA_DIR}/${BASENAME}_fp.html -j ${SRA_DIR}/${BASENAME}_fp.json

    # 映射到参考基因组
    echo "Mapping reads to the reference genome with HISAT2..."
    hisat2 -p 8 -x $HISAT2_INDEX -U ${SRA_DIR}/${BASENAME}_1_trimmed.fq -S ${SRA_DIR}/${BASENAME}_aligned.sam --dta >${SRA_DIR}/${BASENAME}_hisat2_log.log 2>&1

    # 将SAM转换为BAM并排序
    echo "Converting SAM to sorted BAM..."
    samtools sort -@ 4 -o ${SRA_DIR}/${BASENAME}_aligned.bam ${SRA_DIR}/${BASENAME}_aligned.sam

    # 计算基因表达量
    echo "Calculating gene expression with StringTie..."
    stringtie -p 10 -G $ANNOTATION_GTF -o ${SRA_DIR}/${BASENAME}.gtf -l ${BASENAME} -A ${SRA_DIR}/${BASENAME}_gene_abundances.tsv -B ${SRA_DIR}/${BASENAME}_aligned.bam

    # 生成VCF文件
    # echo "Generating VCF file for $SRA_FILE..."
    # bcftools mpileup -Ou -f $REF_GENOME ${SRA_DIR}/${BASENAME}_aligned.bam | bcftools call -Ov -mv | bcftools filter -i 'QUAL > 30' -o ${SRA_DIR}/${BASENAME}_filtered.vcf

    # 将基因表达量汇总.tsv文件转换为.csv
    echo "Converting gene abundance TSV to CSV format..."
    cat ${SRA_DIR}/${BASENAME}_gene_abundances.tsv | tr '\t' ',' > ${SRA_DIR}/${BASENAME}_gene_abundances.csv

    # 删除文件
    find "${SRA_DIR}" -type f | grep -vP "(${BASENAME}_aligned\.bam|${BASENAME}\.sra|\.vcf$|\.tsv$|\.csv$|\.log$)" | xargs rm

    echo "Analysis complete for $SRA_FILE."
done

echo "All analyses complete."
