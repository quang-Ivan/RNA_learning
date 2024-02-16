#!/user/bin/bash
# bash解释器可能要根据自己的调一下
# 更新包列表
sudo apt-get update

# 设置工作目录（根据自己的目录给就行）
WORKDIR="~/workspace/"
cd $WORKDIR

# 下载参考基因组和注释文件（这些步骤只需执行一次）
echo "Downloading reference genome and annotation..."
wget -nc ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget -nc ftp://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz

# 解压参考基因组和注释文件（如果尚未解压）
echo "Decompressing reference genome and annotation..."
gunzip -k -f Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip -k -f Homo_sapiens.GRCh38.104.gtf.gz

# 下载并解压HISAT2索引（这些步骤只需执行一次）
echo "Downloading and extracting HISAT2 index..."
wget -nc https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
tar -zxvf grch38_genome.tar.gz

# 循环遍历所有SRA文件
for SRA_FILE in *.sra; do
    # 检测meta信息
    # echo "Checking SRA metadata for $SRA_FILE..."
    # sra-stat --meta $SRA_FILE

    # 转换SRA文件为fastq格式
    echo "Converting $SRA_FILE to fastq format..."
    fastq-dump --split-files $SRA_FILE

    # 获取不带扩展名的文件名
    BASENAME=$(basename $SRA_FILE .sra)

    # 质量控制
    echo "Performing quality control with FastQC..."
    fastqc ${BASENAME}_*.fastq

    # 修剪fastq文件，去除adapter
    echo "Trimming fastq files with fastp..."
    fastp -i ${BASENAME}_1.fastq -o ${BASENAME}_1_trimmed.fq

    # 映射到参考基因组
    echo "Mapping reads to the reference genome with HISAT2..."
    hisat2 -p 4 -x grch38/genome -1 ${BASENAME}_1_trimmed.fq -S ${BASENAME}_aligned.sam --dta

    # 将SAM转换为BAM并排序
    echo "Converting SAM to sorted BAM..."
    samtools sort -@ 4 -o ${BASENAME}_aligned.bam ${BASENAME}_aligned.sam

    # 计算基因表达量
    echo "Calculating gene expression with StringTie..."
    stringtie -p 4 -G Homo_sapiens.GRCh38.104.gtf -o ${BASENAME}.gtf -l ${BASENAME} -A ${BASENAME}_gene_abundances.tsv -B ${BASENAME}_aligned.bam

    # 生成VCF文件
    echo "Generating VCF file for $SRA_FILE..."
    bcftools mpileup -Ou -f Homo_sapiens.GRCh38.dna.primary_assembly.fa ${BASENAME}_aligned.bam -o ${BASENAME}.bcf
    bcftools call -Ov -mv -o ${BASENAME}.vcf ${BASENAME}.bcf
    bcftools filter -i 'QUAL > 30' ${BASENAME}.vcf -o ${BASENAME}_filtered.vcf

    # 将基因表达量汇总.tsv文件转换为.csv
    echo "Converting gene abundance TSV to CSV format..."
    cat ${BASENAME}_gene_abundances.tsv | tr '\t' ',' > ${BASENAME}_gene_abundances.csv

    echo "Analysis complete for $SRA_FILE."
done

echo "All analyses complete."
