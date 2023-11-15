# SRA文件
（以下均在docker容器中进行，建议将本地文件夹挂载到docker中。以下步骤都在同一个文件夹发生。）
我们从获取的SRA文件开始，这是一种NCBI的压缩格式，需要使用特定的工具将其转换为更常见的如fastq格式。工具是sra-toolkit

```
# 更新包，Ubuntu
sudo apt-get update 
# 下载sra-toolkit
sudo apt-get install libxml2 libxml2-dev libxslt1-dev zlib1g-dev ncurses-dev

wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz

tar -vxzf sratoolkit.current-ubuntu64.tar.gz

# 添加环境变量，到~/.bashrc，改版本号
export PATH=$PATH:/path/to/sratoolkit.3.0.7-ubuntu64/bin

source ~/.bashrc

# 验证安装
fastq-dump --version
```
## 格式转换
cd到sra的文件夹
```
# 检测meta信息，看看是不是paired-end
sra-stat --meta SRR7961253.sra

fastq-dump --split-files SRR7961253.sra

```

## 质量控制，看看这个数据集的质量
```
fastqc *.fastq

# 单端测序
fastp -i path.fastq 

# 双端测序
fastp -i path1.fastq -I path2.fastq
```

# 参考基因组
有了测序文件后，我们需要参考基因组以及基因注释
本例中，我们需要人类基因组的所有染色体
从ensembl库下载GRCh38的整体fasta文件，大概840MB，以及完整的注释文件，约47MB。注意release可能已经过时，建议先搜索最新版，更改命令。
```
wget ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

wget ftp://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz

# 预览基因组文件
zcat Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz | head -1000 | tail -30
# 预览注释文件
zcat Homo_sapiens.GRCh38.104.gtf.gz |head -1000 |tail -30

zcat Homo_sapiens.GRCh38.104.gtf.gz |head -100 |column -t | less -p exon -S

```
看到第一个基因“ENSG00000284662”

```
zgrep NSG00000284662 Homo_sapiens.GRCh38.104.gtf.gz | column -t |less -p "exon\s" -S
```
需要稍等一会才会出结果。

## 对参考基因组建立索引
为了加快测序文件到参考基因组的映射过程，需要对参考基因组进行索引。可以在本机建立索引，但需要160GB的RAM。但这种大项目早有人做好了索引，我们只需要下载即可，因为之后我们使用HISAT2工具，因此使用HISAT的索引文件。约3.9GB
```
wget https://genome-idx.s3.amazonaws.com/hisat/grch38_genome.tar.gz
# 解压
tar -zxvf grch38_genome.tar.gz
```

## 修剪fastq文件，去除adapter
一般用fastp进行这一步。
```
# -o 表示输出名
fastp -i SRR7961253_1.fastq -o SRR7961253_1_trimmed.fq
```


# 映射
## 建立参考转录组
以下解压过程如果已经完成可以不做。但要确保命名一致。
```
# 解压gtf文件
gunzip -c Homo_sapiens.GRCh38.104.gtf.gz > Homo_sapiens.GRCh38.104.gtf
# 转为genepred文件
gtfToGenePred Homo_sapiens.GRCh38.104.gtf Homo_sapiens.GRCh38.104.genePred
# 转为bed12文件
genePredToBed Homo_sapiens.GRCh38.104.genePred Homo_sapiens.GRCh38.104.bed12
# 解压参考基因组文件
gunzip -c Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > Homo_sapiens.GRCh38.dna.primary_assembly.fa
# 建立参考转录组，输入是参考基因组文件，bed12文件，输出转录组
bedtools getfasta -fi Homo_sapiens.GRCh38.dna.primary_assembly.fa -bed Homo_sapiens.GRCh38.104.bed12 -s -split -name -fo Homo_sapiens.GRCh38.104_transcripts.fa
```
！！！
然而，这种事情当然也有别人做啦，我们只需要下载就好了。注意release版本跟之前一样。
```
wget ftp://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
# 解压
gunzip -c Homo_sapiens.GRCh38.cdna.all.fa.gz > Homo_sapiens.GRCh38.104_transcripts.fa

```

## 建立链特异性（为映射做准备）
如果之前都在docker内运行，那么需要在外面的机子运行以下docker命令。或者搞docker in docker
输入是gtf文件，转录组文件，以及测序文件
```
# 对于双端
docker run -v /YOUR/PATH/TO/RNA-work:/docker_workspace mgibio/checkstrandedness:latest check_strandedness --gtf /docker_workspace/Homo_sapiens.GRCh38.104.gtf --transcripts /docker_workspace/Homo_sapiens.GRCh38.104_transcripts.fa --reads_1 SRR7961253_1.fastq --reads_2 .....
```
但由于本例是单端，因此无法使用该工具。

## 使用HISAT2进行映射
根据前面建立的索引，将测序文件映射到参考基因组上。由于我们后续需要计算表达量因此加上--dta
```
# -x 后是索引文件夹，刚解压的，-U是输入fastq文件，刚修剪的，-S是输出名
hisat2 -p 4 -x grch38/genome -U SRR7961253_1_trimmed.fq -S SRR7961253_aligned.sam --dta

```
若知道链特异性，则加--rna-strandness选项。

# 计算基因表达量
刚刚使用的 hisat2 命令是用于比对到人类基因组的，并且在命令中使用了 --dta 选项，这是为了生成用于转录本组装和表达量估计的数据。这表明我们对转录本级别的分析感兴趣，而不仅仅是基因级别。基于这一点，可以使用 StringTie 来计算表达量。
## 将sam文件变为bam让程序更好处理，以及利于下一步Stringtie

```
# 变得过程中同时排序
samtools sort -@ 4 -o SRR7961253_aligned.bam SRR7961253_aligned.sam
```

## Stringtie
```
stringtie -p 4 -G Homo_sapiens.GRCh38.104.gtf -o SRR7961253.gtf -l SRR7961253 -A SRR7961253_gene_abundances.tsv -B SRR7961253_aligned.bam
```

输出一个转录本的组装gtf文件，以及转录本表达量的数据，在tsv中。
-B表示输出ctab文件，为ballgown做准备。包括e_data.ctab（外显子）、i_data.ctab（内含子）、t_data.ctab（转录本）

## 有了Stringtie的输出就可以用各种工具进行后续分析了

我们尝试直接在IGV中可视化基因表达量。

```
gtf2bed --gtf SRR7961253.gtf --bed SRR7961253.bed
```