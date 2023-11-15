# 使用 Ubuntu 作为基础镜像
FROM ubuntu:latest

# 安装必要的依赖
RUN apt-get update && apt-get install -y wget tar unzip make gcc g++ git cmake python3 python3-pip openjdk-17-jdk

# 安装samtools的依赖
RUN apt-get install -y autoconf automake make gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev vim
RUN apt-get update && apt-get install -y libncurses5-dev


# 添加新用户 'ubuntu'
RUN useradd -m ubuntu

RUN apt-get update && apt-get install -y sudo curl
RUN adduser ubuntu sudo
RUN echo "ubuntu ALL=(ALL) NOPASSWD: ALL" | sudo tee /etc/sudoers.d/ubuntu

USER ubuntu

# 设置工作目录
WORKDIR /home/ubuntu

SHELL ["/bin/bash", "-c"]

ENV RNA_HOME=/home/ubuntu/workspace/rnaseq \
    RNA_EXT_DATA_DIR=/home/ubuntu/CourseData/RNA_data \
    RNA_DATA_DIR=/home/ubuntu/workspace/rnaseq/data \
    RNA_DATA_TRIM_DIR=/home/ubuntu/workspace/rnaseq/data/trimmed \
    RNA_REFS_DIR=/home/ubuntu/workspace/rnaseq/refs \
    RNA_REF_INDEX=/home/ubuntu/workspace/rnaseq/refs/chr22_with_ERCC92 \
    RNA_REF_FASTA=/home/ubuntu/workspace/rnaseq/refs/chr22_with_ERCC92.fa \
    RNA_REF_GTF=/home/ubuntu/workspace/rnaseq/refs/chr22_with_ERCC92.gtf \
    RNA_ALIGN_DIR=/home/ubuntu/workspace/rnaseq/alignments/hisat2 \
    WORKSPACE=/home/ubuntu/workspace \
    HOME=/home/ubuntu \
    PICARD=/home/ubuntu/bin/picard.jar \
    GATK_REGIONS='-L chr17'

# 设置环境变量
RUN cd ~ && \
    wget http://genomedata.org/rnaseq-tutorial/bashrc_copy && \
    mv bashrc_copy ~/.bashrc && \
    echo "export PATH=/home/ubuntu/bin/tophat-2.1.1.Linux_x86_64:\$PATH" >> ~/.bashrc && \
    source ~/.bashrc && \
    mkdir -p ~/workspace/rnaseq/

# 设置目录
RUN mkdir -p $RNA_HOME/student_tools

# 安装samtools
RUN mkdir -p ~/bin/ && \
    cd ~/bin/ && \
    wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 && \
    bunzip2 samtools-1.18.tar.bz2 && \
    tar -xvf samtools-1.18.tar && \
    cd samtools-1.18 && \
    ./configure && \
    make

# 安装bam-readcount
RUN cd ~/bin/ && \
    git clone https://github.com/genome/bam-readcount  && \
    cd bam-readcount && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make

# 安装HISAT2
RUN cd ~/bin/ && \
    wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download && \
    mv download hisat2-2.2.1-Linux_x86_64.zip && \
    unzip hisat2-2.2.1-Linux_x86_64.zip && \
    rm hisat2-2.2.1-Linux_x86_64.zip
# 安装StringTie
RUN cd ~/bin/ && \
    wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.2.1.Linux_x86_64.tar.gz && \
    tar -xzvf stringtie-2.2.1.Linux_x86_64.tar.gz && \
    mv stringtie-2.2.1.Linux_x86_64 stringtie-2.2.1
    
# 安装gffcompare
RUN cd ~/bin/ && \
    wget http://ccb.jhu.edu/software/stringtie/dl/gffcompare-0.12.6.Linux_x86_64.tar.gz && \
    tar -xzvf gffcompare-0.12.6.Linux_x86_64.tar.gz

# 安装htseq-count
RUN pip install HTSeq

# 安装TopHat
RUN cd ~/bin/ && \
    wget http://genomedata.org/rnaseq-tutorial/tophat-2.1.1.Linux_x86_64.tar.gz && \
    tar -zxvf tophat-2.1.1.Linux_x86_64.tar.gz

# 安装kallisto
RUN cd ~/bin/ && \
    wget https://github.com/pachterlab/kallisto/releases/download/v0.44.0/kallisto_linux-v0.44.0.tar.gz && \
    tar -zxvf kallisto_linux-v0.44.0.tar.gz

# 安装FastQC
RUN cd ~/bin/ && \
    wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip && \
    unzip fastqc_v0.12.1.zip

# 安装Fastp
RUN cd ~/bin/ && \
    mkdir fastp && \
    cd fastp && \
    wget http://opengene.org/fastp/fastp && \
    chmod a+x ./fastp

# 安装MultiQC
RUN pip install multiqc

# 安装Picard
RUN cd ~/bin/ && \
    wget https://github.com/broadinstitute/picard/releases/download/3.1.0/picard.jar

# 安装flexbar
RUN cd ~/bin/ && \
    wget https://github.com/seqan/flexbar/releases/download/v3.5.0/flexbar-3.5.0-linux.tar.gz && \
    tar -zxvf flexbar-3.5.0-linux.tar.gz && \
    echo "LD_LIBRARY_PATH=/home/ubuntu/bin/flexbar-3.5.0-linux:\$LD_LIBRARY_PATH" >> ~/.bashrc && \
    source ~/.bashrc

# 安装RegTools
RUN cd ~/bin/ && \
    git clone https://github.com/griffithlab/regtools  && \
    cd regtools/ && \
    mkdir build && \
    cd build/ && \
    cmake .. && \
    make

# 安装RSeQC
RUN pip install RSeQC

# 安装bedops
RUN cd ~/bin/ && \
    mkdir bedops_linux_x86_64-v2.4.41 && \
    cd bedops_linux_x86_64-v2.4.41 && \
    wget https://github.com/bedops/bedops/releases/download/v2.4.41/bedops_linux_x86_64-v2.4.41.tar.bz2 && \
    tar -jxvf bedops_linux_x86_64-v2.4.41.tar.bz2

# 安装gtfToGenePred
RUN cd ~/bin/ && \
    mkdir gtfToGenePred && \
    cd gtfToGenePred && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred && \
    chmod a+x gtfToGenePred

# 安装genePredToBed
RUN cd ~/bin/ && \
    mkdir genePredToBed && \
    cd genePredToBed && \
    wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToBed && \
    chmod a+x genePredToBed

# 安装how_are_we_stranded_here
RUN pip install how_are_we_stranded_here

# 安装CellRanger
RUN cd ~/bin/ && \
    wget https://alist.ivanq.top/d/whu-share/cellranger-7.2.0.tar.gz
RUN cd ~/bin/ && \
    tar -xzf cellranger-7.2.0.tar.gz


USER root
ARG TZ=Asia/Shanghai
ARG DEBIAN_FRONTEND=noninteractive
# 安装R
# 安装需要的软件
RUN apt-get install --no-install-recommends -y software-properties-common dirmngr && \
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc && \
    add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

# 安装 R
RUN apt-get update -qq && apt-get install -y --no-install-recommends r-base

# 为了devtools以及edgeR的安装，需要安装一些依赖
RUN apt-get install -y libxml2-dev build-essential libharfbuzz-dev libfribidi-dev libcairo2-dev libmagick++-dev gfortran libeigen3-dev liblapack-dev libblas-dev


# 设置 R 环境以非交互式方式安装包，使用 TUNA 镜像源
RUN Rscript -e "options(repos = list(CRAN = 'https://mirrors.tuna.tsinghua.edu.cn/CRAN/'));install.packages(c('devtools'),dependencies = TRUE)"

RUN Rscript -e "options(repos = list(CRAN = 'https://mirrors.tuna.tsinghua.edu.cn/CRAN/'));install.packages(c('dplyr', 'gplots', 'tidyverse'))"

# 安装 Bioconductor 包，使用 TUNA 镜像源
RUN Rscript -e 'options(repos = list(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")); if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); BiocManager::install(c("genefilter"),dependencies = TRUE,force = TRUE)'
RUN Rscript -e 'options(repos = list(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")); BiocManager::install(c("edgeR","ballgown"),dependencies = TRUE)'
RUN Rscript -e 'options(repos = list(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")); BiocManager::install(c("GenomicRanges", "rhdf5", "biomaRt"),dependencies = TRUE,force = TRUE)'
# 安装 devtools 并通过 devtools 安装 sleuth
RUN Rscript -e "options(repos = list(CRAN = 'https://mirrors.tuna.tsinghua.edu.cn/CRAN/')); BiocManager::install('pachterlab/sleuth')"

ARG DEBIAN_FRONTEND=dialog
USER ubuntu

# 设置默认命令
CMD ["/bin/bash"]