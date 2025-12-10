# download fallopian tube
# https://www.encodeproject.org/multiomics-series/ENCSR420IUS/

# gene expression (SRA)
mkdir -p fastq_rna
cd fastq_rna
prefetch SRR32032620
prefetch SRR32032621
fastq-dump --split-files --gzip SRR32032620/SRR32032620.sra
fastq-dump --split-files --gzip SRR32032621/SRR32032621.sra

# rename fastq rna
mv SRR32032620_1.fastq.gz fallopian_tube_S1_L001_R1_001.fastq.gz
mv SRR32032620_2.fastq.gz fallopian_tube_S1_L001_R2_001.fastq.gz
mv SRR32032621_1.fastq.gz fallopian_tube_S1_L002_R1_001.fastq.gz
mv SRR32032621_2.fastq.gz fallopian_tube_S1_L002_R2_001.fastq.gz
cd .. 

# chromatin accessibility (fastq link)
mkdir -p fastq_atac
cd fastq_atac
wget https://www.encodeproject.org/files/ENCFF261UVP/@@download/ENCFF261UVP.fastq.gz
wget https://www.encodeproject.org/files/ENCFF257KOW/@@download/ENCFF257KOW.fastq.gz
wget https://www.encodeproject.org/files/ENCFF935EJL/@@download/ENCFF935EJL.fastq.gz

wget https://www.encodeproject.org/files/ENCFF764VSQ/@@download/ENCFF764VSQ.fastq.gz
wget https://www.encodeproject.org/files/ENCFF559TRX/@@download/ENCFF559TRX.fastq.gz
wget https://www.encodeproject.org/files/ENCFF869QLB/@@download/ENCFF869QLB.fastq.gz

# rename fastq atac
mv ENCFF261UVP.fastq.gz fallopian_tube_S1_L001_R1_001.fastq.gz
mv ENCFF257KOW.fastq.gz fallopian_tube_S1_L001_R3_001.fastq.gz
mv ENCFF935EJL.fastq.gz fallopian_tube_S1_L001_R2_001.fastq.gz

mv ENCFF764VSQ.fastq.gz fallopian_tube_S1_L002_R1_001.fastq.gz
mv ENCFF559TRX.fastq.gz fallopian_tube_S1_L002_R3_001.fastq.gz
mv ENCFF869QLB.fastq.gz fallopian_tube_S1_L002_R2_001.fastq.gz