# download muscle
# https://www.encodeproject.org/multiomics-series/ENCSR851GBP/

# gene expression (SRA)
mkdir -p fastq_rna
cd fastq_rna
prefetch SRR32042519
prefetch SRR32042520
fastq-dump --split-files --gzip SRR32042519/SRR32042519.sra
fastq-dump --split-files --gzip SRR32042520/SRR32042520.sra

# rename fastq rna
mv SRR32042519_1.fastq.gz muscle_S1_L001_R1_001.fastq.gz
mv SRR32042519_2.fastq.gz muscle_S1_L001_R2_001.fastq.gz
mv SRR32042520_1.fastq.gz muscle_S1_L002_R1_001.fastq.gz
mv SRR32042520_2.fastq.gz muscle_S1_L002_R2_001.fastq.gz
cd .. 

# chromatin accessibility (fastq link)
mkdir -p fastq_atac
cd fastq_atac
wget https://www.encodeproject.org/files/ENCFF747SKZ/@@download/ENCFF747SKZ.fastq.gz
wget https://www.encodeproject.org/files/ENCFF124DBI/@@download/ENCFF124DBI.fastq.gz
wget https://www.encodeproject.org/files/ENCFF734JYC/@@download/ENCFF734JYC.fastq.gz

wget https://www.encodeproject.org/files/ENCFF383ZFT/@@download/ENCFF383ZFT.fastq.gz
wget https://www.encodeproject.org/files/ENCFF723EHZ/@@download/ENCFF723EHZ.fastq.gz
wget https://www.encodeproject.org/files/ENCFF518BEG/@@download/ENCFF518BEG.fastq.gz

# rename fastq atac
mv ENCFF747SKZ.fastq.gz muscle_S1_L001_R1_001.fastq.gz
mv ENCFF124DBI.fastq.gz muscle_S1_L001_R3_001.fastq.gz
mv ENCFF734JYC.fastq.gz muscle_S1_L001_R2_001.fastq.gz

mv ENCFF383ZFT.fastq.gz muscle_S1_L002_R1_001.fastq.gz
mv ENCFF723EHZ.fastq.gz muscle_S1_L002_R3_001.fastq.gz
mv ENCFF518BEG.fastq.gz muscle_S1_L002_R2_001.fastq.gz