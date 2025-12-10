# download lung
# https://www.encodeproject.org/multiomics-series/ENCSR264JIX/

# gene expression (SRA)
mkdir -p fastq_rna
cd fastq_rna
prefetch SRR32072050
prefetch SRR32072051
fastq-dump --split-files --gzip SRR32072050/SRR32072050.sra
fastq-dump --split-files --gzip SRR32072051/SRR32072051.sra

# rename fastq rna
mv SRR32072050_1.fastq.gz lung_S1_L001_R1_001.fastq.gz
mv SRR32072050_2.fastq.gz lung_S1_L001_R2_001.fastq.gz
mv SRR32072051_1.fastq.gz lung_S1_L002_R1_001.fastq.gz
mv SRR32072051_2.fastq.gz lung_S1_L002_R2_001.fastq.gz
cd .. 

# chromatin accessibility (fastq link)
mkdir -p fastq_atac
cd fastq_atac
wget https://www.encodeproject.org/files/ENCFF373DBG/@@download/ENCFF373DBG.fastq.gz
wget https://www.encodeproject.org/files/ENCFF223EZM/@@download/ENCFF223EZM.fastq.gz
wget https://www.encodeproject.org/files/ENCFF057YPM/@@download/ENCFF057YPM.fastq.gz

wget https://www.encodeproject.org/files/ENCFF605ZBO/@@download/ENCFF605ZBO.fastq.gz
wget https://www.encodeproject.org/files/ENCFF990BFG/@@download/ENCFF990BFG.fastq.gz
wget https://www.encodeproject.org/files/ENCFF113VDA/@@download/ENCFF113VDA.fastq.gz

# rename fastq atac
mv ENCFF373DBG.fastq.gz lung_S1_L001_R1_001.fastq.gz
mv ENCFF223EZM.fastq.gz lung_S1_L001_R3_001.fastq.gz
mv ENCFF057YPM.fastq.gz lung_S1_L001_R2_001.fastq.gz

mv ENCFF605ZBO.fastq.gz lung_S1_L002_R1_001.fastq.gz
mv ENCFF990BFG.fastq.gz lung_S1_L002_R3_001.fastq.gz
mv ENCFF113VDA.fastq.gz lung_S1_L002_R2_001.fastq.gz