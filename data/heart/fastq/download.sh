# download heart 
# https://www.encodeproject.org/multiomics-series/ENCSR302EOG/

# gene expression (SRA)
mkdir -p fastq_rna
cd fastq_rna
prefetch SRR32048549
prefetch SRR32048550
fastq-dump --split-files --gzip SRR32048549/SRR32048549.sra
fastq-dump --split-files --gzip SRR32048550/SRR32048550.sra

# rename fastq rna
mv SRR32048549_1.fastq.gz heart_S1_L001_R1_001.fastq.gz
mv SRR32048549_2.fastq.gz heart_S1_L001_R2_001.fastq.gz
mv SRR32048550_1.fastq.gz heart_S1_L002_R1_001.fastq.gz
mv SRR32048550_2.fastq.gz heart_S1_L002_R2_001.fastq.gz

cd .. 

# chromatin accessibility (fastq link)
mkdir -p fastq_atac
cd fastq_atac
wget https://www.encodeproject.org/files/ENCFF352YRT/@@download/ENCFF352YRT.fastq.gz
wget https://www.encodeproject.org/files/ENCFF018QRL/@@download/ENCFF018QRL.fastq.gz
wget https://www.encodeproject.org/files/ENCFF003FBC/@@download/ENCFF003FBC.fastq.gz

wget https://www.encodeproject.org/files/ENCFF871WLY/@@download/ENCFF871WLY.fastq.gz
wget https://www.encodeproject.org/files/ENCFF468VAM/@@download/ENCFF468VAM.fastq.gz
wget https://www.encodeproject.org/files/ENCFF692YAN/@@download/ENCFF692YAN.fastq.gz

# rename fastq atac
mv ENCFF352YRT.fastq.gz heart_S1_L001_R1_001.fastq.gz
mv ENCFF018QRL.fastq.gz heart_S1_L001_R3_001.fastq.gz
mv ENCFF003FBC.fastq.gz heart_S1_L001_R2_001.fastq.gz

mv ENCFF871WLY.fastq.gz heart_S1_L002_R1_001.fastq.gz
mv ENCFF468VAM.fastq.gz heart_S1_L002_R3_001.fastq.gz
mv ENCFF692YAN.fastq.gz heart_S1_L002_R2_001.fastq.gz