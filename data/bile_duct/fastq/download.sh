# download bile duct
# https://www.encodeproject.org/multiomics-series/ENCSR871JTA/

# gene expression (SRA)
mkdir -p fastq_rna
cd fastq_rna
prefetch SRR32034006
prefetch SRR32034007
fastq-dump --split-files --gzip SRR32034006/SRR32034006.sra
fastq-dump --split-files --gzip SRR32034007/SRR32034007.sra

# rename fastq rna
mv SRR32034006_1.fastq.gz bile_duct_S1_L001_R1_001.fastq.gz
mv SRR32034006_2.fastq.gz bile_duct_S1_L001_R2_001.fastq.gz
mv SRR32034007_1.fastq.gz bile_duct_S1_L002_R1_001.fastq.gz
mv SRR32034007_2.fastq.gz bile_duct_S1_L002_R2_001.fastq.gz
cd .. 

# chromatin accessibility (fastq link)
mkdir -p fastq_atac
cd fastq_atac
wget https://www.encodeproject.org/files/ENCFF961PHG/@@download/ENCFF961PHG.fastq.gz
wget https://www.encodeproject.org/files/ENCFF096RMH/@@download/ENCFF096RMH.fastq.gz
wget https://www.encodeproject.org/files/ENCFF937TOU/@@download/ENCFF937TOU.fastq.gz

wget https://www.encodeproject.org/files/ENCFF669FSJ/@@download/ENCFF669FSJ.fastq.gz
wget https://www.encodeproject.org/files/ENCFF935JUA/@@download/ENCFF935JUA.fastq.gz
wget https://www.encodeproject.org/files/ENCFF420NYB/@@download/ENCFF420NYB.fastq.gz

# rename fastq atac
mv ENCFF961PHG.fastq.gz bile_duct_S1_L001_R1_001.fastq.gz
mv ENCFF096RMH.fastq.gz bile_duct_S1_L001_R3_001.fastq.gz
mv ENCFF937TOU.fastq.gz bile_duct_S1_L001_R2_001.fastq.gz

mv ENCFF669FSJ.fastq.gz bile_duct_S1_L002_R1_001.fastq.gz
mv ENCFF935JUA.fastq.gz bile_duct_S1_L002_R3_001.fastq.gz
mv ENCFF420NYB.fastq.gz bile_duct_S1_L002_R2_001.fastq.gz