# download placenta
# https://www.encodeproject.org/multiomics-series/ENCSR694BTU/

# gene expression (SRA)
mkdir -p fastq_rna
cd fastq_rna
prefetch SRR32032790
prefetch SRR32032791
fastq-dump --split-files --gzip SRR32032790/SRR32032790.sra
fastq-dump --split-files --gzip SRR32032791/SRR32032791.sra

# rename fastq rna
mv SRR32032790_1.fastq.gz placenta_S1_L001_R1_001.fastq.gz
mv SRR32032790_2.fastq.gz placenta_S1_L001_R2_001.fastq.gz
mv SRR32032791_1.fastq.gz placenta_S1_L002_R1_001.fastq.gz
mv SRR32032791_2.fastq.gz placenta_S1_L002_R2_001.fastq.gz
cd .. 

# chromatin accessibility (fastq link)
mkdir -p fastq_atac
cd fastq_atac
wget https://www.encodeproject.org/files/ENCFF790RLP/@@download/ENCFF790RLP.fastq.gz
wget https://www.encodeproject.org/files/ENCFF358DNO/@@download/ENCFF358DNO.fastq.gz
wget https://www.encodeproject.org/files/ENCFF322UEB/@@download/ENCFF322UEB.fastq.gz

wget https://www.encodeproject.org/files/ENCFF927HGM/@@download/ENCFF927HGM.fastq.gz
wget https://www.encodeproject.org/files/ENCFF996ZIW/@@download/ENCFF996ZIW.fastq.gz
wget https://www.encodeproject.org/files/ENCFF411UNG/@@download/ENCFF411UNG.fastq.gz

# rename fastq atac
mv ENCFF790RLP.fastq.gz placenta_S1_L001_R1_001.fastq.gz
mv ENCFF358DNO.fastq.gz placenta_S1_L001_R3_001.fastq.gz
mv ENCFF322UEB.fastq.gz placenta_S1_L001_R2_001.fastq.gz

mv ENCFF927HGM.fastq.gz placenta_S1_L002_R1_001.fastq.gz
mv ENCFF996ZIW.fastq.gz placenta_S1_L002_R3_001.fastq.gz
mv ENCFF411UNG.fastq.gz placenta_S1_L002_R2_001.fastq.gz