# download liver
# https://www.encodeproject.org/multiomics-series/ENCSR728OVE/

# gene expression (SRA)
mkdir -p fastq_rna
cd fastq_rna
prefetch SRR32037199
prefetch SRR32037200
fastq-dump --split-files --gzip SRR32037199/SRR32037199.sra
fastq-dump --split-files --gzip SRR32037200/SRR32037200.sra

# rename fastq rna
mv SRR32037199_1.fastq.gz liver_S1_L001_R1_001.fastq.gz
mv SRR32037199_2.fastq.gz liver_S1_L001_R2_001.fastq.gz
mv SRR32037200_1.fastq.gz liver_S1_L002_R1_001.fastq.gz
mv SRR32037200_2.fastq.gz liver_S1_L002_R2_001.fastq.gz
cd .. 

# chromatin accessibility (fastq link)
mkdir -p fastq_atac
cd fastq_atac
wget https://www.encodeproject.org/files/ENCFF414GGW/@@download/ENCFF414GGW.fastq.gz
wget https://www.encodeproject.org/files/ENCFF290TFC/@@download/ENCFF290TFC.fastq.gz
wget https://www.encodeproject.org/files/ENCFF714FWE/@@download/ENCFF714FWE.fastq.gz

wget https://www.encodeproject.org/files/ENCFF451MFE/@@download/ENCFF451MFE.fastq.gz
wget https://www.encodeproject.org/files/ENCFF777ORW/@@download/ENCFF777ORW.fastq.gz
wget https://www.encodeproject.org/files/ENCFF969GSH/@@download/ENCFF969GSH.fastq.gz

wget https://www.encodeproject.org/files/ENCFF665CUM/@@download/ENCFF665CUM.fastq.gz
wget https://www.encodeproject.org/files/ENCFF657GDU/@@download/ENCFF657GDU.fastq.gz
wget https://www.encodeproject.org/files/ENCFF055GDP/@@download/ENCFF055GDP.fastq.gz

# rename fastq atac
mv ENCFF414GGW.fastq.gz liver_S1_L001_R1_001.fastq.gz
mv ENCFF290TFC.fastq.gz liver_S1_L001_R3_001.fastq.gz
mv ENCFF714FWE.fastq.gz liver_S1_L001_R2_001.fastq.gz

mv ENCFF451MFE.fastq.gz liver_S1_L002_R1_001.fastq.gz
mv ENCFF777ORW.fastq.gz liver_S1_L002_R3_001.fastq.gz
mv ENCFF969GSH.fastq.gz liver_S1_L002_R2_001.fastq.gz

mv ENCFF665CUM.fastq.gz liver_S1_L003_R1_001.fastq.gz
mv ENCFF657GDU.fastq.gz liver_S1_L003_R3_001.fastq.gz
mv ENCFF055GDP.fastq.gz liver_S1_L003_R2_001.fastq.gz