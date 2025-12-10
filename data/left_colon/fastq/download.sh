# download left colon
# https://www.encodeproject.org/multiomics-series/ENCSR925IHI/

# gene expression (SRA)
mkdir -p fastq_rna
cd fastq_rna
prefetch SRR32042721
prefetch SRR32042722
fastq-dump --split-files --gzip SRR32042721/SRR32042721.sra
fastq-dump --split-files --gzip SRR32042722/SRR32042722.sra

# rename fastq rna
mv SRR32042721_1.fastq.gz left_colon_S1_L001_R1_001.fastq.gz
mv SRR32042721_2.fastq.gz left_colon_S1_L001_R2_001.fastq.gz
mv SRR32042722_1.fastq.gz left_colon_S1_L002_R1_001.fastq.gz
mv SRR32042722_2.fastq.gz left_colon_S1_L002_R2_001.fastq.gz

cd .. 

# chromatin accessibility (fastq link)
mkdir -p fastq_atac
cd fastq_atac
wget https://www.encodeproject.org/files/ENCFF043NLM/@@download/ENCFF043NLM.fastq.gz
wget https://www.encodeproject.org/files/ENCFF980MGO/@@download/ENCFF980MGO.fastq.gz
wget https://www.encodeproject.org/files/ENCFF892TQK/@@download/ENCFF892TQK.fastq.gz

wget https://www.encodeproject.org/files/ENCFF555DAJ/@@download/ENCFF555DAJ.fastq.gz
wget https://www.encodeproject.org/files/ENCFF251XBS/@@download/ENCFF251XBS.fastq.gz
wget https://www.encodeproject.org/files/ENCFF424GPE/@@download/ENCFF424GPE.fastq.gz

# rename fastq atac
mv ENCFF043NLM.fastq.gz left_colon_S1_L001_R1_001.fastq.gz
mv ENCFF980MGO.fastq.gz left_colon_S1_L001_R3_001.fastq.gz
mv ENCFF892TQK.fastq.gz left_colon_S1_L001_R2_001.fastq.gz

mv ENCFF555DAJ.fastq.gz left_colon_S1_L002_R1_001.fastq.gz
mv ENCFF251XBS.fastq.gz left_colon_S1_L002_R3_001.fastq.gz
mv ENCFF424GPE.fastq.gz left_colon_S1_L002_R2_001.fastq.gz