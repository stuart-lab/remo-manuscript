# download fastq brain cancer I3
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6257100

prefetch SRR19768589
prefetch SRR19768590
fasterq-dump --split-files --gzip SRR19768589
fasterq-dump --split-files --gzip SRR19768590

# I3  
mv SRR19768589_1.fastq.gz sampleI3_S1_L001_I1_001.fastq.gz
mv SRR19768589_2.fastq.gz sampleI3_S1_L001_R1_001.fastq.gz
mv SRR19768589_3.fastq.gz sampleI3_S1_L001_R2_001.fastq.gz
mv SRR19768589_4.fastq.gz sampleI3_S1_L001_R3_001.fastq.gz

mv SRR19768590_1.fastq.gz sampleI3_S1_L002_I1_001.fastq.gz
mv SRR19768590_2.fastq.gz sampleI3_S1_L002_R1_001.fastq.gz
mv SRR19768590_3.fastq.gz sampleI3_S1_L002_R2_001.fastq.gz
mv SRR19768590_4.fastq.gz sampleI3_S1_L002_R3_001.fastq.gz