# download fastq brain cancer M8
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6257103

prefetch SRR19768583
prefetch SRR19768584
fasterq-dump --split-files --gzip SRR19768583
fasterq-dump --split-files --gzip SRR19768584

# M8
mv SRR19768583_1.fastq.gz sampleM8_S1_L001_I1_001.fastq.gz
mv SRR19768583_2.fastq.gz sampleM8_S1_L001_R1_001.fastq.gz
mv SRR19768583_3.fastq.gz sampleM8_S1_L001_R2_001.fastq.gz
mv SRR19768583_4.fastq.gz sampleM8_S1_L001_R3_001.fastq.gz

mv SRR19768584_1.fastq.gz sampleM8_S1_L002_I1_001.fastq.gz
mv SRR19768584_2.fastq.gz sampleM8_S1_L002_R1_001.fastq.gz
mv SRR19768584_3.fastq.gz sampleM8_S1_L002_R2_001.fastq.gz
mv SRR19768584_4.fastq.gz sampleM8_S1_L002_R3_001.fastq.gz