# download fastq brain cancer M7
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6257104

prefetch SRR19768581
prefetch SRR19768582
fasterq-dump --split-files --gzip SRR19768581
fasterq-dump --split-files --gzip SRR19768582

# M7
mv SRR19768581_1.fastq.gz sampleM7_S1_L001_I1_001.fastq.gz
mv SRR19768581_2.fastq.gz sampleM7_S1_L001_R1_001.fastq.gz
mv SRR19768581_3.fastq.gz sampleM7_S1_L001_R2_001.fastq.gz
mv SRR19768581_4.fastq.gz sampleM7_S1_L001_R3_001.fastq.gz

mv SRR19768582_1.fastq.gz sampleM7_S1_L002_I1_001.fastq.gz
mv SRR19768582_2.fastq.gz sampleM7_S1_L002_R1_001.fastq.gz
mv SRR19768582_3.fastq.gz sampleM7_S1_L002_R2_001.fastq.gz
mv SRR19768582_4.fastq.gz sampleM7_S1_L002_R3_001.fastq.gz