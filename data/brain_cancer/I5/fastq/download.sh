# download fastq brain cancer I5
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6257102

prefetch SRR19768585
prefetch SRR19768586
fasterq-dump --split-files --gzip SRR19768585
fasterq-dump --split-files --gzip SRR19768586

# I5 
mv SRR19768585_1.fastq.gz sampleI5_S1_L001_I1_001.fastq.gz
mv SRR19768585_2.fastq.gz sampleI5_S1_L001_R1_001.fastq.gz
mv SRR19768585_3.fastq.gz sampleI5_S1_L001_R2_001.fastq.gz
mv SRR19768585_4.fastq.gz sampleI5_S1_L001_R3_001.fastq.gz

mv SRR19768586_1.fastq.gz sampleI5_S1_L002_I1_001.fastq.gz
mv SRR19768586_2.fastq.gz sampleI5_S1_L002_R1_001.fastq.gz
mv SRR19768586_3.fastq.gz sampleI5_S1_L002_R2_001.fastq.gz
mv SRR19768586_4.fastq.gz sampleI5_S1_L002_R3_001.fastq.gz