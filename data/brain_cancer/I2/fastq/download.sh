# download fastq brain cancer I2
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6257099

prefetch SRR19768591
prefetch SRR19768592
fasterq-dump --split-files --gzip SRR19768591
fasterq-dump --split-files --gzip SRR19768592

# I2
mv SRR19768591_1.fastq.gz sampleI2_S1_L001_I1_001.fastq.gz
mv SRR19768591_2.fastq.gz sampleI2_S1_L001_R1_001.fastq.gz
mv SRR19768591_3.fastq.gz sampleI2_S1_L001_R2_001.fastq.gz
mv SRR19768591_4.fastq.gz sampleI2_S1_L001_R3_001.fastq.gz

mv SRR19768592_1.fastq.gz sampleI2_S1_L002_I1_001.fastq.gz
mv SRR19768592_2.fastq.gz sampleI2_S1_L002_R1_001.fastq.gz
mv SRR19768592_3.fastq.gz sampleI2_S1_L002_R2_001.fastq.gz
mv SRR19768592_4.fastq.gz sampleI2_S1_L002_R3_001.fastq.gz