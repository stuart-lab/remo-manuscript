# download fastq brain cancer I4
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM6257101

prefetch SRR19768587
prefetch SRR19768588
fasterq-dump --split-files --gzip SRR19768587
fasterq-dump --split-files --gzip SRR19768588

# I4   
mv SRR19768587_1.fastq.gz sampleI4_S1_L001_I1_001.fastq.gz
mv SRR19768587_2.fastq.gz sampleI4_S1_L001_R1_001.fastq.gz
mv SRR19768587_3.fastq.gz sampleI4_S1_L001_R2_001.fastq.gz
mv SRR19768587_4.fastq.gz sampleI4_S1_L001_R3_001.fastq.gz

mv SRR19768588_1.fastq.gz sampleI4_S1_L002_I1_001.fastq.gz
mv SRR19768588_2.fastq.gz sampleI4_S1_L002_R1_001.fastq.gz
mv SRR19768588_3.fastq.gz sampleI4_S1_L002_R2_001.fastq.gz
mv SRR19768588_4.fastq.gz sampleI4_S1_L002_R3_001.fastq.gz