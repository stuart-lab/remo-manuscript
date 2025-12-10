# cellranger commands for processing brain cancer FASTQ files
# software: cellranger-atac v2.1.0
# reference: refdata-cellranger-arc-GRCh38-2020-A-2.0.0

# note: --reference needs to be realpath

cellranger-atac count --id=sampleI2 \
  --reference=refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
  --fastqs=/data/disease/brain_cancer/I2/fastq \
  --sample=sampleI2 \
  --localcores=8 --localmem=128

cellranger-atac count --id=sampleI3 \
  --reference=refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
  --fastqs=/data/disease/brain_cancer/I3/fastq \
  --sample=sampleI3 \
  --localcores=8 --localmem=128

cellranger-atac count --id=sampleI4 \
  --reference=refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
  --fastqs=/data/disease/brain_cancer/I4/fastq \
  --sample=sampleI4 \
  --localcores=8 --localmem=128

cellranger-atac count --id=sampleI5 \
  --reference=refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
  --fastqs=/data/disease/brain_cancer/I5/fastq \
  --sample=sampleI5 \
  --localcores=8 --localmem=128

cellranger-atac count --id=sampleM7 \
  --reference=refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
  --fastqs=/data/disease/brain_cancer/M7/fastq \
  --sample=sampleM7 \
  --localcores=8 --localmem=128

cellranger-atac count --id=sampleM8 \
  --reference=refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
  --fastqs=/data/disease/brain_cancer/M8/fastq \
  --sample=sampleM8 \
  --localcores=8 --localmem=128