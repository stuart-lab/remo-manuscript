# cellranger commands for processing FASTQ files
# software: cellranger-ARC v2.0.2
# reference: refdata-cellranger-arc-GRCh38-2020-A-2.0.0

# note: --reference & --libraries needs to be realpath

## heart
cellranger-arc count --id=heart \
            --reference=refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
            --libraries=/data/heart/fastq/libraries_heart.csv \
            --localcores=16 \
            --localmem=128

## left_colon
cellranger-arc count --id=left_colon \
            --reference=refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
            --libraries=/data/left_colon/fastq/libraries_left_colon.csv \
            --localcores=16 \
            --localmem=128

## liver
cellranger-arc count --id=liver \
            --reference=refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
            --libraries=/data/liver/fastq/libraries_liver.csv \
            --localcores=16 \
            --localmem=128

## lung
cellranger-arc count --id=lung \
            --reference=refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
            --libraries=/data/lung/fastq/libraries_lung.csv \
            --localcores=16 \
            --localmem=128

## muscle
cellranger-arc count --id=muscle \
            --reference=refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
            --libraries=/data/muscle/fastq/libraries_muscle.csv \
            --localcores=16 \
            --localmem=128

## pancreas
cellranger-arc count --id=pancreas \
            --reference=refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
            --libraries=/data/pancreas/fastq/libraries_pancreas.csv \
            --localcores=16 \
            --localmem=128

## bile duct
cellranger-arc count --id=bile_duct \
            --reference=refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
            --libraries=/data/bile_duct/fastq/libraries_bile_duct.csv \
            --localcores=16 \
            --localmem=128

## fallopian tube
cellranger-arc count --id=fallopian_tube \
            --reference=refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
            --libraries=/data/fallopian_tube/fastq/libraries_fallopian_tube.csv \
            --localcores=16 \
            --localmem=128

## placenta
cellranger-arc count --id=placenta \
            --reference=refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
            --libraries=/data/placenta/fastq/libraries_placenta.csv \
            --localcores=16 \
            --localmem=128

## islet
cellranger-arc count --id=HPAP-175 \
            --reference=refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
            --libraries=/data/term_enrichment/islet/fastq/libraries_islet.csv \
            --localcores=16 \
            --localmem=128