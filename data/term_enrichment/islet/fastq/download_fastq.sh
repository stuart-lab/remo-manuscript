#!/usr/bin/env bash
#
# This script downloads HPAP data/metadata into your working directory:
#   * Data (if requested) will be saved in "hpapdata" subdirectory;
#   * Metadata (if requested) will be saved in "metadata" subdirectory.

DATA_SERVER="https://hpapdata.faryabilab.com"

FILES="
HPAP-175/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-175_ATAC_NSS0030_5000911_S5_L001_I1_001.fastq.gz
HPAP-175/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-175_ATAC_NSS0030_5000911_S5_L001_R1_001.fastq.gz
HPAP-175/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-175_ATAC_NSS0030_5000911_S5_L001_R2_001.fastq.gz
HPAP-175/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-175_ATAC_NSS0030_5000911_S5_L001_R3_001.fastq.gz
HPAP-175/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-175_ATAC_NSS0030_5000911_S5_L002_I1_001.fastq.gz
HPAP-175/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-175_ATAC_NSS0030_5000911_S5_L002_R1_001.fastq.gz
HPAP-175/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-175_ATAC_NSS0030_5000911_S5_L002_R2_001.fastq.gz
HPAP-175/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-175_ATAC_NSS0030_5000911_S5_L002_R3_001.fastq.gz
HPAP-175/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-175_GEX_NSS0032_5000912_S5_L001_I1_001.fastq.gz
HPAP-175/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-175_GEX_NSS0032_5000912_S5_L001_I2_001.fastq.gz
HPAP-175/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-175_GEX_NSS0032_5000912_S5_L001_R1_001.fastq.gz
HPAP-175/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-175_GEX_NSS0032_5000912_S5_L001_R2_001.fastq.gz
HPAP-175/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-175_GEX_NSS0032_5000912_S5_L002_I1_001.fastq.gz
HPAP-175/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-175_GEX_NSS0032_5000912_S5_L002_I2_001.fastq.gz
HPAP-175/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-175_GEX_NSS0032_5000912_S5_L002_R1_001.fastq.gz
HPAP-175/Islet Studies/Islet molecular phenotyping studies/Single-cell Multiome (ATAC+RNA)/Upenn_Multiome/HPAP-175_GEX_NSS0032_5000912_S5_L002_R2_001.fastq.gz
"

# Set IFS (Internal Field Separator)
IFS_BAK=$IFS
IFS=$'\n'

mkdir -p ./hpapdata
cd ./hpapdata

for f in $FILES; do
    echo "[$(date -Iseconds)] downloading $(basename $f) ..."
    encoded_f="$(echo $f | sed 's/ /%20/g')"
    curl --silent --create-dirs --output $f ${DATA_SERVER}/${encoded_f}
done

# Recover original IFS
IFS=$IFS_BAK
unset IFS_BAK

cd ..
echo; echo "[$(date -Iseconds)] experiment data downloaded"; echo

# Download metadata files into ./metadata
META_URL="https://hpap-public.s3.amazonaws.com/assets/metadata/latest"
mkdir -p ./metadata
cd ./metadata

curl --silent "$META_URL/PancDB_Donors.xlsx" -O
curl --silent "$META_URL/PancDB_snMultiome(ATAC+RNA)_Islets_metadata_2024-12-31.xlsx" -O
curl --silent "$META_URL/README.xlsx" -O

cd ..
echo "[$(date -Iseconds)] metadata downloaded"
