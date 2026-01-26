#!/bin/bash

#SBATCH -J linkpeaks
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=10G
#SBATCH --output=./log/%x-%j.out
#SBATCH --error=./log/%x-%j.err
#SBATCH -p interactive

# ------------------------------------------------------------------------------
PROJECT_ROOT="/Projects/linkage_analysis" # "path/to/renv_proj"
OUTPUT_DIR="${PROJECT_ROOT}/data/multiome_links/jejunum" # "/data/multiome_links/dataset_name"
SEURAT_RDS="${PROJECT_ROOT}/resources/multiome/jejunum_multiome_linkpeaks.rds" # "/resources/multiome/your_data.rds"
# ------------------------------------------------------------------------------

R_SCRIPT="${PROJECT_ROOT}/scripts/01_link_peaks.R"

REMO_BED="${PROJECT_ROOT}/resources/REMOv1_hg38.bed"
GENCODE_GTF="${PROJECT_ROOT}/resources/gencode.v32.basic.annotation.gtf.gz"

WORKERS="${SLURM_CPUS_PER_TASK:-1}"

## Activate env with correct vers of renv installed in HPC
# module load mambaforge
# conda deactivate
# mamba activate mamba_env

cd "${PROJECT_ROOT}" || exit 1   # <-- ensures renv auto-loads via .Rprofile
test -f renv/activate.R || { echo "Missing renv/activate.R in $(pwd)"; exit 2; }

cmd=(
  Rscript "${R_SCRIPT}"
  --seurat_rds "${SEURAT_RDS}"
  --remo_bed "${REMO_BED}"
  --gencode_gtf "${GENCODE_GTF}"
  --output_dir "${OUTPUT_DIR}"
  --workers "${WORKERS}"
)

echo "Running command:"
printf '  %q ' "${cmd[@]}"
echo

"${cmd[@]}"
