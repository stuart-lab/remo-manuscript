#!/bin/bash

#SBATCH -J linkage_meta
#SBATCH -t 1-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --output=./log/%x-%j.out
#SBATCH --error=./log/%x-%j.err
#SBATCH -p cpu

# ------------------------------------------------------------------------------
PROJECT_ROOT="/Projects/linkage_analysis" # "path/to/renv_proj"

LINKS_DIR_CSV="${PROJECT_ROOT}/scripts/config/objects_and_links.csv" # modify example table of seurat obj fpaths and their mutiome_link dirs
FILTER_CONDS_CSV="${PROJECT_ROOT}/scripts/config/link_filtering_conds.csv" # specify pvals and scores for links to keep (multiple conds possible)
# ------------------------------------------------------------------------------

R_SCRIPT="${PROJECT_ROOT}/scripts/02_aggregate_links.R"

REMO_BED="${PROJECT_ROOT}/resources/REMOv1_GRCh38.bed.gz"
GENCODE_GTF="${PROJECT_ROOT}/resources/gencode.v32.basic.annotation.gtf.gz"

COMBINED_OUTDIR="${PROJECT_ROOT}/data/combined_metadata"
FILTERED_PARENT_DIR="${PROJECT_ROOT}/data/filtered_links"

DISTANCE_TO_TSS="500"
OVERWRITE="1"   # set to 1 to pass --overwrite

## Activate env with correct vers of renv installed in HPC
# module load mambaforge
# conda deactivate
# mamba activate mamba_env

mkdir -p ./log
cd "${PROJECT_ROOT}" || exit 1  # ensures renv auto-loads via .Rprofile, plus paths are consistent

cmd=(
  Rscript "${R_SCRIPT}"
  --links_dir_csv "${LINKS_DIR_CSV}"
  --filter_conds_csv "${FILTER_CONDS_CSV}"
  --combined_outdir "${COMBINED_OUTDIR}"
  --filtered_parent_dir "${FILTERED_PARENT_DIR}"
  --gencode_gtf "${GENCODE_GTF}"
  --remo_bed "${REMO_BED}"
  --distance_to_tss "${DISTANCE_TO_TSS}"
)

if [[ "${OVERWRITE}" == "1" ]]; then
  cmd+=(--overwrite)
fi

echo "Running command:"
printf '  %q ' "${cmd[@]}"
echo

"${cmd[@]}"
