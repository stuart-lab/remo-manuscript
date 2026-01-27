#!/bin/bash

#SBATCH -J prob_analysis
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

# use the combined pgl_metadata for the condition of choice (e.g. sig_pvalue)
# to reproduce analysis use example file in resources (see README)
SIG_PGL_FPATH="${PROJECT_ROOT}/data/filtered_links/sig_pvalue/pgl_metadata/pgls_metadata_all_chr.tsv"  
# ------------------------------------------------------------------------------

SCRIPT_COND_PROBS="${PROJECT_ROOT}/scripts/03_get_conditional_probs.R"
SCRIPT_GLM_COEFFS="${PROJECT_ROOT}/scripts/03a_get_glm_coeffs.R"
SCRIPT_GLM_SUMMARY="${PROJECT_ROOT}/scripts/03b_build_glm_summary.R"

# Shared output dir
OUTPUT_DIR="${PROJECT_ROOT}/data/probability_analysis"
mkdir -p "${OUTPUT_DIR}"

# Params for 03b_build_glm_summary.R
GLM_DIR="${OUTPUT_DIR}/logistic_regression_results"
OUTCOME="same_remo"
TERM="same_gene"

## Activate env with correct vers of renv installed in HPC
# module load mambaforge
# conda deactivate
# mamba activate mamba_env

mkdir -p ./log
cd "${PROJECT_ROOT}" || exit 1

echo "Running 03_get_conditional_probs.R"
cmd1=(
	Rscript "${SCRIPT_COND_PROBS}"
	--renv_dir "${PROJECT_ROOT}"
	--sig_pgl_fpath "${SIG_PGL_FPATH}"
	--output_dir "${OUTPUT_DIR}"
)
printf '  %q ' "${cmd1[@]}"; echo
"${cmd1[@]}"
echo

echo "Running 03a_get_glm_coeffs.R"
cmd2=(
	Rscript "${SCRIPT_GLM_COEFFS}"
	--renv_dir "${PROJECT_ROOT}"
	--sig_pgl_fpath "${SIG_PGL_FPATH}"
	--output_dir "${OUTPUT_DIR}"
)
printf '  %q ' "${cmd2[@]}"; echo
"${cmd2[@]}"
echo

echo "Running 03b_build_glm_summary.R"
cmd3=(
	Rscript "${SCRIPT_GLM_SUMMARY}"
	--renv_dir "${PROJECT_ROOT}"
	--glm_dir "${GLM_DIR}"
	--outcome "${OUTCOME}"
	--term "${TERM}"
)
printf '  %q ' "${cmd3[@]}"; echo
"${cmd3[@]}"
echo
