# Prerequisites
### Clone projects and restore environment:
```r
git clone https://github.com/stuart-lab/metapeak-analysis.git
cd metapeak-analysis/linkage_analysis
renv::restore()
```
For restoring renv refer to https://rstudio.github.io/renv/articles/renv.html#collaboration

### Prepare resources folder
See below for expected structure:
```bash
./linkage_analysis/resources
├── example
│   └── pgls_metadata_all_chr.tsv (from supplementary data)
├── gencode.v32.basic.annotation.gtf.gz (from https://www.gencodegenes.org/human/release_32.html)
└── multiome (your Seurat multiome datasets[.rds])
```

# Script setup
##### Modify slurm parameters as needed, e.g.(`cpus-per-task`, `mem-per-cpu`)
- `/scripts/01_link_peaks.sh`
- `/scripts/02_run_aggregate_links.sh`
- `/scripts/03_run_probability_analysis.sh`

##### copy pwd to `PROJECT_ROOT` of .sh files
```bash
pwd
$ /metapeak-analysis/linkage_analysis/
```

# Reproduce results
In `/scripts/03_run_probability_analysis.sh`, change `SIG_PGL_FPATH` variable to `"${PROJECT_ROOT}/resources/example/pgl_metadata_all_chr.tsv"`

Run the following script
```bash
./scripts/03_run_probability_analysis.sh
```

# Run full analysis workflow on user-supplied data
### 1. Replace examples files in `/resources/multiome/` with your intended Seurat multiome objects (`.rds`)

### 2. Modify the following variables in `01_link_peaks.sh` according to each seurat multiome object
- `OUTPUT_DIR` (links for each chromosome are stored here, in `/data/multiome_links/dataset_name`)
- `SEURAT_RDS` (filepath for Seurat multiome object, in `/resources/multiome`)

E.g. 
```bash
OUTPUT_DIR="${PROJECT_ROOT}/data/multiome_links/jejunum"
SEURAT_RDS="${PROJECT_ROOT}/resources/multiome/jejunum_multiome_linkpeaks.rds"
```

### 3. Run the follow script:
```bash
./01_link_peaks.sh
```
Repeat for each Seurat multiome object.

### 4. Create the following config files in `/scripts/configs`
##### `objects_and_links.csv` (modify example table of seurat obj filepaths and their mutiome_link dirs)
```csv
seurat_obj_fpath,links_dir
/resources/multiome/jejunum_multiome_linkpeaks.rds,/data/multiome_links/jejunum/
/resources/multiome/bile_duct_multiome_linkpeaks.rds,/data/multiome_links/bile_duct/
```

##### `link_filtering_conds.csv` (specify p-value and score thresholds to filter peak-gene links for downstream analysis)
```csv
condition,pvalue,score
sig_pvalue,0.05,0
hsig_pvalue_high_coeff,0.01,0.1
```

### 5. Run the following script:
```bash
./scripts/02_run_aggregate_links.sh
```

### 6. Modify `/scripts/03_run_probability_analysis.sh`
Specify the `SIG_PGL_FPATH` variable as the path to the combined pgl_metadata for the condition of choice (e.g. `sig_pvalue`), example below:
```bash
${PROJECT_ROOT}/data/filtered_links/sig_pvalue/pgl_metadata/pgls_metadata_all_chr.tsv
```

### 7. Run the following script:
```bash
./scripts/03_run_probability_analysis.sh
```

# Results

Output structure:
```bash
./linkage_analysis/data
├── combined_metadata
├── filtered_links
├── multiome_links
└── probability_analysis
```

A table of conditional probabilities can be found at:
```
data/probability_analysis/conditional_prob_unlinked_singleCRE.tsv
```

The per-chromosome odds ratio for the likelihood that CRE pairs linked to the same gene are within the same REMO module can be found at:
```
data/probability_analysis/logistic_regression_results/glm_results_all.tsv
```