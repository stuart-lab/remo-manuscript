# Output Files
wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_unsorted_10k/pbmc_unsorted_10k_per_barcode_metrics.csv
wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_unsorted_10k/pbmc_unsorted_10k_filtered_feature_bc_matrix.h5
wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_unsorted_10k/pbmc_unsorted_10k_atac_fragments.tsv.gz
wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_unsorted_10k/pbmc_unsorted_10k_atac_fragments.tsv.gz.tbi

# rename frag files
mv pbmc_unsorted_10k_atac_fragments.tsv.gz fragments.tsv.gz
mv pbmc_unsorted_10k_atac_fragments.tsv.gz.tbi fragments.tsv.gz.tbi

# extract information from h5 file
Rscript ../../code/extract_matrix.R pbmc_unsorted_10k_filtered_feature_bc_matrix.h5 ./