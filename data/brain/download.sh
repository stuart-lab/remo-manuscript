
# download raw data
# https://www.10xgenomics.com/datasets/frozen-human-healthy-brain-tissue-3-k-1-standard-2-0-0
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_per_barcode_metrics.csv
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_filtered_feature_bc_matrix.h5
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_atac_fragments.tsv.gz
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/human_brain_3k/human_brain_3k_atac_fragments.tsv.gz.tbi

# rename frag file
mv human_brain_3k_atac_fragments.tsv.gz fragments.tsv.gz
mv human_brain_3k_atac_fragments.tsv.gz.tbi fragments.tsv.gz.tbi

# extract information from h5 file
Rscript ../../code/extract_matrix.R human_brain_3k_filtered_feature_bc_matrix.h5 ./
