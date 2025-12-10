# download pbmc

wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz
mv pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz fragments.tsv.gz

wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz.tbi
mv pbmc_granulocyte_sorted_3k_atac_fragments.tsv.gz.tbi fragments.tsv.gz.tbi

wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_atac_peaks.bed
mv pbmc_granulocyte_sorted_3k_atac_peaks.bed peaks.bed

wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/pbmc_granulocyte_sorted_3k/pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5
mv pbmc_granulocyte_sorted_3k_filtered_feature_bc_matrix.h5 counts.h5

Rscript ../../../code/extract_matrix.R counts.h5 ./