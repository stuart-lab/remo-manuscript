# https://www.10xgenomics.com/datasets/human-jejunum-nuclei-isolated-with-chromium-nuclei-isolation-kit-saltyez-protocol-and-10x-complex-tissue-dp-ct-sorted-and-ct-unsorted-1-standard
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.2/M_Jejunum_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP/M_Jejunum_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_filtered_feature_bc_matrix.h5
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.2/M_Jejunum_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP/M_Jejunum_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz
wget https://cf.10xgenomics.com/samples/cell-arc/2.0.2/M_Jejunum_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP/M_Jejunum_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz.tbi

# rename frag file
mv M_Jejunum_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz fragments.tsv.gz
mv M_Jejunum_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_atac_fragments.tsv.gz.tbi fragments.tsv.gz.tbi

# extract information from h5 file
Rscript ../../code/extract_matrix.R M_Jejunum_Chromium_Nuc_Isolation_vs_SaltyEZ_vs_ComplexTissueDP_filtered_feature_bc_matrix.h5 ./
