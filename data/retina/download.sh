# download raw data
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8295nnn/GSM8295596/suppl/GSM8295596%5FMulti%5FFetal%5F23w1d%5FNR%5Fatac%5Ffragments.tsv.gz # atac fragments
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM8295nnn/GSM8295596/suppl/GSM8295596%5FMulti%5FFetal%5F23w1d%5FNR%5Fatac%5Ffragments.tsv.gz.tbi.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE268nnn/GSE268630/suppl/GSE268630%5FMulti%5FFetal%5F23w1d%5FNR%5Ffiltered%5Ffeature%5Fbc%5Fmatrix.h5 # rna counts

# rename frag file
mv GSM8295596_Multi_Fetal_23w1d_NR_atac_fragments.tsv.gz fragments.tsv.gz
zcat GSM8295596_Multi_Fetal_23w1d_NR_atac_fragments.tsv.gz.tbi.gz > GSM8295596_Multi_Fetal_23w1d_NR_atac_fragments.tsv.gz.tbi
mv GSM8295596_Multi_Fetal_23w1d_NR_atac_fragments.tsv.gz.tbi fragments.tsv.gz.tbi

# extract information from h5 file
Rscript ../../code/extract_matrix.R GSE268630_Multi_Fetal_23w1d_NR_filtered_feature_bc_matrix.h5 ./