# download raw data
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE232nnn/GSE232222/suppl/GSE232222%5FAJDV174%5Fatac%5Ffragments.tsv.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE232nnn/GSE232222/suppl/GSE232222%5FAJDV174%5Ffiltered%5Ffeature%5Fbc%5Fmatrix.h5

tabix -p bed GSE232222_AJDV174_atac_fragments.tsv.gz

# rename frag file
mv GSE232222_AJDV174_atac_fragments.tsv.gz fragments.tsv.gz
mv GSE232222_AJDV174_atac_fragments.tsv.gz.tbi fragments.tsv.gz.tbi 

# extract information from h5 file
Rscript ../../code/extract_matrix.R GSE232222_AJDV174_filtered_feature_bc_matrix.h5 ./