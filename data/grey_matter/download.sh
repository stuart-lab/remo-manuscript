### download raw data
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE193nnn/GSE193240/suppl/GSE193240%5FPD003%5Fatac%5Ffragments.tsv.gz 
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE193nnn/GSE193240/suppl/GSE193240%5FPD003%5Fatac%5Ffragments.tsv.gz.tbi.gz 
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE193nnn/GSE193240/suppl/GSE193240%5FPD003%5Ffiltered%5Ffeature%5Fbc%5Fmatrix.h5 

# rename frag file
mv GSE193240_PD003_atac_fragments.tsv.gz fragments.tsv.gz
zcat GSE193240_PD003_atac_fragments.tsv.gz.tbi.gz > GSE193240_PD003_atac_fragments.tsv.gz.tbi
mv GSE193240_PD003_atac_fragments.tsv.gz.tbi fragments.tsv.gz.tbi

# extract information from h5 file
Rscript ../../code/extract_matrix.R GSE193240_PD003_filtered_feature_bc_matrix.h5 ./