# download raw data
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5391nnn/GSM5391776/suppl/GSM5391776%5FHF1M%2DCryo%5Fatac%5Ffragments.tsv.gz
wget https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM5391nnn/GSM5391776/suppl/GSM5391776%5FHF1M%2DCryo%5Ffiltered%5Ffeature%5Fbc%5Fmatrix.h5

tabix -p bed GSM5391776_HF1M-Cryo_atac_fragments.tsv.gz

# rename frag file
mv GSM5391776_HF1M-Cryo_atac_fragments.tsv.gz fragments.tsv.gz
mv GSM5391776_HF1M-Cryo_atac_fragments.tsv.gz.tbi fragments.tsv.gz.tbi

# extract information from h5 file
Rscript ../../code/extract_matrix.R GSM5391776_HF1M-Cryo_filtered_feature_bc_matrix.h5 ./