# fallopian tube
# get files from S3
aws s3 cp s3://stuartlab/ENCODE/fallopian_tube/atac_fragments.tsv.gz ./
aws s3 cp s3://stuartlab/ENCODE/fallopian_tube/atac_fragments.tsv.gz.tbi ./
aws s3 cp s3://stuartlab/ENCODE/fallopian_tube/filtered_feature_bc_matrix.h5 ./

# rename frag file
mv atac_fragments.tsv.gz fragments.tsv.gz
mv atac_fragments.tsv.gz.tbi fragments.tsv.gz.tbi

# extract information from h5 file
Rscript ../../code/extract_matrix.R filtered_feature_bc_matrix.h5 ./
