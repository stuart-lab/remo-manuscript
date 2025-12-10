# download islet from S3

aws s3 cp s3://stuartlab/REMO/annotate/islet/ ./ --recursive

# extract information from h5 file
Rscript ../../../code/extract_matrix.R filtered_feature_bc_matrix.h5 ./