# brain cancer 
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206579

# seurat object download - https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE206579
wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/GSE206579/suppl/GSE206579%5FsnATAC%5FPF%5FEPN%5FSeurat.rds.gz
gzip -d GSE206579_snATAC_PF_EPN_Seurat.rds.gz

# fragment files, processed using cellranger-atac
mkdir -p I2
cd I2
aws s3 cp s3://stuartlab/REMO/disease/brain_cancer/I2_fragments.tsv.gz ./
aws s3 cp s3://stuartlab/REMO/disease/brain_cancer/I2_fragments.tsv.gz.tbi ./
mv I2_fragments.tsv.gz fragments.tsv.gz
mv I2_fragments.tsv.gz.tbi fragments.tsv.gz.tbi
cd ..

mkdir -p I3
cd I3
aws s3 cp s3://stuartlab/REMO/disease/brain_cancer/I3_fragments.tsv.gz ./
aws s3 cp s3://stuartlab/REMO/disease/brain_cancer/I3_fragments.tsv.gz.tbi ./
mv I3_fragments.tsv.gz fragments.tsv.gz
mv I3_fragments.tsv.gz.tbi fragments.tsv.gz.tbi
cd ..

mkdir -p I4
cd I4
aws s3 cp s3://stuartlab/REMO/disease/brain_cancer/I4_fragments.tsv.gz ./
aws s3 cp s3://stuartlab/REMO/disease/brain_cancer/I4_fragments.tsv.gz.tbi ./
mv I4_fragments.tsv.gz fragments.tsv.gz
mv I4_fragments.tsv.gz.tbi fragments.tsv.gz.tbi
cd ..

mkdir -p I5
cd I5
aws s3 cp s3://stuartlab/REMO/disease/brain_cancer/I5_fragments.tsv.gz ./
aws s3 cp s3://stuartlab/REMO/disease/brain_cancer/I5_fragments.tsv.gz.tbi ./
mv I5_fragments.tsv.gz fragments.tsv.gz
mv I5_fragments.tsv.gz.tbi fragments.tsv.gz.tbi
cd ..

mkdir -p M7
cd M7
aws s3 cp s3://stuartlab/REMO/disease/brain_cancer/M7_fragments.tsv.gz ./
aws s3 cp s3://stuartlab/REMO/disease/brain_cancer/M7_fragments.tsv.gz.tbi ./
mv M7_fragments.tsv.gz fragments.tsv.gz
mv M7_fragments.tsv.gz.tbi fragments.tsv.gz.tbi
cd ..

mkdir -p M8
cd M8
aws s3 cp s3://stuartlab/REMO/disease/brain_cancer/M8_fragments.tsv.gz ./
aws s3 cp s3://stuartlab/REMO/disease/brain_cancer/M8_fragments.tsv.gz.tbi ./
mv M8_fragments.tsv.gz fragments.tsv.gz
mv M8_fragments.tsv.gz.tbi fragments.tsv.gz.tbi
cd ..