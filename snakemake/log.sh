CODEDIR=$1
OUTDIR=$2

# create conda environment from env yml file 
conda env create -f abc-max.yml

# replace directories with own directories
# replace projectDir 
sed -i s+/oak/stanford/groups/akundaje/kmualim/test/ABC-Max-pipeline/+$CODEDIR+ ABC-Max.example.json
# replace data directory 
sed -i s+/oak/stanford/groups/akundaje/projects/ABC_links/GWAS_test/Test_data/+$OUTDIR+ ABC-Max.example.json 
sed -i s+/oak/stanford/groups/akundaje/projects/ABC_links/GWAS_test/Test_data/+$OUTDIR+ ABC-Max.config-preds.tsv
sed -i s+/oak/stanford/groups/akundaje/projects/ABC_links/GWAS_test/Test_data/+$OUTDIR+ ABC-Max.config-traits.tsv
