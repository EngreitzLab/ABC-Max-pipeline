OUTDIR=$1

# create conda environment from env yml file 
conda env create -f abc-max.yml

# replace directories with own directories 
sed -i s+/oak/stanford/groups/akundaje/projects/ABC_links/GWAS_test/Test_data/+$OUTDIR+ ABC-Max.example.json 

