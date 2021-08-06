TWAS Analysis
---------------------
GWAS studies have successfully identified
genetic loci for several trait, however, the target
molecules on which they exert their effects are still
unknown. Transcriptome-wide association studies (TWAS)
have been developed to address this problem by
identifying putative causal genes. In brief, TWAS tests
the effects of gene expression on phenotypes by using
transcriptome prediction models and GWAS summary statistics.

For further information, read: https://www.nature.com/articles/s41467-018-03621-1


Prediction Models
--------------------
The following steps build prediction models
based on expression and genotype data using
an Elastic Net algorithm with a 10-fold
cross-validation. This is performed in
2 steps, model training and make db file step.

To build prediction models you will need:
1.- Model Training script (see https://github.com/gamazonlab/MR-JTI/blob/master/model_training/predixcan/src/predixcan_r.r)
2.- "--plink_file_name" prefix for 3 genotype files (.bed, .bim and .fam)
3.- "--expression_file_name" An expression matrix with donor id as columns and genes as rows (see twas_example_data/prediction_models/expression/NCM.txt)
4.- "--annotation_file_name" A gene annotation file (see twas_example_data/prediction_models/annotation/gene_annotation.txt)

*Model Training*

mkdir -p output/prediction_models/NCM;

for ID in $(seq 1 $(echo $(($(wc -l twas_example_data/prediction_models/expression/NCM.txt | cut -f1 -d' ')/100))));
        do
        Rscript predixcan_r.r \
        --model_training \
                --main_dir output/prediction_models/NCM \
                --plink_file_name {prefix} \
                --expression_file_name twas_example_data/prediction_models/expression/NCM.txt \
                --subjob_id $ID \
                --n_genes_for_each_subjob 100 \
                --annotation_file_name twas_example_data/prediction_models/annotation/gene_annotation.txt \
                --n_folds 10 \
                --gene_type protein_coding,lincRNA,miRNA;
        done

*Make db file*

Rscript predixcan_r.r \
         --generate_db_and_cov \
         --main_dir  output/prediction_models/NCM \
         --plink_file_name {prefix} \
         --expression_file_name twas_example_data/prediction_models/expression/NCM.txt \
         --annotation_file_name twas_example_data/prediction_models/annotation/gene_annotation.txt \
         --output_file_name NCM;




TWAS
-------------------
To perform TWAS analysis the following files
and software are required:
1.- MetaXcan software (see https://github.com/hakyimlab/MetaXcan)
2.- "--model_db_path" A db file of the prediction model (see output/prediction_models/NCM/output/NCM.db, also see steps above)
3.- "--covariance" A covariance file of the prediction model (see see output/prediction_models/NCM/output/NCM.cov, also see steps above)
4.- "--gwas_file" A gwas summary statistics file (see ../example_data/gwas/C2_ALL_eur_leave_23andme.tsv)

*Perform TWAS*

python SPrediXcan.py \
        --model_db_path output/prediction_models/NCM/output/NCM.db \
        --covariance output/prediction_models/NCM/output/NCM.cov \
        --gwas_file ../example_data/gwas/C2_ALL_eur_leave_23andme.tsv \
        --snp_column rsid \
        --effect_allele_column alt \
        --non_effect_allele_column ref \
        --beta_column beta \
        --pvalue_column pval \
        --output_file output/TWAS/C2_ALL_eur_leave_23andme/NCM.asso.csv;


*Filter results*

To filter TWAS results, by the effective
number of gene-model pairs tested, run the
command below. 'twas_folder' is the path
to all TWAS results, 'models_folder' is
the path where the prediction models are
saved and 'models' a comma separated list
of models used for the TWAS analyses.

python scripts/filter_twas.py --twas_folder output/TWAS/C2_ALL_eur_leave_23andme --models_folder output/prediction_models --models NCM --out output/filtered_TWAS/C2_ALL_eur_leave_23andme.csv
