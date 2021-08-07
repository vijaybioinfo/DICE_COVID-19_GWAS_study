Colocalization Analysis
-------------------------------------
Colocalization analysis infers whether
two phenotypes (having summary statistics) are likely to be under
the influence of the same causal genetic
variant in a given region. *coloc* is the most cited tool to perform
colocalization analysis. This repository contains a working implementation of *coloc*

For further reading see:
https://pubmed.ncbi.nlm.nih.gov/19039033/
https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004383
https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1008720

Colocalization
----------------------------
To perform colocalization analysis
using eQTLs and GWAS summary statistics,

*first install the following R libraries* (type the following command in R prompt):

install.packages(c('data.table', 'stringr', 'snpStats', 'ggplot2', 'ggplotify', 'argparse', 'coloc'))

*Then run the commands below*

Parameters description:

1.- "--celltype": celltype name

2.- "--gwas": gwas summary statistics (see output/gwas_input/C2_ALL_eur_leave_23andme_round5.tsv)
3.- "--snpinfo": file containing chromosome. position. reference and alternative allele for each snp (see coloc_example_data/snpinfo/snpinfo_chr12.txt)
4.- "--chr": chromosome
5.- "--outdir": output directory
6.- "--ref-eqtls": file containg eqtls for the celltype specified (see ../example_data/refeqtls/NCM.bed)
7.- "--matrix-eqtl": Matrix eqtl output file for celltype and chromosome specified (see coloc_example_data/matrix_eqtl/Output_Chr12_all_cis.txt)
8.- "--use-cc": Use 0 to run the analysis using pvalue and standard error (se), 1 to use case control ratio, or 2 to use pval, se and case control ratio. We recommend using 0 for this parameter.
9.- "--cc-ratio": Case control ratio (only used if --use-cc is 1 or 2), obtained from the reference GWAS studies.
10.- "--p12": prior probability of a variant being an eQTL and a GWAS-associated variant. This is an optional parameter. By default p12 = 1e-5 (recommended). Users are advised to not alter this parameter unless they are absolutely sure.


*get GWAS significant snps*
mkdir -p output/gwas_input/;
mkdir -p output/chr12/;

awk '($5<5e-08)' ../example_data/gwas/C2_ALL_eur_leave_23andme.tsv  > output/gwas_input/C2_ALL_eur_leave_23andme.tsv

*Run Coloc*
Rscript scripts/Colocalization_Analysis_GWAS_Script.R \
        --celltype NCM \
        --gwas output/gwas_input/C2_ALL_eur_leave_23andme.tsv \
        --chr chr12 \
        --outdir output/chr12 \
        --ref-eqtls ../example_data/eqtls/NCM.bed \
        --matrix-eqtl coloc_example_data/matrix_eqtl/Output_Chr12_all_cis.txt \
        --use-cc 0 \
        --snpinfo coloc_example_data/snpinfo/snpinfo_chr12.txt;


*Summarize coloc results*
Parameters description:
1.- "--coloc-dir" output folder specified in colocalization analysis (see command above)
2.- "--ref-eqtls" file containg eqtls for the celltype specified (see example_data/refeqtls/NCM.bed)


Rscript scripts/Colocalization_Analysis_GWAS_Summary_Script.R \
        --coloc-dir output/chr12 \
        --ref-eqtls ../example_data/eqtls/NCM.bed;


*Filter results*
To keep only colocalized snps that
are eQTLs or in LD with an eQTL run
the following commands:

Parameters description:
1.- "--bfile" prefix for .bed .bim and .fam genotype files.
2.- "--keep"  file containing samples to keep
3.- ld_snp_list: a file with the list of colocalized snps for one gene

*get colocalized snps for gene OAS1*
mkdir -p output/plink/;

grep OAS1 output/chr12/Out_Summary_Coloc_Gene_SNP_Pairs_Filtered_PP4.bed | cut -f10 > output/plink/input.txt;

*get LD snps for colocalized snps (for further information check https://www.cog-genomics.org/plink/)*
plink --bfile {prefix} --keep {your_file} --r2 --ld-snp-list output/plink/input.txt --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0.8 --out output/plink/LD_out

*get colocalized snps that are eQTLs or in LD with an eQTL*
eqtls=$(grep -wFf <(cut -f3 ../example_data/eqtls/NCM.bed | sort | uniq ) output/plink/LD_out.ld | sed -r "s/\s+/\t/g" | cut -f4 | sort | uniq)

(head -n1 output/chr12/Out_Summary_Coloc_Gene_SNP_Pairs_Filtered_PP4.bed && grep -f <(echo $eqtls | sed -r 's/\s+/\n/g') output/chr12/Out_Summary_Coloc_Gene_SNP_Pairs_Filtered_PP4.bed) > output/chr12/Coloc_results.bed
