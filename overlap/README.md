Overlap of GWAS-associated variants with eQTLs
--------------------------------------------------
The following steps show how to perform a simple
overlap analysis between GWAS and eQTLs
using GWAS-associated snps + snps in LD.


*Input files description*

Plink files/arguments:
1.- "--bfile" prefix for .bed .bim and .fam genotype files.
2.- "--keep"  file containing samples to keep
3.- "--ld-snp-list" a file listing input snps
3.- "--ld-window-r2" r2 threshold
For further information check https://www.cog-genomics.org/plink/2.0/

Overlap.py files/arguments:
1.- "--gwas" gwas summary statistics file (see ../example_data/gwas)
2.- "--ld" plink output file (see output/ld/chr12_out.ld)
3.- "--eqtls" file containin eqtl info (see ../example_data/eqtls/NCM.bed)


*Output Files*

Plink output:
output/ld/chr12/chr12_out.ld contains LD snps for each input snp provided

Overlap.py output:
output/overlap/chr12.csv lists all (Lead + LD) snps that are eqtls along with some additional information


*Commands* (to be executed in bash terminal)

mkdir -p output/ld/input;
mkdir -p output/overlap;
chromosomes=$(awk '($5<5e-08)' ../example_data/gwas/C2_ALL_eur_leave_23andme.tsv | cut -f1 | sort | uniq);

for chrom in $chromosomes; do

   mkdir -p output/ld/$chrom;

   ##keep significant gwas snps for one chromosome
   grep $chrom ../example_data/gwas/C2_ALL_eur_leave_23andme.tsv | awk '($5<5e-08)' | cut -f3 > output/ld/input/$chrom\_input.txt;

   ##Run LD analysis using plink
   plink --bfile {prefix} \
         --keep {your_file} \
         --r2 \
         --ld-snp-list output/ld/input/$chrom\_input.txt \
         --ld-window-kb 250 \
         --ld-window 999999 \
         --ld-window-r2 0.8 \
         --out output/ld/$chrom/$chrom\_out;

   ##Overlap eqtls with gwas + LD snps
   python scripts/Overlap.py --gwas ../example_data/gwas/C2_ALL_eur_leave_23andme.tsv --ld output/ld/$chrom/$chrom\_out.ld --eqtls ../example_data/eqtls/NCM.bed --out output/overlap/$chrom.csv;

done
