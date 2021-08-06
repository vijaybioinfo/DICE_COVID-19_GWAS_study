Inferring potential causal genes for diseases, and a case study on COVID-19.
------------------------------------------------------------------------------

Developers: Job H. Rocha (jrocha@lji.org) and Sourya Bhattacharyya (sourya@lji.org).

Vijayanand and Ay Labs
Division of Vaccine Discovery
La Jolla Institute for Immunology
La Jolla, CA 92037, USA

**************************
This repository shows how to identify potential causal genes relevant
in diseases/traits by taking advantage of 2 statistical methods,
transcriptome-wide association studies (TWAS) and colocalization.

Both TWAS and colocalization use expression quantitative trait loci (eQTL) data of a specific cell type or tissue and
genome wide association studies (GWAS) with respect to the disease/trait of interest
and find the causal genes (colocalization also returns the corresponding variants or SNPs)
jointly associated with both eQTL and GWAS summary statistics.

This repository primarily focusses on analyzing the causal genes of COVID-19.

Expression and eQTL data provided in this repository were taken
from the DICE database of immune cell gene expression, epigenomics
and expression quantitative trait loci (https://dice-database.org).

GWAS summary statistics was downloaded from The COVID-19 Host
Genetics Initiative (https://www.covid19hg.org), release 5 (January 18 2021).
**************************

Description
------------------------
First of all, git clone this repository so you can run
this tutorial in your own system.

You'll find the next 4 main folders:
1.- *example_data* which contains GWAS summary statistics, gene expression and eQTLs data.
2.- *overlap* contains scripts to run the overlap analysis between GWAS and eQTL, and a README file for instructions.
3.- *colocalization* contains data, scripts and a README with detailed instructions for running the colocalization analysis.
4.- *twas* contains scripts and a README file to perform the TWAS analysis.


*Overlap analysis*
This analysis consists in finding overlaping snps that
are both eQTLs and GWAS-associated variants or in LD with
a GWAS-associated variant.

Steps:
1.- Find significant GWAS snps (pval<5e-08).
2.- Get LD snps for all GWAS significant snps.
3.- Overlap GWAS significant + LD snps with eQTL,
taking into account chromosome and position.

To see the detailed information for this analysis, check the
*overlap* folder and underlying README.


*Colocalization analysis*
Colocalization tries to identify causal genetic variants
in a given region that are associated with both input summary statistics, i.e.
both eQTL and GWAS statistics. We have implemented the tool *coloc*
(described in https://github.com/chr1swallace/coloc) which is the most widely used
package for colocalization analysis.

Steps:
1.- Run coloc.abf function using GWAS significant snps and eQTL information.
2.- Summarize and filter results using PP4 > 0.8 and PP4/PP3 ratio > 5.
3.- Keep only colocalized variant that are eQTLs or in LD with an eQTL.

For detailed description, check the *colocalization* folder and underlying README.


*TWAS*
TWAS analysis infers potential causal genes by asking if the effects
in gene expression is relevant in a given disease. TWAS uses
transcriptome prediction models and GWAS information to predict the target genes.

Steps:
1.- Build transcriptome prediction models.
2.- Perform TWAS analysis using prediction models and GWAS summary statistics.
3.- Filter TWAS results by the number of tests performed.

Detailed information for this analysis can be found in
README.md file in the *twas* folder.


Citation
--------------
Please cite the following manuscript if you are using this repository:

Benjamin J. Schmiedel, Vivek Chandra, Job Rocha, Cristian Gonzalez-Colin, Sourya Bhattacharyya, Ariel Madrigal, Christian H. Ottensmeier,  Ferhat Ay, and Pandurangan Vijayanand,
COVID-19 genetic risk variants are associated with expression of multiple genes in diverse immune cell types,
bioRxiv (https://www.biorxiv.org/content/10.1101/2020.12.01.407429v1),
https://doi.org/10.1101/2020.12.01.407429


Contact
--------------
Please email to Job H. Rocha (jrocha@lji.org) and / or Sourya Bhattacharyya (sourya@lji.org) for any questions or suggestions.
