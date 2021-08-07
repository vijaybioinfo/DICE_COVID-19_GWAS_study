#!/usr/bin/env Rscript

#==========================
# summary of colocalization analysis with respect to a particular GWAS data, and for all the chromosomes
# uses output of Colocalization_Analysis_GWAS_Script.R
#==========================
library(data.table)
library(argparse)

options(scipen = 10)
options(datatable.fread.datatable=FALSE)

# set of chromosomes which would be individually analyzed for testing overlap
CHRLIST_NAMENUM <- c(paste("chr", seq(1,22), sep=""), "chrX", "chrY")


#==========================
args <- args<-ArgumentParser()
args$add_argument("--coloc-dir",default=NULL,type="character",help="Dir where coloc results are saved")
args$add_argument("--ref-eqtls",default=NULL,type="character",help="Reference eQTL file used for colocalization")
args$add_argument("--pp4-threshold",default=0.5,type="double",help="PP H4 threshold")
args<-args$parse_args()

BaseDir <- args$coloc_dir
RefEQTLFile <- args$ref_eqtls
THR_POST_PROB <- args$pp4_threshold


textfile <- paste0(BaseDir, '/Out_Summary_coloc_stat.log')


##===== output summary file to contain the colocalized gene - SNP pairs - first draft
ColocSNPInfoFile <- paste0(BaseDir, '/Out_Summary_Coloc_Gene_SNP_Pairs.bed')


##==== output summary file according to the posterior probability condition
ColocSNPInfoFile_Filt <- paste0(BaseDir, '/Out_Summary_Coloc_Gene_SNP_Pairs_Filtered_PP4.bed')

##==== boolean variable indicating the summary output construction
bool_DF <- FALSE


##====== reference EQTL file
RefEQTLData <- data.table::fread(RefEQTLFile, header=T)
dump_Ref_eQTL_Data <- RefEQTLData[, c("rsnumber", "Ens_gene_id", "chr", "pos", "Distance.to.TSS", "pvalue", "FDR", "beta", "Mean.TPM", "gene_name")]

list_eGenes<-list.files(args$coloc_dir)
list_eGenes<-list_eGenes[grepl("GeneID",list_eGenes)]
list_eGenes<-gsub("GeneID_","",list_eGenes)

dump_Ref_eQTL_Data <- dump_Ref_eQTL_Data[which(dump_Ref_eQTL_Data[, 2] %in% list_eGenes), ]
colnames(dump_Ref_eQTL_Data) <- c('rs_id', 'gene_id', 'chr', 'pos', 'TSSDist', 'pvalue', 'FDR', 'slope', 'meanTPM', 'geneName')

list_eGenes<-unique(dump_Ref_eQTL_Data[,c(2,10,3)])


##====== process each eGene - analyze colocalization results
for (geneIdx in 1:nrow(list_eGenes)) {
	currGeneID <- list_eGenes[geneIdx, 1]
	currGeneName <- list_eGenes[geneIdx, 2]
	CurrGeneDir <- paste0(BaseDir, '/GeneID_', currGeneID)

	##========= TSS of current gene, from the eQTL data - eQTL position + distance from TSS
	idx <- which(dump_Ref_eQTL_Data[,2] == currGeneID)
	TSS_position_currgene <- (dump_Ref_eQTL_Data[idx[1], 4] + dump_Ref_eQTL_Data[idx[1], 5])

	##========= colocalization summary output files
	coloc_post_prob_summary_file <- paste0(CurrGeneDir, '/coloc_abf_summary_DF1.bed')
	coloc_summary_file <- paste0(CurrGeneDir, '/coloc_abf_results_detailed.bed')

	if ((file.exists(coloc_post_prob_summary_file) == TRUE) & (file.exists(coloc_summary_file) == TRUE)) {
			
		##===== posterior probability (of shared variants)
		pp_H0 <- as.numeric(system(paste("awk -F\'[\t]\' \'{if (NR==3) {print $2}}\' ", coloc_post_prob_summary_file), intern = TRUE))
		pp_H1 <- as.numeric(system(paste("awk -F\'[\t]\' \'{if (NR==4) {print $2}}\' ", coloc_post_prob_summary_file), intern = TRUE))
		pp_H2 <- as.numeric(system(paste("awk -F\'[\t]\' \'{if (NR==5) {print $2}}\' ", coloc_post_prob_summary_file), intern = TRUE))
		pp_H3 <- as.numeric(system(paste("awk -F\'[\t]\' \'{if (NR==6) {print $2}}\' ", coloc_post_prob_summary_file), intern = TRUE))
		pp_H4 <- as.numeric(system(paste("awk -F\'[\t]\' \'{if (NR==7) {print $2}}\' ", coloc_post_prob_summary_file), intern = TRUE))

		##====== invalid posterior probability values
		if (is.na(pp_H3) | is.na(pp_H4)) {
			next
		}
                
		##====== gene-SNP pairs from "coloc_summary_file" - pp.H4 > THR_POST_PROB
		##====== they indicate potential causal variant
		coloc_SNP_Set_File <- paste0(CurrGeneDir, '/coloc_SNP_Set.txt')
		system(paste0("awk -F\'[\t]\' \'((NR==1) || (sprintf(\"%0.400f\",$NF)>", THR_POST_PROB, "))\' ", coloc_summary_file, " > ", coloc_SNP_Set_File))
		n <- as.integer(system(paste("cat", coloc_SNP_Set_File, "| wc -l"), intern = TRUE))
		if (n <= 1) {
			next
		}

		##======= for at least one potential causal variant
		##======= list the corresponding variant in the final summary file
		tempDF <- data.table::fread(coloc_SNP_Set_File, header=T)
		SNPID <- as.vector(tempDF[, 2])
		SNPPos <- as.vector(tempDF[, 3])
		SNP_pp_H4 <- as.numeric(tempDF[, ncol(tempDF)])

		# construct the summary data frame containing possible colocalized gene and SNP pairs
		# along with their summary statistics, posterior probabilities
		chr<-list_eGenes[geneIdx, 3]
		currDF <- data.frame(chr=rep(chr,length(SNPID)),geneID=rep(currGeneID, length(SNPID)), geneName=rep(currGeneName, length(SNPID)), TSS=rep(TSS_position_currgene, length(SNPID)), pp_H0_Coloc_Summary=rep(pp_H0, length(SNPID)), pp_H1_Coloc_Summary=rep(pp_H1, length(SNPID)), pp_H2_Coloc_Summary=rep(pp_H2, length(SNPID)), pp_H3_Coloc_Summary=rep(pp_H3, length(SNPID)), pp_H4_Coloc_Summary=rep(pp_H4, length(SNPID)), ColocSNP=SNPID, ColocSNPPos=SNPPos, ColocSNPppH4=SNP_pp_H4)

		if (bool_DF == FALSE) {						
			SNPSummaryDF <- currDF
			bool_DF <- TRUE
		} else {
			SNPSummaryDF <- rbind.data.frame(SNPSummaryDF, currDF)
		}

	}	# end file exist condition

}	# end gene loop



# dump the final summary file
SNPSummaryDF["PP4/PP3"]<-SNPSummaryDF["pp_H4_Coloc_Summary"]/SNPSummaryDF["pp_H3_Coloc_Summary"]
SNPSummaryDF<-SNPSummaryDF[SNPSummaryDF["PP4/PP3"]>5,]
write.table(SNPSummaryDF, ColocSNPInfoFile, row.names=F, col.names=T, sep="\t", quote=F, append=F)

##=== now filter the summary file based on the posterior probability (PP4) threshold (field 9)
##=== this will be the final colocalized set of SNPs
system(paste0("awk -F\'\t\' \'((NR==1) || (sprintf(\"%0.400f\",$9) > ", THR_POST_PROB, "))\' ", ColocSNPInfoFile, " > ", ColocSNPInfoFile_Filt))



