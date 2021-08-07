#!/usr/bin/env Rscript

#========================
# this script is for colocalization between eQTL and GWAS summary statistics
#========================

# libraries specifically for colocalization
library(data.table)
library(stringr)
library(coloc)
library(snpStats)
library(ggplot2)
library(ggplotify)
library(argparse)

options(scipen = 10)
options(datatable.fread.datatable=FALSE) 

##======== pp4 threshold for colocalization
##======== not used here - but still worth to note it
THR_POST_PROB <- 0.75

##===== within chr6, excluded MHC region from consideration
##===== as suggested in https://www.nature.com/articles/nature22403
MHC_chr6_LB <- 28477897
MHC_chr6_UB <- 33448354


#=======================
# this function estimates standard error of SNP using beta, MAF values
# references: 1) https://www.biostars.org/p/276869/ (last paragraph), 
# 2) https://www.nature.com/articles/ng.3538, 3) https://www.biostars.org/p/319584/ (last paragraph)
# here b = z / sqrt(2p(1-p)(n+z^2)); where b = beta (regression coefficient), z = z score, p = MAF n = sample size
# so, z can be estimated as  2p(1-p)(b^2)n / sqrt(1 - 2p(1-p)b^2)
# finally standard error is estimated as beta / z
#=======================
estimate_SE <- function(inpdf, beta_column, MAF_column, size_column) {  
        b <- inpdf[, beta_column]
        p <- inpdf[, MAF_column]
        n <- inpdf[, size_column]
        numerator <- (2 * p * (1-p) * (b ^ 2) * n)
        squared_denominator <- (1 - 2 * p * (1-p) * (b ^ 2))

        Estimated_Z <- rep(0, nrow(inpdf))      
        zero_denom_idx <- which(squared_denominator == 0)
        nonzero_denom_idx <- setdiff(seq(1,nrow(inpdf)), zero_denom_idx)
        Estimated_Z[nonzero_denom_idx] <- numerator[nonzero_denom_idx] / sqrt(abs(squared_denominator[nonzero_denom_idx]))
        Estimated_Z[zero_denom_idx] <- numerator[zero_denom_idx]        # avoid zero division 

        # now estimate the standard error (se) which is basically beta / z
        zero_Z_idx <- which(Estimated_Z == 0)
        nonzero_Z_idx <- setdiff(seq(1,length(Estimated_Z)), zero_Z_idx)
        Estimated_SE <- rep(0, length(Estimated_Z))
        Estimated_SE[nonzero_Z_idx] <- b[nonzero_Z_idx] / Estimated_Z[nonzero_Z_idx]
        Estimated_SE[zero_Z_idx] <- b[zero_Z_idx]

        # standard error
        # outdf <- data.frame(SE=Estimated_SE)
        # testing with variance
        # check https://sciencing.com/calculate-variance-standard-error-6372721.html
        outdf <- data.frame(SE=(Estimated_SE * Estimated_SE * n))

        return(outdf)
}

# *************************************************************
# ******** input parameters ***********
# *************************************************************
args <- args<-ArgumentParser()
args$add_argument("--celltype",default=NULL,type="character",help="celltype name")
args$add_argument("--gwas",default=NULL,type="character",help="GWAS summary statistics file")
args$add_argument("--snpinfo",default=NULL,type="character",help="Snp information file for chromosome selected")
args$add_argument("--chr",default=NULL,type="character",help="Chromosome")

args$add_argument("--outdir",default=NULL,type="character",help="Base output directory")
args$add_argument("--ref-eqtls",default=NULL,type="character",help="Reference eqtls file")
args$add_argument("--matrix-eqtl",default=NULL,type="character",help="Matrix eqtl file for celltype and chromosome specified")

args$add_argument("--use-cc",default=1,type="integer",help="Specify 0: to use pvalue and se only, 1: to use case control ratio (cc) only, 2: to use cc, pvalue and se")
args$add_argument("--cc-ratio",default=0.33,type="double",help="Case control ratio")
args$add_argument("--p12",default=1e-05,type="double",help="Prior probability (p12) value")
args<-args$parse_args()



InpCellType<- args$celltype
inp_GWAS_file <- args$gwas
inpchr <- args$chr

BaseOutDir <- args$outdir
RefEQTLFile <- args$ref_eqtls
Complete_MATEQTL_SNP_Gene_Pair_File <- args$matrix_eqtl

Use_case_control_stat <- args$use_cc
case_control_fraction <- args$cc_ratio
p12_val<-args$p12

################################################################
system(paste("mkdir -p", BaseOutDir))

##==== output summary file
textfile <- paste0(BaseOutDir, '/Out_Summary.log')
outtext <- paste0("\n *** Summary statistics for colocalization (GWAS) analysis **** \n cell type : ", InpCellType, "\n input chromosome : ", inpchr, "\n input GWAS summary statistics file : ", inp_GWAS_file, "\n reference EQTL file : ", RefEQTLFile, "\n Complete MATRIXQTL output : ", Complete_MATEQTL_SNP_Gene_Pair_File, "\n Use_case_control_stat : ", Use_case_control_stat, "\n case_control_fraction : ", case_control_fraction, "\n BaseOutDir : ", BaseOutDir)
cat(outtext, file=textfile, append=FALSE, sep="\n")

##==== dump GWAS info for "inpchr"
Total_GWAS_entries <- as.integer(system(paste("cat", inp_GWAS_file, "| wc -l"), intern = TRUE))-1
GWAS_File_CurrChr <- paste0(BaseOutDir, '/dump_GWAS_Statistic_', inpchr, '.bed')
system(paste0("awk -F\'[\t]\' \'{if ((NR==1) || ($1 == \"", inpchr, "\")) {print $0}}\' ", inp_GWAS_file, " > ", GWAS_File_CurrChr))

n <- as.integer(system(paste("cat", GWAS_File_CurrChr, "| wc -l"), intern = TRUE))
if (n <= 1) {
        stop("STOP !!!!!!! No reference GWAS data for the current chromosome !! quit !!! \n", call.=FALSE)
}


GWAS_Data_CurrChr <- data.table::fread(GWAS_File_CurrChr, header=T)
colnames(GWAS_Data_CurrChr) <- c('chr', 'pos', 'rs_id', 'variant_id', 'pval_nominal', 'slope', 'slope_se')
GWAS_Data_CurrChr<-GWAS_Data_CurrChr[,c('chr', 'pos', 'rs_id', 'variant_id', 'pval_nominal', 'slope', 'slope_se')]

outtext <- paste0("\n *** Number of reference GWAS SNPs for the current chromosome: ", nrow(GWAS_Data_CurrChr))
cat(outtext, file=textfile, append=TRUE, sep="\n")


dump_Ref_eQTL_Data <- data.table::fread(RefEQTLFile, header=T)
dump_Ref_eQTL_Data <- dump_Ref_eQTL_Data[, c("rsnumber", "Ens_gene_id", "chr", "pos", "Distance.to.TSS", "pvalue", "FDR", "beta", "Mean.TPM")]
dump_Ref_eQTL_Data <- dump_Ref_eQTL_Data[which(dump_Ref_eQTL_Data[, 3] == inpchr), ]
if (nrow(dump_Ref_eQTL_Data) == 0) {
        stop("STOP !!!!!!! No reference eQTL for the current chromosome !! quit !!! \n", call.=FALSE)
}


colnames(dump_Ref_eQTL_Data) <- c('rs_id', 'gene_id', 'chr', 'pos', 'TSSDist', 'pvalue', 'FDR', 'slope', 'meanTPM')
outtext <- paste0("\n *** Number of reference EQTLs for the current chromosome: ", nrow(dump_Ref_eQTL_Data))
cat(outtext, file=textfile, append=TRUE, sep="\n")


##======== standard deviation of gene expression : possible use in colocalization
stdev_gene_expr <- sd(unique(dump_Ref_eQTL_Data[, ncol(dump_Ref_eQTL_Data)]),na.rm=T)
outtext <- paste0("\n\n ****** Standard deviation of gene expression for the current chromosome: ", stdev_gene_expr)
cat(outtext, file=textfile, append=TRUE, sep="\n")


##======== list of eGenes - for each gene, colocalization analysis will be carried out separately
list_eGenes <- as.vector(unique(dump_Ref_eQTL_Data[, 2]))
outtext <- paste0("\n *** Number of eGenes for the current chromosome: ", length(list_eGenes))
cat(outtext, file=textfile, append=TRUE, sep="\n")


##======== complete SNP information for the current chromosome
SNPInfoData <- data.table::fread(args$snpinfo, header=T, sep=" ")        # employ space separator
colnames(SNPInfoData) <- c('chromosome', 'position', 'variant_id', 'rs_id', 'ref', 'alt', 'AC', 'AF', 'AN')
cat(sprintf("\n ===>> after reading SNPInfoData "))
outtext <- paste0("\n SNP info file (complete SNP set) : ", args$snpinfo, "\n ==>> Number of SNPs in the reference SNP information file for the current chromosome : ", nrow(SNPInfoData), " number of columns for this reference dataset : ", ncol(SNPInfoData))
cat(outtext, file=textfile, append=TRUE, sep="\n")


##======== complete matrixQTL results (significant or not) for every gene-SNP pairs
##======== CIS distance < 1 Mb, default DICE configuration
##======== Note: the file name has "Chr" instead of "chr" - so we used stringr::str_to_title() function
##======== fields: snps    gene    statistic       pvalue  FDR     beta
Complete_MATEQTL_SNP_Gene_Pair_Data <- data.table::fread(Complete_MATEQTL_SNP_Gene_Pair_File, header=T)
colnames(Complete_MATEQTL_SNP_Gene_Pair_Data) <- c('variant_id', 'gene_id', 'statistic', 'pvalue', 'FDR', 'beta')
cat(sprintf("\n ===>> after reading Complete_MATEQTL_SNP_Gene_Pair_Data "))
outtext <- paste0("\n *** Complete MatrixeQTL Association file (matrixQTL output) : ", Complete_MATEQTL_SNP_Gene_Pair_File, "\n Number of MatrixEQTL gene - SNP pairs for the current chromosome : ", nrow(Complete_MATEQTL_SNP_Gene_Pair_Data))
cat(outtext, file=textfile, append=TRUE, sep="\n")


##=========== merged SNP information + matrixeQTL gene-SNP pairs
merge_MATEQTL_SNPInfo_CurrChr_DF <- dplyr::inner_join(Complete_MATEQTL_SNP_Gene_Pair_Data, SNPInfoData, by='variant_id')
cat(sprintf("\n ===>> after creating merge_MATEQTL_SNPInfo_CurrChr_DF "))


##==== extract fields for colocalization
Coloc_CurrChr_DF_SNP <- merge_MATEQTL_SNPInfo_CurrChr_DF[, c("rs_id", "variant_id", "gene_id", "chromosome", "position", "pvalue", "FDR", "beta", "ref", "alt", "AC", "AF", "AN")]
colnames(Coloc_CurrChr_DF_SNP) <- c("rs_id", "variant_id", "gene_id", "chr", "pos", "pvalue", "FDR", "slope", "ref", "alt", "AC", "AF", "AN")   
cat(sprintf("\n ===>> after creating Coloc_CurrChr_DF_SNP"))


##====== estimate standard error using three arguments: beta (column 8), MAF (column 12), size (column 13)
CN <- colnames(Coloc_CurrChr_DF_SNP)
SE_DF <- estimate_SE(Coloc_CurrChr_DF_SNP, 8, 12, 13)
Coloc_CurrChr_DF_SNP <- cbind.data.frame(Coloc_CurrChr_DF_SNP, SE_DF$SE)
colnames(Coloc_CurrChr_DF_SNP) <- c(CN, "slope_se")


outtext <- paste0("\n ==>> Number of SNPs to be used for colocalization for the current chromosome (basically number of rows in Coloc_CurrChr_DF_SNP) : ", nrow(Coloc_CurrChr_DF_SNP))
cat(outtext, file=textfile, append=TRUE, sep="\n")
##====== analyze for each eGene, its associated SNPs and GWAS SNPs for colocalization (i.e. shared variants)
bool_gene_specific_process <- FALSE
for (geneIdx in 1:length(list_eGenes)) {
        currGene <- list_eGenes[geneIdx]
        cat(sprintf("\n processing eGene : %s ", currGene))

        ##==== get TSS from eQTL data = eQTL position + distance from TSS
        idx <- which(dump_Ref_eQTL_Data[,2] == currGene)
        TSS_position_currgene <- (dump_Ref_eQTL_Data[idx[1], 4] + dump_Ref_eQTL_Data[idx[1], 5])
        outtext <- paste0("\n *** Processing gene : ", currGene, "  its TSS : ", TSS_position_currgene)
        cat(outtext, file=textfile, append=TRUE, sep="\n")      

        ##===== SNPs associated with current gene
        ##===== *** IMPORTANT: excluding SNPs belonging to chr6 MHC locus
        Coloc_DF_SNP <- Coloc_CurrChr_DF_SNP[which((Coloc_CurrChr_DF_SNP[,3] == currGene) & (!((Coloc_CurrChr_DF_SNP[,4] == "chr6") & (Coloc_CurrChr_DF_SNP[,5] >= MHC_chr6_LB) & (Coloc_CurrChr_DF_SNP[,5] <= MHC_chr6_UB)))), ]       
        outtext <- paste0("\n *** number of rows in Coloc_DF_SNP structure for the gene : ", currGene, "  is (processed for colocalization) : ", nrow(Coloc_DF_SNP))    
        cat(outtext, file=textfile, append=TRUE, sep="\n")
        if (nrow(Coloc_DF_SNP) == 0) {
                next
        }

        ##===== colocalization between the GWAS SNPs and the dumped SNPs
        merge_SNP_GWAS_Data <- merge(Coloc_DF_SNP, GWAS_Data_CurrChr, by=c("chr", "pos"), all=FALSE, suffixes=c("_snp","_gwas"))                
                
        if (nrow(merge_SNP_GWAS_Data) == 0) {
                next
        }

        ##===== check for non-NA fields
        NA_idx <- which((is.na(merge_SNP_GWAS_Data$pval_nominal)) | (is.na(merge_SNP_GWAS_Data$pvalue)) | (is.na(merge_SNP_GWAS_Data$slope_snp)) | (is.na(merge_SNP_GWAS_Data$slope_se_snp)) | (is.na(merge_SNP_GWAS_Data$AF)))
        if (length(NA_idx) > 0) {
                merge_SNP_GWAS_Data <- merge_SNP_GWAS_Data[-c(NA_idx), ]
        }
        if (nrow(merge_SNP_GWAS_Data) == 0) {
                next
        }

        ##=== dump the input SNPs for colocalization analysis in the specified output directory 
        ##=== for processing the current gene
        CurrGene_OutDir <- paste0(BaseOutDir, '/GeneID_', currGene)
        system(paste("mkdir -p", CurrGene_OutDir))
        write.table(Coloc_DF_SNP, paste0(CurrGene_OutDir, '/dump_SNPs_for_colocalization.txt'), row.names=F, col.names=T, sep="\t", quote=F, append=F)
        write.table(merge_SNP_GWAS_Data, paste0(CurrGene_OutDir, '/merged_SNPs_GWAS_input_colocalization.txt'), row.names=F, col.names=T, sep="\t", quote=F, append=F)
        

        ##=================
        ## main colocalization routine
        ##=================     
        if (Use_case_control_stat == 1) {
                ##====== cc for dataset1, 
                ##====== N: number of SNPs for this gene
                dataset1=list(pvalues=merge_SNP_GWAS_Data$pval_nominal, type="cc", s=case_control_fraction, N=nrow(GWAS_Data_CurrChr))
        } else if (Use_case_control_stat == 2) {
                ##====== cc for dataset1, 
                ##====== N: number of SNPs for this gene
                ##====== but here we also specify the GWAS data beta and varbeta stats
                dataset1=list(pvalues=merge_SNP_GWAS_Data$pval_nominal, type="cc", s=case_control_fraction, N=nrow(GWAS_Data_CurrChr), beta=merge_SNP_GWAS_Data$slope_gwas, varbeta=(merge_SNP_GWAS_Data$slope_se_gwas * merge_SNP_GWAS_Data$slope_se_gwas))
        } else {
                ##==== quant for dataset1
                ##====== N: number of SNPs for this gene
                ## from the standard error, compute the variance, 
                ## by Multiplying the square of the standard error by the sample size
                dataset1=list(pvalues=merge_SNP_GWAS_Data$pval_nominal, type="quant", N=nrow(GWAS_Data_CurrChr), beta=merge_SNP_GWAS_Data$slope_gwas, varbeta=(merge_SNP_GWAS_Data$slope_se_gwas * merge_SNP_GWAS_Data$slope_se_gwas))
        }

        ## quant for dataset2, but using stdev of gene expression
        ## N: number of SNPs for this chromosome (used for colocalization)
        dataset2 <- list(pvalues=merge_SNP_GWAS_Data$pvalue, type="quant", N=nrow(Coloc_DF_SNP), beta=merge_SNP_GWAS_Data$slope_snp, varbeta=(merge_SNP_GWAS_Data$slope_se_snp * merge_SNP_GWAS_Data$slope_se_snp), sdY=stdev_gene_expr)

        ## colocalization function according to prior probability settings
        outtext <- paste0("\n *** Using prior probability p12_val : ", p12_val)  
        cat(outtext, file=textfile, append=TRUE, sep="\n")              
        result <- coloc.abf(dataset1=dataset1, dataset2=dataset2, MAF=merge_SNP_GWAS_Data$AF, p12=p12_val)
 

 	##======== colocalization results - plot sensitivity of the posterior probability computation
 	if (!is.null(result)){
        	 pdf(paste0(CurrGene_OutDir, "/coloc_abf_res_sensitivity_plot.pdf"))
         	try(sensitivity(result,rule="H4 > 0.5"))
         	dev.off()
 	}

	##======== colocalization results - two data frames
 	##======== let N = nrow(merge_eQTL_GWAS_Data)
	##======== result$summary: summary statistic of posterior probability (considering all SNPs)
 	##======== result$results: detailed statistic of posterior probability per SNP
 	##======== Column 1: SNP.NUM where NUM = index number of SNP - varies from 1 to N
 	write.table(result$summary, paste0(CurrGene_OutDir, '/coloc_abf_summary_DF1.bed'), row.names=T, col.names=T, sep="\t", quote=F, append=F)
 	write.table(result$results, paste0(CurrGene_OutDir, '/coloc_abf_results_DF2.bed'), row.names=F, col.names=T, sep="\t", quote=F, append=F)

 	##======== modify the results (DF2) file to add SNP IDs in addition to the SNP.NUM format
 	resDF <- result$results         
 	CN <- colnames(resDF)
 	temp_SNP_vec <- as.vector(paste0("SNP.", seq(1,nrow(merge_SNP_GWAS_Data))))
 	m <- match(resDF[,1], temp_SNP_vec)
 	idx_r <- which(!is.na(m))
 	idx_t <- m[!is.na(m)]
 
 	# when data frame merging was done by chr and pos fields
 	resDF <- cbind.data.frame(resDF[idx_r, 1], merge_SNP_GWAS_Data[idx_t, c("rs_id_snp", "pos")], resDF[idx_r, 2:ncol(resDF)])
 
 	colnames(resDF) <- c(CN[1], 'rs_id', 'pos', CN[2:length(CN)])
 	write.table(resDF, paste0(CurrGene_OutDir, '/coloc_abf_results_detailed.bed'), row.names=F, col.names=T, sep="\t", quote=F, append=F)
 
 	outtext <- paste0("\n *** creating detailed results -- entries in resDF : ", nrow(resDF), "  matching SNP entries : ", length(idx_r))
 	cat(outtext, file=textfile, append=TRUE, sep="\n")      

}

##===== delete temporary files
if (file.exists(GWAS_File_CurrChr)) {
        system(paste("rm", GWAS_File_CurrChr))
}


