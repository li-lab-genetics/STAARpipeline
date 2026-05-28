#' Gene-centric analysis of long noncoding RNA (ncRNA) category using STAAR procedure
#'
#' The \code{ncRNA} function takes in chromosome, gene name,
#' the object of opened annotated GDS file, and the object from fitting the null model to analyze the association between a
#' quantitative/dichotomous phenotype (including imbalanced case-control design) and the exonic and splicing category of an ncRNA gene by using STAAR procedure.
#' For each ncRNA category, the STAAR-O p-value is a p-value from an omnibus test
#' that aggregated SKAT(1,25), SKAT(1,1), Burden(1,25), Burden(1,1), ACAT-V(1,25),
#' and ACAT-V(1,1) together with p-values of each test weighted by each annotation
#' using Cauchy method. For imbalance case-control setting, the results correspond to the STAAR-B p-value, which is a p-value from
#' an omnibus test that aggregated Burden(1,25) and Burden(1,1) together with p-values of each test weighted by each annotation using Cauchy method.
#' For multiple phenotype analysis (\code{obj_nullmodel$n.pheno > 1}),
#' the results correspond to multi-trait association p-values (e.g. MultiSTAAR-O) by leveraging
#' the correlation structure between multiple phenotypes.
#' For ancestry-informed analysis, the results correspond to ensemble p-values across base tests, 
#' with the option to return a list of base weights and p-values for each base test.
#' @param chr chromosome.
#' @param gene_name name of the ncRNA gene to be analyzed using STAAR procedure.
#' @param genofile an object of opened annotated GDS (aGDS) file.
#' @param obj_nullmodel an object from fitting the null model, which is either the output from \code{\link{fit_nullmodel}} function,
#' or the output from \code{fitNullModel} function in the \code{GENESIS} package and transformed using the \code{\link{genesis2staar_nullmodel}} function.
#' @param rare_maf_cutoff the cutoff of maximum minor allele frequency in
#' defining rare variants (default = 0.01).
#' @param rv_num_cutoff the cutoff of minimum number of variants of analyzing
#' a given variant-set (default = 2).
#' @param rv_num_cutoff_max the cutoff of maximum number of variants of analyzing
#' a given variant-set (default = 1e+09).
#' @param rv_num_cutoff_max_prefilter the cutoff of maximum number of variants
#' before extracting the genotype matrix (default = 1e+09).
#' @param QC_label channel name of the QC label in the GDS/aGDS file (default = "annotation/filter").
#' @param variant_type type of variant included in the analysis. Choices include "SNV", "Indel", or "variant" (default = "SNV").
#' @param geno_missing_imputation method of handling missing genotypes. Either "mean" or "minor" (default = "mean").
#' @param Annotation_dir channel name of the annotations in the aGDS file \cr (default = "annotation/info/FunctionalAnnotation").
#' @param Annotation_name_catalog a data frame containing the name and the corresponding channel name in the aGDS file.
#' @param Use_annotation_weights use annotations as weights or not (default = TRUE).
#' @param Annotation_name a vector of annotation names used in STAAR (default = NULL).
#' @param SPA_p_filter logical: are only the variants with a normal approximation based p-value smaller than a pre-specified threshold use the SPA method to recalculate the p-value, only used for imbalanced case-control setting (default = TRUE).
#' @param p_filter_cutoff threshold for the p-value recalculation using the SPA method, only used for imbalanced case-control setting (default = 0.05).
#' @param use_ancestry_informed logical: is ancestry-informed association analysis used to estimate p-values (default = FALSE).
#' @param find_weight logical: should the ancestry group-specific weights and weighting scenario-specific p-values for each base test be saved as output (default = FALSE).
#' @param silent logical: should the report of error messages be suppressed (default = FALSE).
#' @return A data frame containing the STAAR p-values (including STAAR-O), or AI-STAAR p-values under ancestry-informed analysis, corresponding to the exonic and splicing category of the given ncRNA gene.
#' If \code{find_weight} is TRUE, returns a list containing the AI-STAAR p-values corresponding to the exonic and splicing category of the given ncRNA gene, as well as the ensemble weights under two sampling scenarios 
#' and p-values under scenarios 1, 2, and combined for each base test.
#' @references Li, Z., Li, X., et al. (2022). A framework for detecting
#' noncoding rare-variant associations of large-scale whole-genome sequencing
#' studies. \emph{Nature Methods}, \emph{19}(12), 1599-1611.
#' (\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})
#' @references Li, X., Li, Z., et al. (2020). Dynamic incorporation of multiple
#' in silico functional annotations empowers rare variant association analysis of
#' large whole-genome sequencing studies at scale. \emph{Nature Genetics}, \emph{52}(9), 969-983.
#' (\href{https://doi.org/10.1038/s41588-020-0676-4}{pub})
#' @export

ncRNA <- function(chr,gene_name,genofile,obj_nullmodel,
                  rare_maf_cutoff=0.01,rv_num_cutoff=2,rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                  QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                  Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                  Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                  SPA_p_filter=TRUE,p_filter_cutoff=0.05,use_ancestry_informed=FALSE,find_weight=FALSE,silent=FALSE){

	## evaluate choices
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	phenotype.id <- as.character(obj_nullmodel$id_include)
	n_pheno <- obj_nullmodel$n.pheno

	## SPA status
	if(!is.null(obj_nullmodel$use_SPA))
	{
		use_SPA <- obj_nullmodel$use_SPA
	}else
	{
		use_SPA <- FALSE
	}

	## get SNV id
	filter <- seqGetData(genofile, QC_label)
	if(variant_type=="variant")
	{
		SNVlist <- filter == "PASS"
	}

	if(variant_type=="SNV")
	{
		SNVlist <- (filter == "PASS") & isSNV(genofile)
	}

	if(variant_type=="Indel")
	{
		SNVlist <- (filter == "PASS") & (!isSNV(genofile))
	}

	variant.id <- seqGetData(genofile, "variant.id")

	rm(filter)
	gc()

	## ncRNA SNVs
	GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
	is.in <- ((GENCODE.Category=="ncRNA_exonic")|(GENCODE.Category=="ncRNA_exonic;splicing")|(GENCODE.Category=="ncRNA_splicing"))&(SNVlist)

	variant.id.ncRNA <- variant.id[is.in]

	rm(GENCODE.Category)
	gc()

	seqSetFilter(genofile,variant.id=variant.id.ncRNA,sample.id=phenotype.id)

	rm(variant.id.ncRNA)
	gc()

	GENCODE.Info <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Info")]))
	GENCODE.Info.split <- strsplit(GENCODE.Info, split = "[;]")
	Gene <- as.character(sapply(GENCODE.Info.split,function(z) gsub("\\(.*\\)","",z[1])))

	Gene_list_1 <- as.character(sapply(strsplit(Gene,','),'[',1))
	Gene_list_2 <- as.character(sapply(strsplit(Gene,','),'[',2))
	Gene_list_3 <- as.character(sapply(strsplit(Gene,','),'[',3))

	rm(GENCODE.Info)
	gc()

	rm(GENCODE.Info.split)
	gc()

	variant.id.ncRNA <- seqGetData(genofile, "variant.id")

	seqResetFilter(genofile)

	### Gene
	is.in <- union(which(Gene_list_1==gene_name),which(Gene_list_2==gene_name))
	is.in <- union(is.in,which(Gene_list_3==gene_name))

	variant.is.in <- variant.id.ncRNA[is.in]

	seqSetFilter(genofile,variant.id=variant.is.in,sample.id=phenotype.id)

	## genotype id
	id.genotype <- seqGetData(genofile,"sample.id")
	# id.genotype.match <- rep(0,length(id.genotype))

	id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
	phenotype.id.merge <- data.frame(phenotype.id)
	phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
	id.genotype.match <- phenotype.id.merge$index

	## Genotype
	Geno <- NULL
	if(length(seqGetData(genofile, "variant.id"))<rv_num_cutoff_max_prefilter)
	{
		Geno <- seqGetData(genofile, "$dosage")
		Geno <- Geno[id.genotype.match,,drop=FALSE]
	}

	## impute missing
	if(!is.null(dim(Geno)))
	{
		if(dim(Geno)[2]>0)
		{
			if(geno_missing_imputation=="mean")
			{
				Geno <- matrix_flip_mean(Geno)$Geno
			}
			if(geno_missing_imputation=="minor")
			{
				Geno <- matrix_flip_minor(Geno)$Geno
			}
		}
	}

	## Annotation
	Anno.Int.PHRED.sub <- NULL
	Anno.Int.PHRED.sub.name <- NULL

	if(!is.null(Annotation_name))
	{
		if(variant_type=="SNV")
		{
			if(Use_annotation_weights)
			{
				for(k in 1:length(Annotation_name))
				{
					if(Annotation_name[k]%in%Annotation_name_catalog$name)
					{
						Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,Annotation_name[k])
						Annotation.PHRED <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name==Annotation_name[k])]))

						if(Annotation_name[k]=="CADD")
						{
							Annotation.PHRED[is.na(Annotation.PHRED)] <- 0
						}

						if(Annotation_name[k]=="aPC.LocalDiversity")
						{
							Annotation.PHRED.2 <- -10*log10(1-10^(-Annotation.PHRED/10))
							Annotation.PHRED <- cbind(Annotation.PHRED,Annotation.PHRED.2)
							Anno.Int.PHRED.sub.name <- c(Anno.Int.PHRED.sub.name,paste0(Annotation_name[k],"(-)"))
						}
						Anno.Int.PHRED.sub <- cbind(Anno.Int.PHRED.sub,Annotation.PHRED)
					}
				}

				Anno.Int.PHRED.sub <- data.frame(Anno.Int.PHRED.sub)
				colnames(Anno.Int.PHRED.sub) <- Anno.Int.PHRED.sub.name
			}
		}
	}else{
		Anno.Int.PHRED.sub <- NULL
	}

	pvalues <- 0
	if(n_pheno == 1)
	{
		if(!use_SPA)
		{
			if(use_ancestry_informed == FALSE){
				try(pvalues <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
			}else{
				try(pvalues <- AI_STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,find_weight=find_weight),silent=silent)
			}
		}else{
			if(use_ancestry_informed == FALSE){
				try(pvalues <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
			}else{
				try(pvalues <- AI_STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,find_weight=find_weight),silent=silent)
			}
		}
	}else
	{
		try(pvalues <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
	}

	results <- c()
	if(inherits(pvalues, "list"))
	{
		results_temp <- rep(NA,4)
		results_temp[3] <- "ncRNA"
		results_temp[2] <- chr
		results_temp[1] <- as.character(gene_name)
		results_temp[4] <- pvalues$num_variant

		if(!use_SPA)
		{
			results_temp <- c(results_temp,pvalues$cMAC,pvalues$results_STAAR_S_1_25,pvalues$results_STAAR_S_1_1,
			                  pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_A_1_25,
			                  pvalues$results_STAAR_A_1_1,pvalues$results_ACAT_O,pvalues$results_STAAR_O)
		}else
		{
			results_temp <- c(results_temp,pvalues$cMAC,
			                  pvalues$results_STAAR_B_1_25,pvalues$results_STAAR_B_1_1,pvalues$results_STAAR_B)
		}

		results <- rbind(results,results_temp)
	}

	if(!is.null(results))
	{
		if(!use_SPA)
		{
			colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
			colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
			colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("ACAT-O","STAAR-O")
		}else
		{
			colnames(results) <- colnames(results, do.NULL = FALSE, prefix = "col")
			colnames(results)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
			colnames(results)[dim(results)[2]] <- c("STAAR-B")
		}

		if(use_ancestry_informed == TRUE & find_weight == TRUE & !use_SPA){
			results_weight <- results_weight1 <- results_weight2 <- c()
			for(i in 1:ncol(pvalues$results_weight)){
				results_weight_temp <- pvalues$results_weight[-c(1,2),i]
				results_weight_temp <- unlist(pvalues$results_weight[,i][c(5:length(pvalues$results_weight[,i]), 4,3)])
				names(results_weight_temp) <- colnames(results)[-c(1:5)]

				results_weight <- cbind(results_weight, results_weight_temp)
				colnames(results_weight)[i] <- c(i-1)
			}

			for(i in 1:ncol(pvalues$results_weight1)){
				results_weight_temp1 <- pvalues$results_weight1[-c(1,2),i]
				results_weight_temp1 <- unlist(pvalues$results_weight1[,i][c(5:length(pvalues$results_weight1[,i]), 4,3)])
				names(results_weight_temp1) <- colnames(results)[-c(1:5)]

				results_weight1 <- cbind(results_weight1, results_weight_temp1)
				colnames(results_weight1)[i] <- c(i-1)
			}

			for(i in 1:ncol(pvalues$results_weight2)){
				results_weight_temp2 <- pvalues$results_weight2[-c(1,2),i]
				results_weight_temp2 <- unlist(pvalues$results_weight2[,i][c(5:length(pvalues$results_weight2[,i]), 4,3)])
				names(results_weight_temp2) <- colnames(results)[-c(1:5)]

				results_weight2 <- cbind(results_weight2, results_weight_temp2)
				colnames(results_weight2)[i] <- c(i-1)
			}
			rownames(pvalues$weight_all_1) <- rownames(pvalues$weight_all_2) <- unique(obj_nullmodel$pop.groups)

			results <- list(results,
			                weight_all_1 = pvalues$weight_all_1,
			                weight_all_2 = pvalues$weight_all_2,
			                results_weight = results_weight,
			                results_weight1 = results_weight1,
			                results_weight2 = results_weight2)
		}else if(use_ancestry_informed == TRUE & find_weight == TRUE & use_SPA){
			results_weight <- results_weight1 <- results_weight2 <- c()
			for(i in 1:ncol(pvalues$results_weight)){
				results_weight_temp <- unlist(pvalues$results_weight[,i][c(4:length(pvalues$results_weight[,i]), 3)])
				names(results_weight_temp) <- colnames(results)[-c(1:5)]

				results_weight <- cbind(results_weight, results_weight_temp)
				colnames(results_weight)[i] <- c(i-1)
			}

			for(i in 1:ncol(pvalues$results_weight1)){
				results_weight_temp1 <- unlist(pvalues$results_weight1[,i][c(4:length(pvalues$results_weight[,i]), 3)])
				names(results_weight_temp1) <- colnames(results)[-c(1:5)]

				results_weight1 <- cbind(results_weight1, results_weight_temp1)
				colnames(results_weight1)[i] <- c(i-1)
			}

			for(i in 1:ncol(pvalues$results_weight2)){
				results_weight_temp2 <- unlist(pvalues$results_weight2[,i][c(4:length(pvalues$results_weight[,i]), 3)])
				names(results_weight_temp2) <- colnames(results)[-c(1:5)]

				results_weight2 <- cbind(results_weight2, results_weight_temp2)
				colnames(results_weight2)[i] <- c(i-1)
			}
			rownames(pvalues$weight_all_1) <- rownames(pvalues$weight_all_2) <- unique(obj_nullmodel$pop.groups)

			results <- list(results,
							weight_all_1 = pvalues$weight_all_1,
							weight_all_2 = pvalues$weight_all_2,
							results_weight = results_weight,
							results_weight1 = results_weight1,
							results_weight2 = results_weight2)
		}
	}

	seqResetFilter(genofile)

	return(results)
}
