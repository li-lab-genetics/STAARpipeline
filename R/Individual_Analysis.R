#' Individual-variant analysis using score test
#'
#' The \code{Individual_Analysis} function takes in chromosome, starting location, ending location, an user-defined variant list for
#' ancestry-informed analyses, the object of opened annotated GDS file, and the object from fitting the null model to analyze the association between a
#' quantitative/dichotomous phenotype (including imbalanced case-control design) and each individual variant in a genetic region by using score test.
#' For multiple phenotype analysis (\code{obj_nullmodel$n.pheno > 1}),
#' the results correspond to multi-trait score test p-values by leveraging
#' the correlation structure between multiple phenotypes.
#' For ancestry-informed analysis, the results correspond to ensemble p-values across base tests,
#' with the option to return a list of base weights and p-values for each base test.
#' @param chr chromosome.
#' @param start_loc starting location (position) of the genetic region for each individual variant to be analyzed using score test.
#' @param end_loc ending location (position) of the genetic region for each individual variant to be analyzed using score test.
#' @param individual_results the data frame of (significant) individual variants of interest for ancestry-informed analysis.
#' The first 4 columns should correspond to chromosome (CHR), position (POS), reference allele (REF), and alternative allele (ALT).
#' @param genofile an object of opened annotated GDS (aGDS) file.
#' @param obj_nullmodel an object from fitting the null model, which is either the output from \code{\link{fit_nullmodel}} function,
#' or the output from \code{fitNullModel} function in the \code{GENESIS} package and transformed using the \code{\link{genesis2staar_nullmodel}} function.
#' @param mac_cutoff the cutoff of minimum minor allele count in
#' defining individual variants (default = 20).
#' @param subset_variants_num the number of variants to run per subset for each time (default = 5e3).
#' @param QC_label channel name of the QC label in the GDS/aGDS file (default = "annotation/filter").
#' @param variant_type type of variant included in the analysis. Choices include "variant", "SNV", or "Indel" (default = "variant").
#' @param geno_missing_imputation method of handling missing genotypes. Either "mean" or "minor" (default = "mean").
#' @param tol a positive number specifying tolerance, the difference threshold for parameter
#' estimates in saddlepoint approximation algorithm below which iterations should be stopped (default = ".Machine$double.eps^0.25").
#' @param max_iter a positive integer specifying the maximum number of iterations for applying the saddlepoint approximation algorithm (default = "1000").
#' @param SPA_p_filter logical: are only the variants with a score-test-based p-value smaller than a pre-specified threshold use the SPA method to recalculate the p-value, only used for imbalanced case-control setting (default = TRUE).
#' @param p_filter_cutoff threshold for the p-value recalculation using the SPA method, only used for imbalanced case-control setting (default = 0.05)
#' @param use_ancestry_informed logical: is ancestry-informed association analysis used to estimate p-values (default = FALSE).
#' @param find_weight logical: should the ancestry group-specific weights and weighting scenario-specific p-values for each base test be saved as output (default = FALSE).
#' @return A data frame containing the score test p-value and the estimated effect size of the minor allele for each individual variant in the given genetic region, or as provided in \code{individual_results}
#' for ancestry-informed variant analysis. The first 4 columns correspond to chromosome (CHR), position (POS), reference allele (REF), and alternative allele (ALT).
#' If \code{find_weight} is TRUE, returns a list containing the ancestry-informed score test p-values and the estimated effect size of the minor allele for each individual variant provided in \code{individual_results}.
#' The ensemble weights under two sampling scenarios and p-values under scenarios 1, 2, and combined for each base test are saved as well.
#' @references Chen, H., et al. (2016). Control for population structure and relatedness for binary traits
#' in genetic association studies via logistic mixed models. \emph{The American Journal of Human Genetics}, \emph{98}(4), 653-666.
#' (\href{https://doi.org/10.1016/j.ajhg.2016.02.012}{pub})
#' @references Li, Z., Li, X., et al. (2022). A framework for detecting
#' noncoding rare-variant associations of large-scale whole-genome sequencing
#' studies. \emph{Nature Methods}, \emph{19}(12), 1599-1611.
#' (\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})
#' @export

Individual_Analysis <- function(chr,start_loc=NULL,end_loc=NULL,individual_results=NULL,genofile,obj_nullmodel,mac_cutoff=20,subset_variants_num=5e3,
                                QC_label="annotation/filter",variant_type=c("variant","SNV","Indel"),geno_missing_imputation=c("mean","minor"),
                                tol=.Machine$double.eps^0.25,max_iter=1000,SPA_p_filter=TRUE,p_filter_cutoff=0.05,
                                use_ancestry_informed=FALSE,find_weight=FALSE){

	## evaluate choices
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	## Null Model
	phenotype.id <- as.character(obj_nullmodel$id_include)
	samplesize <- length(phenotype.id)
	n_pheno <- obj_nullmodel$n.pheno

	if(!is.null(obj_nullmodel$use_SPA))
	{
		use_SPA <- obj_nullmodel$use_SPA
	}else
	{
		use_SPA <- FALSE
	}

	## residuals and cov
	residuals.phenotype <- as.vector(obj_nullmodel$scaled.residuals)
	if(SPA_p_filter)
	{
		### dense GRM
		if(!obj_nullmodel$sparse_kins)
		{
			P <- obj_nullmodel$P
		}

		### sparse GRM
		if(obj_nullmodel$sparse_kins)
		{
			Sigma_i <- obj_nullmodel$Sigma_i
			Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
			cov <- obj_nullmodel$cov
		}
	}

	## SPA
	if(use_SPA)
	{
		muhat <- obj_nullmodel$fitted.values

		if(obj_nullmodel$relatedness)
		{
			if(!obj_nullmodel$sparse_kins)
			{
				XW <- obj_nullmodel$XW
				XXWX_inv <- obj_nullmodel$XXWX_inv
			}else
			{
				XW <- as.matrix(obj_nullmodel$XSigma_i)
				XXWX_inv <- as.matrix(obj_nullmodel$XXSigma_iX_inv)
			}
		}else
		{
			XW <- obj_nullmodel$XW
			XXWX_inv <- obj_nullmodel$XXWX_inv
		}
	}else
	{
		### dense GRM
		if(!obj_nullmodel$sparse_kins)
		{
			P <- obj_nullmodel$P
		}

		### sparse GRM
		if(obj_nullmodel$sparse_kins)
		{
			Sigma_i <- obj_nullmodel$Sigma_i
			Sigma_iX <- as.matrix(obj_nullmodel$Sigma_iX)
			cov <- obj_nullmodel$cov
		}
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

	results <- c()

	if(use_ancestry_informed)
	{
		results <- AI_Individual_Analysis(chr=chr,start_loc=start_loc,end_loc=end_loc,individual_results=individual_results,
		                                  genofile=genofile,obj_nullmodel=obj_nullmodel,QC_label=QC_label,variant_type=variant_type,
		                                  geno_missing_imputation=geno_missing_imputation,find_weight=find_weight)
		return(results)
	}

	position <- as.numeric(seqGetData(genofile, "position"))

	variant.id <- seqGetData(genofile, "variant.id")
	is.in <- (SNVlist)&(position>=start_loc)&(position<=end_loc)
	SNV.id <- variant.id[is.in]

	subset.num <- ceiling(length(SNV.id)/subset_variants_num)

	if(subset.num == 0)
	{
		return(results)
	}

	for(kk in 1:subset.num)
	{
		if(kk < subset.num)
		{
			is.in <- ((kk-1)*subset_variants_num+1):(kk*subset_variants_num)
			seqSetFilter(genofile,variant.id=SNV.id[is.in],sample.id=phenotype.id)
		}
		if(kk == subset.num)
		{
			is.in <- ((kk-1)*subset_variants_num+1):length(SNV.id)
			seqSetFilter(genofile,variant.id=SNV.id[is.in],sample.id=phenotype.id)
		}

		## genotype id
		id.genotype <- seqGetData(genofile,"sample.id")

		id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
		phenotype.id.merge <- data.frame(phenotype.id)
		phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("phenotype.id"="id.genotype"))
		id.genotype.match <- phenotype.id.merge$index


		Geno <- seqGetData(genofile, "$dosage")
		Geno <- Geno[id.genotype.match,,drop=FALSE]

		if(geno_missing_imputation=="mean")
		{
			Geno <- matrix_flip_mean(Geno)
		}
		if(geno_missing_imputation=="minor")
		{
			Geno <- matrix_flip_minor(Geno)
		}

		MAF <- Geno$MAF
		ALT_AF <- 1 - Geno$AF

		CHR <- as.numeric(seqGetData(genofile, "chromosome"))
		position <- as.numeric(seqGetData(genofile, "position"))
		REF <- as.character(seqGetData(genofile, "$ref"))
		ALT <- as.character(seqGetData(genofile, "$alt"))
		N <- rep(samplesize,length(CHR))

		if(!all(CHR==chr))
		{
			warning("chr does not match the chromosome of genofile (the opened aGDS)!")
		}

		if((use_SPA)&!SPA_p_filter)
		{
			if(sum(MAF>(mac_cutoff-0.5)/samplesize/2)>=1)
			{
				Geno <- Geno$Geno

				## Common_variants
				Geno_common <- Geno[,(MAF>(mac_cutoff-0.5)/samplesize/2),drop=FALSE]

				CHR_common <- CHR[(MAF>(mac_cutoff-0.5)/samplesize/2)]
				position_common <- position[(MAF>(mac_cutoff-0.5)/samplesize/2)]
				REF_common <- REF[(MAF>(mac_cutoff-0.5)/samplesize/2)]
				ALT_common <- ALT[(MAF>(mac_cutoff-0.5)/samplesize/2)]
				MAF_common <- MAF[(MAF>(mac_cutoff-0.5)/samplesize/2)]
				ALT_AF_common <- ALT_AF[(MAF>(mac_cutoff-0.5)/samplesize/2)]
				N_common <- N[(MAF>(mac_cutoff-0.5)/samplesize/2)]

				rm(Geno)
				gc()

				pvalue <- Individual_Score_Test_SPA(Geno_common,XW,XXWX_inv,residuals.phenotype,muhat,tol,max_iter)

				results_temp <- data.frame(CHR=CHR_common,POS=position_common,REF=REF_common,ALT=ALT_common,ALT_AF=ALT_AF_common,MAF=MAF_common,N=N_common,
				                           pvalue=pvalue)

				results <- rbind(results,results_temp)
			}
		}else
		{
			## Common_variants
			if(sum(MAF>=0.05)>=1)
			{
				Geno_common <- Geno$Geno[,MAF>=0.05]

				CHR_common <- CHR[MAF>=0.05]
				position_common <- position[MAF>=0.05]
				REF_common <- REF[MAF>=0.05]
				ALT_common <- ALT[MAF>=0.05]
				MAF_common <- MAF[MAF>=0.05]
				ALT_AF_common <- ALT_AF[MAF>=0.05]
				N_common <- N[MAF>=0.05]

				if(sum(MAF>=0.05)==1)
				{
					Geno_common <- as.matrix(Geno_common,ncol=1)
				}

				## sparse GRM
				if(obj_nullmodel$sparse_kins)
				{
					if(n_pheno == 1)
					{
						Score_test <- Individual_Score_Test(Geno_common, Sigma_i, Sigma_iX, cov, residuals.phenotype)
					}else
					{
						Geno_common <- Diagonal(n = n_pheno) %x% Geno_common
						Score_test <- Individual_Score_Test_sp_multi(Geno_common, Sigma_i, Sigma_iX, cov, residuals.phenotype, n_pheno)
					}
				}

				## dense GRM
				if(!obj_nullmodel$sparse_kins)
				{
					if(n_pheno == 1)
					{
						Score_test <- Individual_Score_Test_denseGRM(Geno_common, P, residuals.phenotype)
					}else
					{
						Geno_common <- Diagonal(n = n_pheno) %x% Geno_common
						Score_test <- Individual_Score_Test_sp_denseGRM_multi(Geno_common, P, residuals.phenotype, n_pheno)
					}
				}

				## SPA approximation for small p-values
				if(use_SPA)
				{
					pvalue <- exp(-Score_test$pvalue_log)

					if(sum(pvalue < p_filter_cutoff)>=1)
					{
						Geno_common_SPA <- Geno_common[,pvalue < p_filter_cutoff,drop=FALSE]
						pvalue_SPA <- Individual_Score_Test_SPA(Geno_common_SPA,XW,XXWX_inv,residuals.phenotype,muhat,tol,max_iter)

						pvalue[pvalue < p_filter_cutoff] <- pvalue_SPA
					}
				}

				if(use_SPA)
				{
					results_temp <- data.frame(CHR=CHR_common,POS=position_common,REF=REF_common,ALT=ALT_common,ALT_AF=ALT_AF_common,MAF=MAF_common,N=N_common,
					                           pvalue=pvalue)
				}else
				{
					if(n_pheno == 1)
					{
						results_temp <- data.frame(CHR=CHR_common,POS=position_common,REF=REF_common,ALT=ALT_common,ALT_AF=ALT_AF_common,MAF=MAF_common,N=N_common,
						                           pvalue=exp(-Score_test$pvalue_log),pvalue_log10=Score_test$pvalue_log/log(10),
						                           Score=Score_test$Score,Score_se=Score_test$Score_se,
						                           Est=Score_test$Est,Est_se=Score_test$Est_se)
					}else
					{
						results_temp <- data.frame(CHR=CHR_common,POS=position_common,REF=REF_common,ALT=ALT_common,ALT_AF=ALT_AF_common,MAF=MAF_common,N=N_common,
						                           pvalue=exp(-Score_test$pvalue_log),pvalue_log10=Score_test$pvalue_log/log(10))
						results_temp <- cbind(results_temp,matrix(Score_test$Score,ncol=n_pheno))
						colnames(results_temp)[10:(10+n_pheno-1)] <- paste0("Score",seq_len(n_pheno))
					}
				}
				results <- rbind(results,results_temp)
			}


			## Rare_variants
			if(sum((MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05))>=1)
			{
				Geno_rare <- Geno$Geno[,(MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05)]

				CHR_rare <- CHR[(MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05)]
				position_rare <- position[(MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05)]
				REF_rare <- REF[(MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05)]
				ALT_rare <- ALT[(MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05)]
				MAF_rare <- MAF[(MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05)]
				ALT_AF_rare <- ALT_AF[(MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05)]
				N_rare <- N[(MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05)]

				## sparse GRM
				if(obj_nullmodel$sparse_kins)
				{
					if(sum((MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05))>=2)
					{
						if(n_pheno == 1)
						{
							Geno_rare <- as(Geno_rare,"dgCMatrix")
							Score_test <- Individual_Score_Test_sp(Geno_rare, Sigma_i, Sigma_iX, cov, residuals.phenotype)
						}else
						{
							Geno_rare <- Diagonal(n = n_pheno) %x% Geno_rare
							Score_test <- Individual_Score_Test_sp_multi(Geno_rare, Sigma_i, Sigma_iX, cov, residuals.phenotype, n_pheno)
						}
					}else
					{
						if(n_pheno == 1)
						{
							Geno_rare <- as.matrix(Geno_rare,ncol=1)
							Score_test <- Individual_Score_Test(Geno_rare, Sigma_i, Sigma_iX, cov, residuals.phenotype)
						}
						else
						{
							Geno_rare <- as.matrix(Diagonal(n = n_pheno) %x% Geno_rare)
							Score_test <- Individual_Score_Test_multi(Geno_rare, Sigma_i, Sigma_iX, cov, residuals.phenotype, n_pheno)
						}
					}
				}

				## dense GRM
				if(!obj_nullmodel$sparse_kins)
				{
					if(sum((MAF>(mac_cutoff-0.5)/samplesize/2)&(MAF<0.05))>=2)
					{
						if(n_pheno == 1)
						{
							Geno_rare <- as(Geno_rare,"dgCMatrix")
							Score_test <- Individual_Score_Test_sp_denseGRM(Geno_rare, P, residuals.phenotype)
						}
						else{
							Geno_rare <- Diagonal(n = n_pheno) %x% Geno_rare
							Score_test <- Individual_Score_Test_sp_denseGRM_multi(Geno_rare, P, residuals.phenotype, n_pheno)
						}
					}else
					{
						if(n_pheno == 1)
						{
							Geno_rare <- as.matrix(Geno_rare,ncol=1)
							Score_test <- Individual_Score_Test_denseGRM(Geno_rare, P, residuals.phenotype)
						}
						else
						{
							Geno_rare <- as.matrix(Diagonal(n = n_pheno) %x% Geno_rare)
							Score_test <- Individual_Score_Test_denseGRM_multi(Geno_rare, P, residuals.phenotype, n_pheno)
						}
					}
				}

				## SPA approximation for small p-values
				if(use_SPA)
				{
					pvalue <- exp(-Score_test$pvalue_log)

					if(sum(pvalue < p_filter_cutoff)>=2)
					{
						Geno_rare_SPA <- as.matrix(Geno_rare)[,pvalue < p_filter_cutoff]
					}

					if(sum(pvalue < p_filter_cutoff)==1)
					{
						Geno_rare_SPA <- as.matrix(Geno_rare)[,pvalue < p_filter_cutoff]
						Geno_rare_SPA <- as.matrix(Geno_rare_SPA,ncol=1)
					}

					if(sum(pvalue < p_filter_cutoff)>=1)
					{
						pvalue_SPA <- Individual_Score_Test_SPA(Geno_rare_SPA,XW,XXWX_inv,residuals.phenotype,muhat,tol,max_iter)

						pvalue[pvalue < p_filter_cutoff] <- pvalue_SPA
					}

				}

				if(use_SPA)
				{
					results_temp <- data.frame(CHR=CHR_rare,POS=position_rare,REF=REF_rare,ALT=ALT_rare,ALT_AF=ALT_AF_rare,MAF=MAF_rare,N=N_rare,
					                           pvalue=pvalue)

				}else
				{
					if(n_pheno == 1)
					{
						results_temp <- data.frame(CHR=CHR_rare,POS=position_rare,REF=REF_rare,ALT=ALT_rare,ALT_AF=ALT_AF_rare,MAF=MAF_rare,N=N_rare,
						                           pvalue=exp(-Score_test$pvalue_log),pvalue_log10=Score_test$pvalue_log/log(10),
						                           Score=Score_test$Score,Score_se=Score_test$Score_se,
						                           Est=Score_test$Est,Est_se=Score_test$Est_se)
					}
					else
					{
						results_temp <- data.frame(CHR=CHR_rare,POS=position_rare,REF=REF_rare,ALT=ALT_rare,ALT_AF=ALT_AF_rare,MAF=MAF_rare,N=N_rare,
						                           pvalue=exp(-Score_test$pvalue_log),pvalue_log10=Score_test$pvalue_log/log(10))
						results_temp <- cbind(results_temp,matrix(Score_test$Score,ncol=n_pheno))
						colnames(results_temp)[10:(10+n_pheno-1)] <- paste0("Score",seq_len(n_pheno))
					}
				}

				results <- rbind(results,results_temp)
			}
		}
		seqResetFilter(genofile)
	}

	if(!is.null(results))
	{
		results <- results[order(results[,2]),]
	}

	return(results)
}

