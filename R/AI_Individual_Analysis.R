#' Ancestry-informed individual-variant analysis using score test
#'
#' The \code{AI_Individual_Analysis} function takes in chromosome, an user-defined variant list,
#' the object of opened annotated GDS file, and the object from fitting the null model to analyze the association between a
#' quantitative/dichotomous phenotype and each individual variant by using score test.
#' The results of the ancestry-informed analysis correspond to ensemble p-values across base tests,
#' with the option to return a list of base weights and p-values for each base test.
#' @param chr chromosome.
#' @param individual_results the data frame of (significant) individual variants of interest for ancestry-informed analysis.
#' The first 4 columns should correspond to chromosome (CHR), position (POS), reference allele (REF), and alternative allele (ALT).
#' @param start_loc starting location (position) of the genetic region for each individual variant to be analyzed using score test.
#' @param end_loc ending location (position) of the genetic region for each individual variant to be analyzed using score test.
#' @param genofile an object of opened annotated GDS (aGDS) file.
#' @param obj_nullmodel an object from fitting the null model, which is either the output from \code{\link{fit_nullmodel}} function with two or more specified ancestries in \code{pop.groups},
#' or the output from \code{\link{fit_nullmodel}} function transformed using the \code{\link{staar2aistaar_nullmodel}} function.
#' @param QC_label channel name of the QC label in the GDS/aGDS file (default = "annotation/filter").
#' @param variant_type type of variant included in the analysis. Choices include "variant", "SNV", or "Indel" (default = "variant").
#' @param geno_missing_imputation method of handling missing genotypes. Either "mean" or "minor" (default = "mean").
#' @param tol a positive number specifying tolerance, the difference threshold for parameter
#' estimates in saddlepoint approximation algorithm below which iterations should be stopped (default = ".Machine$double.eps^0.25").
#' @param max_iter a positive integer specifying the maximum number of iterations for applying the saddlepoint approximation algorithm (default = "1000").
#' @param SPA_p_filter logical: are only the variants with a score-test-based p-value smaller than a pre-specified threshold use the SPA method to recalculate the p-value, only used for imbalanced case-control setting (default = TRUE).
#' @param p_filter_cutoff threshold for the p-value recalculation using the SPA method, only used for imbalanced case-control setting (default = 0.05)
#' @param find_weight logical: should the ancestry group-specific weights and weighting scenario-specific p-values for each base test be saved as output (default = FALSE).
#' @return A data frame containing the score test p-value and the estimated effect size of the minor allele for each individual variant in the given genetic region.
#' The first 4 columns correspond to chromosome (CHR), position (POS), reference allele (REF), and alternative allele (ALT).
#' If \code{find_weight} is TRUE, returns a list containing the ancestry-informed score test p-values, estimated effect sizes with corresponding variant characteristics,
#' as well as the ensemble weights under two sampling scenarios and p-values under scenarios 1, 2, and combined for each base test.
#' @references Chen, H., et al. (2016). Control for population structure and relatedness for binary traits
#' in genetic association studies via logistic mixed models. \emph{The American Journal of Human Genetics}, \emph{98}(4), 653-666.
#' (\href{https://doi.org/10.1016/j.ajhg.2016.02.012}{pub})
#' @references Li, Z., Li, X., et al. (2022). A framework for detecting
#' noncoding rare-variant associations of large-scale whole-genome sequencing
#' studies. \emph{Nature Methods}, \emph{19}(12), 1599-1611.
#' (\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})
#' @export

AI_Individual_Analysis <- function(chr,individual_results, start_loc, end_loc, genofile, obj_nullmodel, QC_label="annotation/filter",
                                   variant_type=c("variant","SNV","Indel"),geno_missing_imputation=c("mean","minor"),
                                   tol=.Machine$double.eps^0.25,max_iter=1000,SPA_p_filter=TRUE,p_filter_cutoff=0.05,
                                   find_weight = TRUE){

	## evaluate choices
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	individual_results_chr <- individual_results[individual_results$CHR == chr, c("CHR", "POS", "REF", "ALT")]

	## Null Model
	phenotype.id <- as.character(obj_nullmodel$id_include)
	samplesize <- length(phenotype.id)

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

	####### Obtain Genotype Information from Genofiles #######
	genotype <- char <- c()

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
	is.in <- SNVlist
	SNV.id <- variant.id[is.in]

	seqSetFilter(genofile,variant.id=SNV.id,sample.id=phenotype.id)
	position <- as.numeric(seqGetData(genofile, "position"))

	#further subset by position, which may not be unique - uniqueness further identified in matching step
	if(!is.null(start_loc) | !is.null(end_loc)){
	    if(is.null(start_loc)){
	      start_loc <- range(position)[1]
	    }
	    if(is.null(end_loc)){
	      end_loc <-  range(position)[2]
	    }
	    range_of_interest <- position[(position >= start_loc) & (position <=
	                                                     end_loc)]
	    is.in <-  range_of_interest %in% individual_results_chr$POS
	    SNV.id <- SNV.id[(position >= start_loc) & (position <= end_loc)]
  }

	seqSetFilter(genofile,variant.id=SNV.id[which(is.in)],sample.id=phenotype.id)
	CHR <- as.numeric(seqGetData(genofile, "chromosome"))
	POS <- as.numeric(seqGetData(genofile, "position"))
	REF <- as.character(seqGetData(genofile, "$ref"))
	ALT <- as.character(seqGetData(genofile, "$alt"))
	N <- rep(samplesize,length(CHR))

	## all variant identifying information from genofile
	ref_group <- data.frame(CHR=CHR,POS=POS,REF=REF,ALT=ALT)
	ref_group$id <- rownames(ref_group)

	## match variant information in provided data to those in genofile
	individual_results_chr <- dplyr::inner_join(individual_results_chr, ref_group, by = c("CHR", "POS", "REF", "ALT"))

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

	## subset variants of interest on indices unique to variants in genofile
	index <- as.numeric(individual_results_chr$id)
	Geno_chr <- Geno$Geno[,index]

	CHR_chr <- CHR[index]
	position_chr <- POS[index]
	REF_chr <- REF[index]
	ALT_chr <- ALT[index]
	N_chr <- N[index]
	Geno_chr <- as.matrix(Geno_chr,ncol=1)

	genotype <- cbind(genotype, Geno_chr)
	char <- rbind(char, cbind(CHR_chr, position_chr, REF_chr, ALT_chr, N_chr))

	####### AI-Individual Analysis #######
	genotype_ref <- genotype
	B <- dim(obj_nullmodel$pop_weights_1_1)[2]

	n_pop <- length(unique(obj_nullmodel$pop.groups))
	pop <- obj_nullmodel$pop.groups
	indices <- list()
	a_p <- matrix(0, nrow = ncol(genotype), ncol = n_pop)

	for(i in 1:n_pop)
	{
		eth <-  unique(pop)[i]
		indices[[i]] <- which(pop %in% eth)
		a_p[,i] <- apply(as.matrix(genotype[indices[[i]],]), 2, function(x){min(mean(x)/2, 1-mean(x)/2)})
	}

	a_p <- ifelse(a_p > 0, dbeta(a_p,1,25), a_p)

	w_b_1 <- w_b_2 <- matrix(0, nrow = ncol(genotype), ncol = n_pop)
	weight_all_1 <- weight_all_2 <- array(NA, dim = c(n_pop,B,ncol(genotype)))
	pvalues_1_all <- pvalues_2_all <- pvalues_12_all <- c()
	pvalues_aggregate_all <- pvalues_aggregate_s1_all <- pvalues_aggregate_s2_all <- rep(NA, ncol(genotype))

	for(g in 1:ncol(genotype))
	{
		pvalues_1_tot <- pvalues_2_tot <- matrix(NA,nrow = ncol(genotype), ncol = B)
		genotype_1_all <- genotype_2_all <- vector("list", B)

		for(b in 1:B)
		{
			w_b_1 <- matrix(rep(obj_nullmodel$pop_weights_1_1[,b], ncol(genotype)),nrow = ncol(genotype),
							ncol = n_pop, byrow = TRUE)[g,]
			w_b_2 <- (a_p%*%diag(obj_nullmodel$pop_weights_1_25[,b]))[g,]

			if(find_weight == T)
			{
				weight_all_1[,b,g] <- w_b_1
				weight_all_2[,b,g] <- w_b_2
			}

			genotype_1 <- genotype_2 <- matrix(0, nrow = nrow(genotype), ncol = 1)
			w_vec_1 <- w_vec_2 <- rep(NA,nrow(genotype))

			for(i in 1:n_pop)
			{
				eth <- unique(pop)[i]
				eth_wt_1 <- w_b_1[i]
				eth_wt_2 <- w_b_2[i]

				w_vec_1[indices[[i]]] <- eth_wt_1
        w_vec_2[indices[[i]]] <- eth_wt_2

				genotype_1[indices[[i]],] <- t(t(genotype[indices[[i]],g])*as.vector(eth_wt_1))
				genotype_2[indices[[i]],] <- t(t(genotype[indices[[i]],g])*as.vector(eth_wt_2))
			}

			genotype_1_all[[b]] <- genotype_1
			genotype_2_all[[b]] <- genotype_2
		}

		genotype_1_all <- do.call(cbind, genotype_1_all)
		genotype_2_all <- do.call(cbind, genotype_2_all)

		if(use_SPA){
		  genotype_1_all <-  as.matrix((Diagonal(x = w_vec_1) %*% genotype_1_all) - (Diagonal(x = w_vec_1) %*% XXWX_inv %*% (XW %*% genotype_1_all)))
		  genotype_2_all <-  as.matrix((Diagonal(x = w_vec_2) %*% genotype_2_all) - (Diagonal(x = w_vec_2) %*% XXWX_inv %*% (XW %*% genotype_2_all)))
		}

		if((use_SPA)&!SPA_p_filter)
		{
			pvalues_1_tot <- t(Individual_Score_Test_SPA_wt(genotype_1_all, residuals.phenotype, muhat, tol, max_iter))
			pvalues_2_tot <- t(Individual_Score_Test_SPA_wt(genotype_2_all, residuals.phenotype, muhat, tol, max_iter))
		}else
		{
			if(obj_nullmodel$sparse_kins)
			{
				Score_test1 <- Individual_Score_Test(genotype_1_all, Sigma_i, Sigma_iX, cov, residuals.phenotype)
				Score_test2 <- Individual_Score_Test(genotype_2_all, Sigma_i, Sigma_iX, cov, residuals.phenotype)
			}else if(!obj_nullmodel$sparse_kins)
			{
				Score_test1 <- Individual_Score_Test_denseGRM(genotype_1_all, P, residuals.phenotype)
				Score_test2 <- Individual_Score_Test_denseGRM(genotype_2_all, P, residuals.phenotype)
			}

		  pvalues_1_tot <- t(exp(-Score_test1$pvalue_log))
		  pvalues_2_tot <- t(exp(-Score_test2$pvalue_log))

			if(use_SPA)
			{
				if(sum(pvalues_1_tot < p_filter_cutoff)>=1)
				{
					genotype_1_all_SPA <- genotype_1_all[,pvalues_1_tot < p_filter_cutoff,drop=FALSE]
					pvalue1_SPA <- Individual_Score_Test_SPA_wt(genotype_1_all_SPA, residuals.phenotype, muhat, tol, max_iter)

					pvalues_1_tot[pvalues_1_tot < p_filter_cutoff] <- pvalue1_SPA
				}
				if(sum(pvalues_2_tot< p_filter_cutoff)>=1)
				{
					genotype_2_all_SPA <- genotype_2_all[,pvalues_2_tot < p_filter_cutoff,drop=FALSE]
					pvalue2_SPA <- Individual_Score_Test_SPA_wt(genotype_2_all_SPA, residuals.phenotype, muhat, tol, max_iter)

					pvalues_2_tot[pvalues_2_tot < p_filter_cutoff] <- pvalue2_SPA
				}
			}
		}

		obj_1 <- matrix(pvalues_1_tot,ncol = 1, nrow = B, byrow = TRUE)
		obj_2 <- matrix(pvalues_2_tot,ncol = 1, nrow = B, byrow = TRUE)

		obj_1 <- t(obj_1)
		obj_2 <- t(obj_2)

		pvalues_tot <- cbind(obj_1,obj_2)
		pvalues_aggregate_tot <- as.vector(as.numeric(formatC(pvalues_tot, format="e", digits=50)))

		if(sum(is.na(pvalues_aggregate_tot))>0){
	      ## all NAs
	      if(sum(is.na(pvalues_aggregate_tot))==length(pvalues_aggregate_tot)){
	        pvalues_aggregate <- 1
	      }else{
	        ## not all NAs
	        pvalues_aggregate_tot_sub <- na.omit(as.vector(pvalues_tot))
	        if(sum(pvalues_aggregate_tot_sub[pvalues_aggregate_tot_sub<1])>0){
	          ## not all ones
	          pvalues_aggregate <- CCT(pvalues_aggregate_tot_sub[pvalues_aggregate_tot_sub<1])
	        }else{
	          pvalues_aggregate <- 1
	        }
	      }
	    }else{
	      if(sum(pvalues_aggregate_tot[pvalues_aggregate_tot<1])>0){
	        pvalues_aggregate <- CCT(pvalues_aggregate_tot[pvalues_aggregate_tot<1])
	      }else{
	        pvalues_aggregate <- 1
	      }
    	}

		if(find_weight == T)
		{
			if(sum(is.na(obj_1))>0){
          	## all NAs
	        	if(sum(is.na(obj_1))==length(obj_1)){
	            	pvalues_aggregate_s1  <- 1
	          	}else{
	            ## not all NAs
	            	pvalues_aggregate_sub <- obj_1[!is.na(obj_1)]
	            	if(sum(pvalues_aggregate_sub[pvalues_aggregate_sub<1])>0){
	              		## not all ones
	              		pvalues_aggregate_s1  <- CCT(pvalues_aggregate_sub[pvalues_aggregate_sub<1])
	            	}else{
	              	pvalues_aggregate_s1  <- 1
	            	}
	          	}
	        }else{
	        	if(sum(obj_1[obj_1<1])>0){
	            	pvalues_aggregate_s1  <- CCT(obj_1[obj_1<1])
	          	}else{
	            	pvalues_aggregate_s1 <- 1
	          	}
        	}

        	if(sum(is.na(obj_2))>0){
          	## all NAs
	        	if(sum(is.na(obj_2))==length(obj_2)){
	            	pvalues_aggregate_s2  <- 1
	          	}else{
	            ## not all NAs
	            	pvalues_aggregate_sub <- obj_2[!is.na(obj_2)]
	            	if(sum(pvalues_aggregate_sub[pvalues_aggregate_sub<1])>0){
	              		## not all ones
	              		pvalues_aggregate_s2  <- CCT(pvalues_aggregate_sub[pvalues_aggregate_sub<1])
	            	}else{
	              	pvalues_aggregate_s2 <- 1
	            	}
	          	}
	        }else{
	        	if(sum(obj_2[obj_2<1])>0){
	            	pvalues_aggregate_s2  <- CCT(obj_2[obj_2<1])
	          	}else{
	            	pvalues_aggregate_s2 <- 1
	          	}
        	}
			pvalues_12_weight <- NULL

			for(i in 1:B)
			{
				pvalues_12_weight <-  cbind(pvalues_12_weight,CCT(pvalues_tot[c(i,i+B)]))
			}

			pvalues_1_all <- rbind(pvalues_1_all,pvalues_1_tot)
			pvalues_2_all <- rbind(pvalues_2_all,pvalues_2_tot)
			pvalues_aggregate_s1_all[g] <- pvalues_aggregate_s1
			pvalues_aggregate_s2_all[g] <- pvalues_aggregate_s2
			pvalues_12_all <- rbind(pvalues_12_all, pvalues_12_weight)
		}

		pvalues_aggregate_all[g] <- pvalues_aggregate
	}

	if(find_weight == TRUE)
	{
		dimnames(weight_all_1)[[1]] <- dimnames(weight_all_2)[[1]] <- unique(obj_nullmodel$pop.groups)
		dimnames(weight_all_1)[[2]] <- dimnames(weight_all_2)[[2]] <- seq(0,c(B-1))
		dimnames(weight_all_1)[[3]] <- dimnames(weight_all_2)[[3]] <- paste0(individual_results_chr$CHR, "_",
		                                                                     individual_results_chr$POS, "_",
		                                                                     individual_results_chr$REF,"_",
		                                                                     individual_results_chr$ALT)

		rownames(pvalues_1_all) <- rownames(pvalues_2_all) <- rownames(pvalues_12_all) <- paste0(individual_results_chr$CHR, "_",
		                                                                                         individual_results_chr$POS, "_",
		                                                                                         individual_results_chr$REF,"_",
		                                                                                         individual_results_chr$ALT)
		colnames(pvalues_1_all) <- colnames(pvalues_2_all) <- colnames(pvalues_12_all) <- seq(0,c(B-1))
	}

	if(use_SPA)
	{
		if(find_weight == TRUE)
		{
			results <- list(data.frame(CHR=CHR_chr,POS=position_chr,REF=REF_chr,ALT=ALT_chr,N=N_chr,
			                           pvalue=pvalues_aggregate_all,pvalue_s1=pvalues_aggregate_s1_all,pvalue_s2=pvalues_aggregate_s2_all,
			                           pvalue_log10=-log10(pvalues_aggregate_all)),
			                weight_all_1=weight_all_1, weight_all_2=weight_all_2, results_weight=pvalues_12_all,
			                results_weight1=pvalues_1_all, results_weight2=pvalues_2_all)
		}else
		{
			results <- data.frame(CHR=CHR_chr,POS=position_chr,REF=REF_chr,ALT=ALT_chr,N=N_chr,
								  pvalue=pvalues_aggregate_all,pvalue_log10=-log10(pvalues_aggregate_all))
		}
	}else
	{
		if(find_weight == TRUE)
		{
			results <- list(data.frame(CHR=CHR_chr,POS=position_chr,REF=REF_chr,ALT=ALT_chr,N=N_chr,
			                           pvalue=pvalues_aggregate_all,pvalue_s1=pvalues_aggregate_s1_all,pvalue_s2=pvalues_aggregate_s2_all,
			                           pvalue_log10=-log10(pvalues_aggregate_all)),
			                weight_all_1=weight_all_1, weight_all_2=weight_all_2, results_weight=pvalues_12_all,
			                results_weight1=pvalues_1_all, results_weight2=pvalues_2_all)
		}else
		{
			results <- data.frame(CHR=CHR_chr,POS=position_chr,REF=REF_chr,ALT=ALT_chr,N=N_chr,
			                      pvalue=pvalues_aggregate_all,pvalue_log10=-log10(pvalues_aggregate_all))
		}
	}

	seqResetFilter(genofile)

	return(results)

}
