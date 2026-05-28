missense <- function(chr,gene_name,genofile,obj_nullmodel,genes,
                     rare_maf_cutoff=0.01,rv_num_cutoff=2,rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                     QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                     Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                     Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                     SPA_p_filter=FALSE,p_filter_cutoff=0.05,use_ancestry_informed=FALSE,find_weight=FALSE,silent=FALSE){

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

	## get SNV id, position, REF, ALT (whole genome)
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

	position <- as.numeric(seqGetData(genofile, "position"))
	variant.id <- seqGetData(genofile, "variant.id")

	### Gene
	kk <- which(genes[,1]==gene_name)

	sub_start_loc <- genes[kk,3]
	sub_end_loc <- genes[kk,4]

	############################################################
	#                      missense
	############################################################

	is.in <- (SNVlist)&(position>=sub_start_loc)&(position<=sub_end_loc)
	variant.id.gene <- variant.id[is.in]

	seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

	## Gencode_Exonic
	GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))

	variant.id.gene <- seqGetData(genofile, "variant.id")
	lof.in.missense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")
	variant.id.gene <- variant.id.gene[lof.in.missense]

	seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

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
		results_temp <- as.vector(genes[kk,])
		results_temp[3] <- "missense"
		results_temp[2] <- chr
		results_temp[1] <- as.character(genes[kk,1])
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

	############################################################
	#                      disruptive missense
	############################################################

	is.in <- (SNVlist)&(position>=sub_start_loc)&(position<=sub_end_loc)
	variant.id.gene <- variant.id[is.in]

	seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

	## Gencode_Exonic
	GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
	## Meta.SVM.Pred
	MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))

	variant.id.gene <- seqGetData(genofile, "variant.id")
	lof.in.dmissense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D")
	variant.id.gene <- variant.id.gene[lof.in.dmissense]

	seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

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

	pvalues_dm <- 0
	if(n_pheno == 1)
	{
		if(!use_SPA)
		{
			if(use_ancestry_informed == FALSE){
				try(pvalues_dm <- STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
			}else{
				try(pvalues_dm <- AI_STAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,find_weight=find_weight),silent=silent)
			}
		}else{
			if(use_ancestry_informed == FALSE){
				try(pvalues_dm <- STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
			}else{
				try(pvalues_dm <- AI_STAAR_Binary_SPA(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,find_weight=find_weight),silent=silent)
			}
		}
	}else
	{
		try(pvalues_dm <- MultiSTAAR(Geno,obj_nullmodel,Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
	}

	if(inherits(pvalues_dm, "list"))
	{
		results_temp <- as.vector(genes[kk,])
		results_temp[3] <- "disruptive_missense"
		results_temp[2] <- chr
		results_temp[1] <- as.character(genes[kk,1])
		results_temp[4] <- pvalues_dm$num_variant

		if(!use_SPA)
		{
			results_temp <- c(results_temp,pvalues_dm$cMAC,pvalues_dm$results_STAAR_S_1_25,pvalues_dm$results_STAAR_S_1_1,
			                  pvalues_dm$results_STAAR_B_1_25,pvalues_dm$results_STAAR_B_1_1,pvalues_dm$results_STAAR_A_1_25,
			                  pvalues_dm$results_STAAR_A_1_1,pvalues_dm$results_ACAT_O,pvalues_dm$results_STAAR_O)
		}else
		{
			results_temp <- c(results_temp,pvalues_dm$cMAC,
			                  pvalues_dm$results_STAAR_B_1_25,pvalues_dm$results_STAAR_B_1_1,pvalues_dm$results_STAAR_B)
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

		if(dim(results)[1]==1)
		{
			if(results[3]!="disruptive_missense")
			{
				if(!use_SPA)
				{
					results <- cbind(results,matrix(1,1,6))
					colnames(results)[(dim(results)[2]-5):dim(results)[2]] <- c("SKAT(1,25)-Disruptive","SKAT(1,1)-Disruptive","Burden(1,25)-Disruptive","Burden(1,1)-Disruptive","ACAT-V(1,25)-Disruptive","ACAT-V(1,1)-Disruptive")
				}else{
					results <- cbind(results,matrix(1,1,2))
					colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("Burden(1,25)-Disruptive","Burden(1,1)-Disruptive")
				}

				if(use_ancestry_informed == TRUE & find_weight == TRUE & !use_SPA)
				{
					results_weight <- results_weight1 <- results_weight2 <- c()
					for(i in 1:ncol(pvalues$results_weight)){
						results_m_weight <- unlist(pvalues$results_weight[,i][c(5:length(pvalues$results_weight[,i]), 4,3)])
						names(results_m_weight) <- colnames(results)[-c(1:5,(dim(results)[2]-5):dim(results)[2])]
						
						results_m_weight <- c(results_m_weight, rep(1,6))
						names(results_m_weight)[(length(results_m_weight)-5):length(results_m_weight)] <- c("SKAT(1,25)-Disruptive","SKAT(1,1)-Disruptive","Burden(1,25)-Disruptive","Burden(1,1)-Disruptive","ACAT-V(1,25)-Disruptive","ACAT-V(1,1)-Disruptive")
						results_weight <- cbind(results_weight, results_m_weight)
						colnames(results_weight)[i] <- c(i-1)
					}

					for(i in 1:ncol(pvalues$results_weight1)){
						results_m_weight1 <- unlist(pvalues$results_weight1[,i][c(5:length(pvalues$results_weight1[,i]), 4,3)])
						names(results_m_weight1) <- colnames(results)[-c(1:5,(dim(results)[2]-5):dim(results)[2])]
						
						results_m_weight1 <- c(results_m_weight1, rep(1,6))
						names(results_m_weight1)[(length(results_m_weight1)-5):length(results_m_weight1)] <- c("SKAT(1,25)-Disruptive","SKAT(1,1)-Disruptive","Burden(1,25)-Disruptive","Burden(1,1)-Disruptive","ACAT-V(1,25)-Disruptive","ACAT-V(1,1)-Disruptive")
						results_weight1 <- cbind(results_weight1, results_m_weight1)
						colnames(results_weight1)[i] <- c(i-1)
					}

					for(i in 1:ncol(pvalues$results_weight2)){
						results_m_weight2 <- unlist(pvalues$results_weight2[,i][c(5:length(pvalues$results_weight2[,i]), 4,3)])
						names(results_m_weight2) <- colnames(results)[-c(1:5,(dim(results)[2]-5):dim(results)[2])]
						
						results_m_weight2 <- c(results_m_weight2, rep(1,6))
						names(results_m_weight2)[(length(results_m_weight2)-5):length(results_m_weight2)] <- c("SKAT(1,25)-Disruptive","SKAT(1,1)-Disruptive","Burden(1,25)-Disruptive","Burden(1,1)-Disruptive","ACAT-V(1,25)-Disruptive","ACAT-V(1,1)-Disruptive")
						results_weight2 <- cbind(results_weight2, results_m_weight2)
						colnames(results_weight2)[i] <- c(i-1)
					}

					rownames(pvalues$weight_all_1) <- rownames(pvalues$weight_all_2) <- unique(obj_nullmodel$pop.groups)
					results <- list(results,
					                weight_all_1 = pvalues$weight_all_1,
					                weight_all_2 = pvalues$weight_all_2,
					                results_weight = results_weight,
					                results_weight1 = results_weight1,
					                results_weight2 = results_weight2)

				}else if(use_ancestry_informed == TRUE & find_weight == TRUE & use_SPA)
				{
					results_weight <- results_weight1 <- results_weight2 <- c()
					for(i in 1:ncol(pvalues$results_weight)){
						results_m_weight <- unlist(pvalues$results_weight[,i][c(4:length(pvalues$results_weight[,i]), 3)])
						names(results_m_weight) <- colnames(results)[6:10]

						results_m_weight <- c(results_m_weight, rep(1,2))
						names(results_m_weight)[(length(results_m_weight)-1):length(results_m_weight)] <- c("Burden(1,25)-Disruptive", "Burden(1,1)-Disruptive")
						results_weight <- cbind(results_weight, results_m_weight)
						colnames(results_weight)[i] <- c(i-1)
					}

					for(i in 1:ncol(pvalues$results_weight1)){
						results_m_weight1 <- unlist(pvalues$results_weight1[,i][c(4:length(pvalues$results_weight1[,i]), 3)])
						names(results_m_weight1) <- colnames(results)[6:10]

						results_m_weight1 <- c(results_m_weight1, rep(1,2))
						names(results_m_weight1)[(length(results_m_weight1)-1):length(results_m_weight1)] <- c("Burden(1,25)-Disruptive", "Burden(1,1)-Disruptive")
						results_weight1 <- cbind(results_weight1, results_m_weight1)
						colnames(results_weight1)[i] <- c(i-1)
					}

					for(i in 1:ncol(pvalues$results_weight2)){
						results_m_weight2 <- unlist(pvalues$results_weight2[,i][c(4:length(pvalues$results_weight2[,i]), 3)])
						names(results_m_weight2) <- colnames(results)[6:10]

						results_m_weight2 <- c(results_m_weight2, rep(1,2))
						names(results_m_weight2)[(length(results_m_weight2)-1):length(results_m_weight2)] <- c("Burden(1,25)-Disruptive", "Burden(1,1)-Disruptive")
						results_weight2 <- cbind(results_weight2, results_m_weight2)
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
			}else
			{
				results <- c()
			}
		}else if(dim(results)[1]==2)
		{
			if(!use_SPA)
			{
				results_m <- c(results[1,],rep(0,6))
				names(results_m)[(length(results_m)-5):length(results_m)] <- c("SKAT(1,25)-Disruptive","SKAT(1,1)-Disruptive","Burden(1,25)-Disruptive","Burden(1,1)-Disruptive","ACAT-V(1,25)-Disruptive","ACAT-V(1,1)-Disruptive")
				results_m[(length(results_m)-5):length(results_m)] <- results[2,c("SKAT(1,25)","SKAT(1,1)","Burden(1,25)","Burden(1,1)","ACAT-V(1,25)","ACAT-V(1,1)")]
				apc_num <- (length(results_m)-19)/6
				p_seq <- c(1:apc_num,1:apc_num+(apc_num+1),1:apc_num+2*(apc_num+1),1:apc_num+3*(apc_num+1),1:apc_num+4*(apc_num+1),1:apc_num+5*(apc_num+1),(6*apc_num+9):(6*apc_num+14))
				results_m["STAAR-O"] <- CCT(as.numeric(results_m[6:length(results_m)][p_seq]))
				results_m["STAAR-S(1,25)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num,6*apc_num+9)]))
				results_m["STAAR-S(1,1)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+(apc_num+1),6*apc_num+10)]))
				results_m["STAAR-B(1,25)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+2*(apc_num+1),6*apc_num+11)]))
				results_m["STAAR-B(1,1)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+3*(apc_num+1),6*apc_num+12)]))
				results_m["STAAR-A(1,25)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+4*(apc_num+1),6*apc_num+13)]))
				results_m["STAAR-A(1,1)"] <- CCT(as.numeric(results_m[6:length(results_m)][c(1:apc_num+5*(apc_num+1),6*apc_num+14)]))

				results <- c()
				results <- rbind(results,results_m)

				#results_weight
				if(use_ancestry_informed == TRUE & find_weight == TRUE){
					results_weight <- results_weight1 <- results_weight2 <- c()

					for(i in 1:ncol(pvalues$results_weight)){
						results_m_weight <- unlist(pvalues$results_weight[,i][c(5:length(pvalues$results_weight[,i]), 4,3)])
						names(results_m_weight) <- names(results_m)[-c(1:5,(length(results_m)-5):length(results_m))]

						results_m_weight <- c(results_m_weight, rep(0,6))
						names(results_m_weight)[(length(results_m_weight)-5):length(results_m_weight)] <- c("SKAT(1,25)-Disruptive","SKAT(1,1)-Disruptive","Burden(1,25)-Disruptive","Burden(1,1)-Disruptive","ACAT-V(1,25)-Disruptive","ACAT-V(1,1)-Disruptive")
						results_m_weight[(length(results_m_weight)-5):length(results_m_weight)] <- unlist(pvalues_dm$results_weight[,i][c("results_STAAR_S_1_25.SKAT(1,25)","results_STAAR_S_1_1.SKAT(1,1)","results_STAAR_B_1_25.Burden(1,25)",
						                                                                                                                  "results_STAAR_B_1_1.Burden(1,1)","results_STAAR_A_1_25.ACAT-V(1,25)","results_STAAR_A_1_1.ACAT-V(1,1)")])
						results_m_weight["STAAR-O"] <- CCT(as.numeric(results_m_weight[1:length(results_m_weight)][p_seq]))
						results_m_weight["STAAR-S(1,25)"] <- CCT(as.numeric(results_m_weight[1:length(results_m_weight)][c(1:apc_num,6*apc_num+9)]))
						results_m_weight["STAAR-S(1,1)"] <- CCT(as.numeric(results_m_weight[1:length(results_m_weight)][c(1:apc_num+(apc_num+1),6*apc_num+10)]))
						results_m_weight["STAAR-B(1,25)"] <- CCT(as.numeric(results_m_weight[1:length(results_m_weight)][c(1:apc_num+2*(apc_num+1),6*apc_num+11)]))
						results_m_weight["STAAR-B(1,1)"] <- CCT(as.numeric(results_m_weight[1:length(results_m_weight)][c(1:apc_num+3*(apc_num+1),6*apc_num+12)]))
						results_m_weight["STAAR-A(1,25)"] <- CCT(as.numeric(results_m_weight[1:length(results_m_weight)][c(1:apc_num+4*(apc_num+1),6*apc_num+13)]))
						results_m_weight["STAAR-A(1,1)"] <- CCT(as.numeric(results_m_weight[1:length(results_m_weight)][c(1:apc_num+5*(apc_num+1),6*apc_num+14)]))

						results_weight <- cbind(results_weight, results_m_weight)
						colnames(results_weight)[i] <- c(i-1)
					}

					#results_weight1
					for(i in 1:ncol(pvalues$results_weight1)){
						results_m_weight1 <- unlist(pvalues$results_weight1[,i][c(5:length(pvalues$results_weight1[,i]), 4,3)])
						names(results_m_weight1) <- names(results_m)[-c(1:5,(length(results_m)-5):length(results_m))]

						results_m_weight1 <- c(results_m_weight1, rep(0,6))
						names(results_m_weight1)[(length(results_m_weight1)-5):length(results_m_weight1)] <- c("SKAT(1,25)-Disruptive","SKAT(1,1)-Disruptive","Burden(1,25)-Disruptive","Burden(1,1)-Disruptive","ACAT-V(1,25)-Disruptive","ACAT-V(1,1)-Disruptive")
						results_m_weight1[(length(results_m_weight1)-5):length(results_m_weight1)] <- unlist(pvalues_dm$results_weight1[,i][c("results_STAAR_S_1_25.SKAT(1,25)","results_STAAR_S_1_1.SKAT(1,1)","results_STAAR_B_1_25.Burden(1,25)",
						                                                                                                                      "results_STAAR_B_1_1.Burden(1,1)","results_STAAR_A_1_25.ACAT-V(1,25)","results_STAAR_A_1_1.ACAT-V(1,1)")])
						results_m_weight1["STAAR-O"] <- CCT(as.numeric(results_m_weight1[1:length(results_m_weight1)][p_seq]))
						results_m_weight1["STAAR-S(1,25)"] <- CCT(as.numeric(results_m_weight1[1:length(results_m_weight1)][c(1:apc_num,6*apc_num+9)]))
						results_m_weight1["STAAR-S(1,1)"] <- CCT(as.numeric(results_m_weight1[1:length(results_m_weight1)][c(1:apc_num+(apc_num+1),6*apc_num+10)]))
						results_m_weight1["STAAR-B(1,25)"] <- CCT(as.numeric(results_m_weight1[1:length(results_m_weight1)][c(1:apc_num+2*(apc_num+1),6*apc_num+11)]))
						results_m_weight1["STAAR-B(1,1)"] <- CCT(as.numeric(results_m_weight1[1:length(results_m_weight1)][c(1:apc_num+3*(apc_num+1),6*apc_num+12)]))
						results_m_weight1["STAAR-A(1,25)"] <- CCT(as.numeric(results_m_weight1[1:length(results_m_weight1)][c(1:apc_num+4*(apc_num+1),6*apc_num+13)]))
						results_m_weight1["STAAR-A(1,1)"] <- CCT(as.numeric(results_m_weight1[1:length(results_m_weight1)][c(1:apc_num+5*(apc_num+1),6*apc_num+14)]))

						results_weight1 <- cbind(results_weight1, results_m_weight1)
						colnames(results_weight1)[i] <- c(i-1)
					}

					#results_weight2
					for(i in 1:ncol(pvalues$results_weight2)){
						results_m_weight2 <- unlist(pvalues$results_weight2[,i][c(5:length(pvalues$results_weight2[,i]), 4,3)])
						names(results_m_weight2) <- names(results_m)[-c(1:5,(length(results_m)-5):length(results_m))]

						results_m_weight2 <- c(results_m_weight2, rep(0,6))
						names(results_m_weight2)[(length(results_m_weight2)-5):length(results_m_weight2)] <- c("SKAT(1,25)-Disruptive","SKAT(1,1)-Disruptive","Burden(1,25)-Disruptive","Burden(1,1)-Disruptive","ACAT-V(1,25)-Disruptive","ACAT-V(1,1)-Disruptive")
						results_m_weight2[(length(results_m_weight2)-5):length(results_m_weight2)] <- unlist(pvalues_dm$results_weight2[,i][c("results_STAAR_S_1_25.SKAT(1,25)","results_STAAR_S_1_1.SKAT(1,1)","results_STAAR_B_1_25.Burden(1,25)",
						                                                                                                                      "results_STAAR_B_1_1.Burden(1,1)","results_STAAR_A_1_25.ACAT-V(1,25)","results_STAAR_A_1_1.ACAT-V(1,1)")])
						results_m_weight2["STAAR-O"] <- CCT(as.numeric(results_m_weight2[1:length(results_m_weight2)][p_seq]))
						results_m_weight2["STAAR-S(1,25)"] <- CCT(as.numeric(results_m_weight2[1:length(results_m_weight2)][c(1:apc_num,6*apc_num+9)]))
						results_m_weight2["STAAR-S(1,1)"] <- CCT(as.numeric(results_m_weight2[1:length(results_m_weight2)][c(1:apc_num+(apc_num+1),6*apc_num+10)]))
						results_m_weight2["STAAR-B(1,25)"] <- CCT(as.numeric(results_m_weight2[1:length(results_m_weight2)][c(1:apc_num+2*(apc_num+1),6*apc_num+11)]))
						results_m_weight2["STAAR-B(1,1)"] <- CCT(as.numeric(results_m_weight2[1:length(results_m_weight2)][c(1:apc_num+3*(apc_num+1),6*apc_num+12)]))
						results_m_weight2["STAAR-A(1,25)"] <- CCT(as.numeric(results_m_weight2[1:length(results_m_weight2)][c(1:apc_num+4*(apc_num+1),6*apc_num+13)]))
						results_m_weight2["STAAR-A(1,1)"] <- CCT(as.numeric(results_m_weight2[1:length(results_m_weight2)][c(1:apc_num+5*(apc_num+1),6*apc_num+14)]))

						results_weight2 <- cbind(results_weight2, results_m_weight2)
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
			}else
			{
				results_m <- c(results[1,],rep(0,2))
				names(results_m)[(length(results_m)-1):length(results_m)] <- c("Burden(1,25)-Disruptive","Burden(1,1)-Disruptive")
				results_m[(length(results_m)-1):length(results_m)] <- results[2,c("Burden(1,25)","Burden(1,1)")]

				## check whether the p-values is NA. If so, set NA equals 1.
				if(is.na(results_m[(length(results_m)-1)]))
				{
					results_m[(length(results_m)-1)] <- 1
				}

				if(is.na(results_m[length(results_m)]))
				{
					results_m[length(results_m)] <- 1
				}

				apc_num <- (length(results_m)-10)/2
				p_seq <- c(1:apc_num,1:apc_num+(apc_num+1),(length(results_m)-6):(length(results_m)-5))

				## calculate STAAR-B
				pvalues_dm_sub <- as.numeric(results_m[6:length(results_m)][p_seq])
				if(sum(is.na(pvalues_dm_sub))>0)
				{
					if(sum(is.na(pvalues_dm_sub))==length(pvalues_dm_sub))
					{
						results_m["STAAR-B"] <- 1
					}else
					{
						## not all NAs
						pvalues_dm_sub <- pvalues_dm_sub[!is.na(pvalues_dm_sub)]
						if(sum(pvalues_dm_sub[pvalues_dm_sub<1])>0)
						{
							## not all ones
							results_m["STAAR-B"] <- CCT(pvalues_dm_sub[pvalues_dm_sub<1])

						}else
						{
							results_m["STAAR-B"] <- 1

						}
					}
				}else
				{
					if(sum(pvalues_dm_sub[pvalues_dm_sub<1])>0)
					{
						results_m["STAAR-B"] <- CCT(pvalues_dm_sub[pvalues_dm_sub<1])
					}else
					{
						results_m["STAAR-B"] <- 1
					}
				}

				## calculate STAAR-B(1,25)
				pvalues_dm_sub <- as.numeric(results_m[6:length(results_m)][c(1:apc_num,(length(results_m)-6))])
				if(sum(is.na(pvalues_dm_sub))>0)
				{
					if(sum(is.na(pvalues_dm_sub))==length(pvalues_dm_sub))
					{
						results_m["STAAR-B(1,25)"] <- 1
					}else
					{
						## not all NAs
						pvalues_dm_sub <- pvalues_dm_sub[!is.na(pvalues_dm_sub)]
						if(sum(pvalues_dm_sub[pvalues_dm_sub<1])>0)
						{
							## not all ones
							results_m["STAAR-B(1,25)"] <- CCT(pvalues_dm_sub[pvalues_dm_sub<1])

						}else
						{
							results_m["STAAR-B(1,25)"] <- 1

						}
					}
				}else
				{
					if(sum(pvalues_dm_sub[pvalues_dm_sub<1])>0)
					{
						results_m["STAAR-B(1,25)"] <- CCT(pvalues_dm_sub[pvalues_dm_sub<1])
					}else
					{
						results_m["STAAR-B(1,25)"] <- 1
					}
				}

				## calculate STAAR-B(1,1)
				pvalues_dm_sub <- as.numeric(results_m[6:length(results_m)][c(1:apc_num+(apc_num+1),(length(results_m)-5))])
				if(sum(is.na(pvalues_dm_sub))>0)
				{
					if(sum(is.na(pvalues_dm_sub))==length(pvalues_dm_sub))
					{
						results_m["STAAR-B(1,1)"] <- 1
					}else
					{
						## not all NAs
						pvalues_dm_sub <- pvalues_dm_sub[!is.na(pvalues_dm_sub)]
						if(sum(pvalues_dm_sub[pvalues_dm_sub<1])>0)
						{
							## not all ones
							results_m["STAAR-B(1,1)"] <- CCT(pvalues_dm_sub[pvalues_dm_sub<1])

						}else
						{
							results_m["STAAR-B(1,1)"] <- 1

						}
					}
				}else
				{
					if(sum(pvalues_dm_sub[pvalues_dm_sub<1])>0)
					{
						results_m["STAAR-B(1,1)"] <- CCT(pvalues_dm_sub[pvalues_dm_sub<1])
					}else
					{
						results_m["STAAR-B(1,1)"] <- 1
					}
				}

				results <- c()
				results <- rbind(results,results_m)

				if(use_ancestry_informed == TRUE & find_weight == TRUE){
					results_weight <- results_weight1 <- results_weight2 <- c()

					#results_weight
					for(i in 1:ncol(pvalues$results_weight)){
						results_m_weight <- c(pvalues$results_weight[-c(1,2,3),i], pvalues$results_weight[c(3),i])
						names(results_m_weight) <- colnames(results)[c(6:(6+apc_num), (7+apc_num):(7+2*apc_num+1))]
						results_m_weight <- c(results_m_weight,rep(0,2))
						names(results_m_weight)[(length(results_m_weight)-1):length(results_m_weight)] <- c("Burden(1,25)-Disruptive","Burden(1,1)-Disruptive")
						results_m_weight[(length(results_m_weight)-1):length(results_m_weight)] <- unlist(pvalues_dm$results_weight[,i][c("results_STAAR_B_1_25.Burden(1,25)","results_STAAR_B_1_1.Burden(1,1)")])

						## check whether the p-values are NA. If so, set NA equals 1.
						if(is.na(results_m_weight[(length(results_m_weight)-1)]))
						{
							results_m_weight[(length(results_m_weight)-1)] <- 1
						}

						if(is.na(results_m_weight[length(results_m_weight)]))
						{
							results_m_weight[length(results_m_weight)] <- 1
						}

						## calculate STAAR-B
						pvalues_dm_sub <- as.numeric(unlist(results_m_weight[1:length(results_m_weight)][p_seq]))
						if(sum(is.na(pvalues_dm_sub))>0)
						{
							if(sum(is.na(pvalues_dm_sub))==length(pvalues_dm_sub))
							{
								results_m_weight["STAAR-B"] <- 1
							}else
							{
								## not all NAs
								pvalues_dm_sub <- pvalues_dm_sub[!is.na(pvalues_dm_sub)]
								if(sum(pvalues_dm_sub[pvalues_dm_sub<1])>0)
								{
									## not all ones
									results_m_weight["STAAR-B"] <- CCT(pvalues_dm_sub[pvalues_dm_sub<1])

								}else
								{
									results_m_weight["STAAR-B"] <- 1

								}
							}
						}else
						{
							if(sum(pvalues_dm_sub[pvalues_dm_sub<1])>0)
							{
								results_m_weight["STAAR-B"] <- CCT(pvalues_dm_sub[pvalues_dm_sub<1])
							}else
							{
								results_m_weight["STAAR-B"] <- 1
							}
						}

						## calculate STAAR-B(1,25)
						pvalues_dm_sub <- as.numeric(unlist(results_m_weight[1:length(results_m_weight)][c(1:apc_num,(length(results_m_weight)-1))]))
						if(sum(is.na(pvalues_dm_sub))>0)
						{
							if(sum(is.na(pvalues_dm_sub))==length(pvalues_dm_sub))
							{
								results_m_weight["STAAR-B(1,25)"] <- 1
							}else
							{
								## not all NAs
								pvalues_dm_sub <- pvalues_dm_sub[!is.na(pvalues_dm_sub)]
								if(sum(pvalues_dm_sub[pvalues_dm_sub<1])>0)
								{
									## not all ones
									results_m_weight["STAAR-B(1,25)"] <- CCT(pvalues_dm_sub[pvalues_dm_sub<1])

								}else
								{
									results_m_weight["STAAR-B(1,25)"] <- 1

								}
							}
						}else
						{
							if(sum(pvalues_dm_sub[pvalues_dm_sub<1])>0)
							{
								results_m_weight["STAAR-B(1,25)"] <- CCT(pvalues_dm_sub[pvalues_dm_sub<1])
							}else
							{
								results_m_weight["STAAR-B(1,25)"] <- 1
							}
						}

						## calculate STAAR-B(1,1)
						pvalues_dm_sub <- as.numeric(unlist(results_m_weight[1:length(results_m_weight)][c(1:apc_num+(apc_num+1),(length(results_m_weight)))]))
						if(sum(is.na(pvalues_dm_sub))>0)
						{
							if(sum(is.na(pvalues_dm_sub))==length(pvalues_dm_sub))
							{
								results_m_weight["STAAR-B(1,1)"] <- 1
							}else
							{
								## not all NAs
								pvalues_dm_sub <- pvalues_dm_sub[!is.na(pvalues_dm_sub)]
								if(sum(pvalues_dm_sub[pvalues_dm_sub<1])>0)
								{
									## not all ones
									results_m_weight["STAAR-B(1,1)"] <- CCT(pvalues_dm_sub[pvalues_dm_sub<1])

								}else
								{
									results_m_weight["STAAR-B(1,1)"] <- 1

								}
							}
						}else
						{
							if(sum(pvalues_dm_sub[pvalues_dm_sub<1])>0)
							{
								results_m_weight["STAAR-B(1,1)"] <- CCT(pvalues_dm_sub[pvalues_dm_sub<1])
							}else
							{
								results_m_weight["STAAR-B(1,1)"] <- 1
							}
						}
						results_weight <- cbind(results_weight, results_m_weight)
						colnames(results_weight)[i] <- c(i-1)

					}

					#results_weight1
					for(i in 1:ncol(pvalues$results_weight1)){
						results_m_weight1 <- c(pvalues$results_weight1[-c(1,2,3),i], pvalues$results_weight1[c(3),i])
						names(results_m_weight1) <- colnames(results)[c(6:(6+apc_num), (7+apc_num):(7+2*apc_num+1))]
						results_m_weight1 <- c(results_m_weight1,rep(0,2))
						names(results_m_weight1)[(length(results_m_weight1)-1):length(results_m_weight1)] <- c("Burden(1,25)-Disruptive","Burden(1,1)-Disruptive")
						results_m_weight1[(length(results_m_weight1)-1):length(results_m_weight1)] <- unlist(pvalues_dm$results_weight1[,i][c("results_STAAR_B_1_25.Burden(1,25)","results_STAAR_B_1_1.Burden(1,1)")])

						## check whether the p-values are NA. If so, set NA equals 1.
						if(is.na(results_m_weight1[(length(results_m_weight1)-1)]))
						{
							results_m_weight1[(length(results_m_weight1)-1)] <- 1
						}

						if(is.na(results_m_weight1[length(results_m_weight1)]))
						{
							results_m_weight1[length(results_m_weight1)] <- 1
						}

						## calculate STAAR-B
						pvalues_dm_sub <- as.numeric(unlist(results_m_weight1[1:length(results_m_weight1)][p_seq]))
						if(sum(is.na(pvalues_dm_sub))>0)
						{
							if(sum(is.na(pvalues_dm_sub))==length(pvalues_dm_sub))
							{
								results_m_weight1["STAAR-B"] <- 1
							}else
							{
								## not all NAs
								pvalues_dm_sub <- pvalues_dm_sub[!is.na(pvalues_dm_sub)]
								if(sum(pvalues_dm_sub[pvalues_dm_sub<1])>0)
								{
									## not all ones
									results_m_weight1["STAAR-B"] <- CCT(pvalues_dm_sub[pvalues_dm_sub<1])

								}else
								{
									results_m_weight1["STAAR-B"] <- 1

								}
							}
						}else
						{
							if(sum(pvalues_dm_sub[pvalues_dm_sub<1])>0)
							{
								results_m_weight1["STAAR-B"] <- CCT(pvalues_dm_sub[pvalues_dm_sub<1])
							}else
							{
								results_m_weight1["STAAR-B"] <- 1
							}
						}

						## calculate STAAR-B(1,25)
						pvalues_dm_sub <- as.numeric(unlist(results_m_weight1[1:length(results_m_weight1)][c(1:apc_num,(length(results_m_weight1)-1))]))
						if(sum(is.na(pvalues_dm_sub))>0)
						{
							if(sum(is.na(pvalues_dm_sub))==length(pvalues_dm_sub))
							{
								results_m_weight1["STAAR-B(1,25)"] <- 1
							}else
							{
								## not all NAs
								pvalues_dm_sub <- pvalues_dm_sub[!is.na(pvalues_dm_sub)]
								if(sum(pvalues_dm_sub[pvalues_dm_sub<1])>0)
								{
									## not all ones
									results_m_weight1["STAAR-B(1,25)"] <- CCT(pvalues_dm_sub[pvalues_dm_sub<1])

								}else
								{
									results_m_weight1["STAAR-B(1,25)"] <- 1

								}
							}
						}else
						{
							if(sum(pvalues_dm_sub[pvalues_dm_sub<1])>0)
							{
								results_m_weight1["STAAR-B(1,25)"] <- CCT(pvalues_dm_sub[pvalues_dm_sub<1])
							}else
							{
								results_m_weight1["STAAR-B(1,25)"] <- 1
							}
						}

						## calculate STAAR-B(1,1)
						pvalues_dm_sub <- as.numeric(unlist(results_m_weight1[1:length(results_m_weight1)][c(1:apc_num+(apc_num+1),(length(results_m_weight1)))]))
						if(sum(is.na(pvalues_dm_sub))>0)
						{
							if(sum(is.na(pvalues_dm_sub))==length(pvalues_dm_sub))
							{
								results_m_weight1["STAAR-B(1,1)"] <- 1
							}else
							{
								## not all NAs
								pvalues_dm_sub <- pvalues_dm_sub[!is.na(pvalues_dm_sub)]
								if(sum(pvalues_dm_sub[pvalues_dm_sub<1])>0)
								{
									## not all ones
									results_m_weight1["STAAR-B(1,1)"] <- CCT(pvalues_dm_sub[pvalues_dm_sub<1])

								}else
								{
									results_m_weight1["STAAR-B(1,1)"] <- 1

								}
							}
						}else
						{
							if(sum(pvalues_dm_sub[pvalues_dm_sub<1])>0)
							{
								results_m_weight1["STAAR-B(1,1)"] <- CCT(pvalues_dm_sub[pvalues_dm_sub<1])
							}else
							{
								results_m_weight1["STAAR-B(1,1)"] <- 1
							}
						}
						results_weight1 <- cbind(results_weight1, results_m_weight1)
						colnames(results_weight1)[i] <- c(i-1)

					}

					#results_weight2
					for(i in 1:ncol(pvalues$results_weight2)){
						results_m_weight2 <- c(pvalues$results_weight2[-c(1,2,3),i], pvalues$results_weight2[c(3),i])
						names(results_m_weight2) <- colnames(results)[c(6:(6+apc_num), (7+apc_num):(7+2*apc_num+1))]
						results_m_weight2 <- c(results_m_weight2,rep(0,2))
						names(results_m_weight2)[(length(results_m_weight2)-1):length(results_m_weight2)] <- c("Burden(1,25)-Disruptive","Burden(1,1)-Disruptive")
						results_m_weight2[(length(results_m_weight2)-1):length(results_m_weight2)] <- unlist(pvalues_dm$results_weight2[,i][c("results_STAAR_B_1_25.Burden(1,25)","results_STAAR_B_1_1.Burden(1,1)")])

						## check whether the p-values is NA. If so, set NA equals 1.
						if(is.na(results_m_weight2[(length(results_m_weight2)-1)]))
						{
							results_m_weight2[(length(results_m_weight2)-1)] <- 1
						}

						if(is.na(results_m_weight2[length(results_m_weight2)]))
						{
							results_m_weight2[length(results_m_weight2)] <- 1
						}

						## calculate STAAR-B
						pvalues_dm_sub <- as.numeric(unlist(results_m_weight2[1:length(results_m_weight2)][p_seq]))
						if(sum(is.na(pvalues_dm_sub))>0)
						{
							if(sum(is.na(pvalues_dm_sub))==length(pvalues_dm_sub))
							{
								results_m_weight2["STAAR-B"] <- 1
							}else
							{
								## not all NAs
								pvalues_dm_sub <- pvalues_dm_sub[!is.na(pvalues_dm_sub)]
								if(sum(pvalues_dm_sub[pvalues_dm_sub<1])>0)
								{
									## not all ones
									results_m_weight2["STAAR-B"] <- CCT(pvalues_dm_sub[pvalues_dm_sub<1])

								}else
								{
									results_m_weight2["STAAR-B"] <- 1

								}
							}
						}else
						{
							if(sum(pvalues_dm_sub[pvalues_dm_sub<1])>0)
							{
								results_m_weight2["STAAR-B"] <- CCT(pvalues_dm_sub[pvalues_dm_sub<1])
							}else
							{
								results_m_weight2["STAAR-B"] <- 1
							}
						}

						## calculate STAAR-B(1,25)
						pvalues_dm_sub <- as.numeric(unlist(results_m_weight2[1:length(results_m_weight2)][c(1:apc_num,(length(results_m_weight2)-1))]))
						if(sum(is.na(pvalues_dm_sub))>0)
						{
							if(sum(is.na(pvalues_dm_sub))==length(pvalues_dm_sub))
							{
								results_m_weight2["STAAR-B(1,25)"] <- 1
							}else
							{
								## not all NAs
								pvalues_dm_sub <- pvalues_dm_sub[!is.na(pvalues_dm_sub)]
								if(sum(pvalues_dm_sub[pvalues_dm_sub<1])>0)
								{
									## not all ones
									results_m_weight2["STAAR-B(1,25)"] <- CCT(pvalues_dm_sub[pvalues_dm_sub<1])

								}else
								{
									results_m_weight2["STAAR-B(1,25)"] <- 1

								}
							}
						}else
						{
							if(sum(pvalues_dm_sub[pvalues_dm_sub<1])>0)
							{
								results_m_weight2["STAAR-B(1,25)"] <- CCT(pvalues_dm_sub[pvalues_dm_sub<1])
							}else
							{
								results_m_weight2["STAAR-B(1,25)"] <- 1
							}
						}

						## calculate STAAR-B(1,1)
						pvalues_dm_sub <- as.numeric(unlist(results_m_weight2[1:length(results_m_weight2)][c(1:apc_num+(apc_num+1),(length(results_m_weight2)))]))
						if(sum(is.na(pvalues_dm_sub))>0)
						{
							if(sum(is.na(pvalues_dm_sub))==length(pvalues_dm_sub))
							{
								results_m_weight2["STAAR-B(1,1)"] <- 1
							}else
							{
								## not all NAs
								pvalues_dm_sub <- pvalues_dm_sub[!is.na(pvalues_dm_sub)]
								if(sum(pvalues_dm_sub[pvalues_dm_sub<1])>0)
								{
									## not all ones
									results_m_weight2["STAAR-B(1,1)"] <- CCT(pvalues_dm_sub[pvalues_dm_sub<1])

								}else
								{
									results_m_weight2["STAAR-B(1,1)"] <- 1

								}
							}
						}else
						{
							if(sum(pvalues_dm_sub[pvalues_dm_sub<1])>0)
							{
								results_m_weight2["STAAR-B(1,1)"] <- CCT(pvalues_dm_sub[pvalues_dm_sub<1])
							}else
							{
								results_m_weight2["STAAR-B(1,1)"] <- 1
							}
						}
						results_weight2 <- cbind(results_weight2, results_m_weight2)
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
		}
	}else
	{
		results <- c()
	}

	seqResetFilter(genofile)

	return(results)
}

