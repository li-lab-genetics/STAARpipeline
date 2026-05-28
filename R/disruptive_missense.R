disruptive_missense <- function(chr,gene_name,genofile,obj_nullmodel,genes,
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

	rm(filter)
	gc()

	### Gene
	kk <- which(genes[,1]==gene_name)

	sub_start_loc <- genes[kk,3]
	sub_end_loc <- genes[kk,4]

	is.in <- (SNVlist)&(position>=sub_start_loc)&(position<=sub_end_loc)
	variant.id.gene <- variant.id[is.in]

	seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

	## disruptive_missense
	## Gencode_Exonic
	GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
	## Meta.SVM.Pred
	MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))

	variant.id.gene <- seqGetData(genofile, "variant.id")
	lof.in.ds <- ((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
	variant.id.gene <- variant.id.gene[lof.in.ds]

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
		results_temp[3] <- "disruptive_missense"
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
				results_weight_temp <- unlist(pvalues$results_weight[,i][c(5:length(pvalues$results_weight[,i]), 4,3)])
				names(results_weight_temp) <- colnames(results)[-c(1:5)]

				results_weight <- cbind(results_weight, results_weight_temp)
				colnames(results_weight)[i] <- c(i-1)
			}

			for(i in 1:ncol(pvalues$results_weight1)){
				results_weight_temp1 <- unlist(pvalues$results_weight1[,i][c(5:length(pvalues$results_weight1[,i]), 4,3)])
				names(results_weight_temp1) <- colnames(results)[-c(1:5)]

				results_weight1 <- cbind(results_weight1, results_weight_temp1)
				colnames(results_weight1)[i] <- c(i-1)
			}

			for(i in 1:ncol(pvalues$results_weight2)){
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

