coding_incl_ptv_PheWAS <- function(chr,gene_name,genofile,obj_nullmodel_list,genes,
                                   rare_maf_cutoff=0.01,rv_num_cutoff=2,rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                                   QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                   Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                   Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                   SPA_p_filter=FALSE,p_filter_cutoff=0.05,silent=FALSE){

	## evaluate choices
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	phenotype.id <- Reduce(union, lapply(obj_nullmodel_list, function(x){
		as.character(x$id_include)
	}))

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

	rm(position)
	gc()

	seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

	## Gencode_Exonic
	GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
	## Gencode
	GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
	## Meta.SVM.Pred
	MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))

	################################################
	#           Coding
	################################################
	variant.id.gene <- seqGetData(genofile, "variant.id")
	lof.in.coding <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|(GENCODE.EXONIC.Category=="nonsynonymous SNV")|(GENCODE.EXONIC.Category=="synonymous SNV")|(GENCODE.EXONIC.Category=="frameshift deletion")|(GENCODE.EXONIC.Category=="frameshift insertion")
	variant.id.gene <- variant.id.gene[lof.in.coding]
	
	if(length(variant.id.gene)==0)
	{
	  results_list_sub <- rep(list(c()), length(obj_nullmodel_list))
	  results_coding_list <- list(plof = results_list_sub,
	                              plof_ds = results_list_sub,
	                              missense = results_list_sub,
	                              disruptive_missense = results_list_sub,
	                              synonymous = results_list_sub,
	                              ptv = results_list_sub,
	                              ptv_ds = results_list_sub)
	  seqResetFilter(genofile)
	  return(results_coding_list)
	}

	seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

	## Gencode_Exonic
	GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))
	## Gencode
	GENCODE.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.Category")]))
	## Meta.SVM.Pred
	MetaSVM_pred <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="MetaSVM")]))

	## Annotation
	Anno.Int.PHRED.sub <- NULL
	Anno.Int.PHRED.sub.name <- NULL

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
	
	## get AF, Missing rate
	AF_AC_Missing <- seqGetAF_AC_Missing(genofile,minor=FALSE,parallel=FALSE)
	REF_AF <- AF_AC_Missing$af
	Missing_rate <- AF_AC_Missing$miss
	variant_maf_cutoff_filter <- ifelse(rare_maf_cutoff<=0.01,0.05,1)

	################################################
	#                  plof_ds
	################################################
	variant.id.gene <- seqGetData(genofile, "variant.id")
	lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")|((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
	variant.id.gene.category <- variant.id.gene[lof.in.plof]

	## Annotation
	Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.plof,]
	
	Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.id.gene.category,
	                                      sample.id=phenotype.id,
	                                      REF_AF=REF_AF[lof.in.plof],variant_maf_cutoff_filter=variant_maf_cutoff_filter,
	                                      Missing_rate=Missing_rate[lof.in.plof],
	                                      rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
	                                      annotation_phred=Anno.Int.PHRED.sub.category)
	Geno <- Genotype_sp$Geno
	Anno.Int.PHRED.sub.category <- Genotype_sp$annotation_phred
	rm(Genotype_sp)
	gc()

	results_list_plof_ds <- rep(list(c()), length(obj_nullmodel_list))
	phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))

	if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
	{
	  for(num in 1:length(obj_nullmodel_list))
	  {
	    id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
	    id.phenotype.num.merge <- data.frame(id.phenotype.num)
	    id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
	    phenotype.id.match <- id.phenotype.num.merge$index
	    samplesize.num <- length(phenotype.id.match)
	    
	    ## Extract the genotype matrix and variant information of the current trait
	    Geno.num <- Geno[phenotype.id.match, , drop = FALSE]
	    MAC.in <- Matrix::colSums(Geno.num,na.rm = TRUE)
	    Missing_num.in <- Missing_num.sp(Geno.num)
	    MAF.in <- MAC.in/(2*(samplesize.num-Missing_num.in))
	    
	    if (geno_missing_imputation == "mean")
	    {
	      Geno.num <- na.replace.sp(Geno.num,m=2*MAF.in)
	    }
	    if (geno_missing_imputation == "minor")
	    {
	      Geno.num <- na.replace.sp(Geno.num,is_NA_to_Zero=TRUE)
	      MAF.in <- MAC.in/(2*samplesize.num)
	    }
	    
	    ## number of traits in analysis
	    n_pheno <- obj_nullmodel_list[[num]]$n.pheno
	    
	    ## SPA status
	    if(!is.null(obj_nullmodel_list[[num]]$use_SPA))
	    {
	      use_SPA <- obj_nullmodel_list[[num]]$use_SPA
	    }else
	    {
	      use_SPA <- FALSE
	    }
	    
	    pvalues <- 0
	    if(n_pheno == 1)
	    {
	      if(!use_SPA)
	      {
	        try(pvalues <- STAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
	      }else{
	        try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
	      }
	    }else
	    {
	      try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
	    }
	    
	    results_plof_ds <- c()
	    if(inherits(pvalues, "list"))
	    {
	      results_temp <- as.vector(genes[kk,])
	      results_temp[3] <- "plof_ds"
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
	      
	      results_plof_ds <- rbind(results_plof_ds,results_temp)
	    }
	    
	    if(!is.null(results_plof_ds))
	    {
	      if(!use_SPA)
	      {
	        colnames(results_plof_ds) <- colnames(results_plof_ds, do.NULL = FALSE, prefix = "col")
	        colnames(results_plof_ds)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
	        colnames(results_plof_ds)[(dim(results_plof_ds)[2]-1):dim(results_plof_ds)[2]] <- c("ACAT-O","STAAR-O")
	      }else
	      {
	        colnames(results_plof_ds) <- colnames(results_plof_ds, do.NULL = FALSE, prefix = "col")
	        colnames(results_plof_ds)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
	        colnames(results_plof_ds)[dim(results_plof_ds)[2]] <- c("STAAR-B")
	      }
	      
	      results_list_plof_ds[[num]] <- results_plof_ds
	    }
	  }
	}

	################################################
	#                  ptv_ds
	################################################
	if(variant_type=="SNV")
	{
		lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
	}

	if(variant_type=="Indel")
	{
		lof.in.plof <- (GENCODE.EXONIC.Category=="frameshift deletion")|(GENCODE.EXONIC.Category=="frameshift insertion")|((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
	}

	if(variant_type=="variant")
	{
		lof.in.plof.snv <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")
		lof.in.plof.indel <- (GENCODE.EXONIC.Category=="frameshift deletion")|(GENCODE.EXONIC.Category=="frameshift insertion")
		lof.in.plof <- lof.in.plof.snv|lof.in.plof.indel|((GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D"))
	}

	variant.id.gene.category <- variant.id.gene[lof.in.plof]
	
	## Annotation
	Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.plof,]
	
	Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.id.gene.category,
	                                      sample.id=phenotype.id,
	                                      REF_AF=REF_AF[lof.in.plof],variant_maf_cutoff_filter=variant_maf_cutoff_filter,
	                                      Missing_rate=Missing_rate[lof.in.plof],
	                                      rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
	                                      annotation_phred=Anno.Int.PHRED.sub.category)
	Geno <- Genotype_sp$Geno
	Anno.Int.PHRED.sub.category <- Genotype_sp$annotation_phred
	rm(Genotype_sp)
	gc()

	results_list_ptv_ds <- rep(list(c()), length(obj_nullmodel_list))
	phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))

	if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
	{
	  for(num in 1:length(obj_nullmodel_list))
	  {
	    id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
	    id.phenotype.num.merge <- data.frame(id.phenotype.num)
	    id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
	    phenotype.id.match <- id.phenotype.num.merge$index
	    samplesize.num <- length(phenotype.id.match)
	    
	    ## Extract the genotype matrix and variant information of the current trait
	    Geno.num <- Geno[phenotype.id.match, , drop = FALSE]
	    MAC.in <- Matrix::colSums(Geno.num,na.rm = TRUE)
	    Missing_num.in <- Missing_num.sp(Geno.num)
	    MAF.in <- MAC.in/(2*(samplesize.num-Missing_num.in))
	    
	    if (geno_missing_imputation == "mean")
	    {
	      Geno.num <- na.replace.sp(Geno.num,m=2*MAF.in)
	    }
	    if (geno_missing_imputation == "minor")
	    {
	      Geno.num <- na.replace.sp(Geno.num,is_NA_to_Zero=TRUE)
	      MAF.in <- MAC.in/(2*samplesize.num)
	    }
	    
	    ## number of traits in analysis
	    n_pheno <- obj_nullmodel_list[[num]]$n.pheno
	    
	    ## SPA status
	    if(!is.null(obj_nullmodel_list[[num]]$use_SPA))
	    {
	      use_SPA <- obj_nullmodel_list[[num]]$use_SPA
	    }else
	    {
	      use_SPA <- FALSE
	    }
	    
	    pvalues <- 0
	    if(n_pheno == 1)
	    {
	      if(!use_SPA)
	      {
	        try(pvalues <- STAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
	      }else{
	        try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
	      }
	    }else
	    {
	      try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
	    }
	    
	    results_ptv_ds <- c()
	    if(inherits(pvalues, "list"))
	    {
	      results_temp <- as.vector(genes[kk,])
	      results_temp[3] <- "ptv_ds"
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
	      
	      results_ptv_ds <- rbind(results_ptv_ds,results_temp)
	    }
	    
	    if(!is.null(results_ptv_ds))
	    {
	      if(!use_SPA)
	      {
	        colnames(results_ptv_ds) <- colnames(results_ptv_ds, do.NULL = FALSE, prefix = "col")
	        colnames(results_ptv_ds)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
	        colnames(results_ptv_ds)[(dim(results_ptv_ds)[2]-1):dim(results_ptv_ds)[2]] <- c("ACAT-O","STAAR-O")
	      }else
	      {
	        colnames(results_ptv_ds) <- colnames(results_ptv_ds, do.NULL = FALSE, prefix = "col")
	        colnames(results_ptv_ds)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
	        colnames(results_ptv_ds)[dim(results_ptv_ds)[2]] <- c("STAAR-B")
	      }
	      
	      results_list_ptv_ds[[num]] <- results_ptv_ds
	    }
	  }
	}

	#####################################################
	#                      plof
	#####################################################
	lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")|(GENCODE.Category=="ncRNA_splicing")|(GENCODE.Category=="ncRNA_exonic;splicing")
	variant.id.gene.category <- variant.id.gene[lof.in.plof]

	## Annotation
	Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.plof,]
	
	Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.id.gene.category,
	                                      sample.id=phenotype.id,
	                                      REF_AF=REF_AF[lof.in.plof],variant_maf_cutoff_filter=variant_maf_cutoff_filter,
	                                      Missing_rate=Missing_rate[lof.in.plof],
	                                      rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
	                                      annotation_phred=Anno.Int.PHRED.sub.category)
	Geno <- Genotype_sp$Geno
	Anno.Int.PHRED.sub.category <- Genotype_sp$annotation_phred
	rm(Genotype_sp)
	gc()

	results_list_plof <- rep(list(c()), length(obj_nullmodel_list))
	phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))

	if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
	{
	  for(num in 1:length(obj_nullmodel_list))
	  {
	    id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
	    id.phenotype.num.merge <- data.frame(id.phenotype.num)
	    id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
	    phenotype.id.match <- id.phenotype.num.merge$index
	    samplesize.num <- length(phenotype.id.match)
	    
	    ## Extract the genotype matrix and variant information of the current trait
	    Geno.num <- Geno[phenotype.id.match, , drop = FALSE]
	    MAC.in <- Matrix::colSums(Geno.num,na.rm = TRUE)
	    Missing_num.in <- Missing_num.sp(Geno.num)
	    MAF.in <- MAC.in/(2*(samplesize.num-Missing_num.in))
	    
	    if (geno_missing_imputation == "mean")
	    {
	      Geno.num <- na.replace.sp(Geno.num,m=2*MAF.in)
	    }
	    if (geno_missing_imputation == "minor")
	    {
	      Geno.num <- na.replace.sp(Geno.num,is_NA_to_Zero=TRUE)
	      MAF.in <- MAC.in/(2*samplesize.num)
	    }
	    
	    ## number of traits in analysis
	    n_pheno <- obj_nullmodel_list[[num]]$n.pheno
	    
	    ## SPA status
	    if(!is.null(obj_nullmodel_list[[num]]$use_SPA))
	    {
	      use_SPA <- obj_nullmodel_list[[num]]$use_SPA
	    }else
	    {
	      use_SPA <- FALSE
	    }
	    
	    pvalues <- 0
	    if(n_pheno == 1)
	    {
	      if(!use_SPA)
	      {
	        try(pvalues <- STAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
	      }else{
	        try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
	      }
	    }else
	    {
	      try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
	    }
	    
	    results_plof <- c()
	    if(inherits(pvalues, "list"))
	    {
	      results_temp <- as.vector(genes[kk,])
	      results_temp[3] <- "plof"
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
	      
	      results_plof <- rbind(results_plof,results_temp)
	    }
	    
	    if(!is.null(results_plof))
	    {
	      if(!use_SPA)
	      {
	        colnames(results_plof) <- colnames(results_plof, do.NULL = FALSE, prefix = "col")
	        colnames(results_plof)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
	        colnames(results_plof)[(dim(results_plof)[2]-1):dim(results_plof)[2]] <- c("ACAT-O","STAAR-O")
	      }else
	      {
	        colnames(results_plof) <- colnames(results_plof_ds, do.NULL = FALSE, prefix = "col")
	        colnames(results_plof)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
	        colnames(results_plof)[dim(results_plof)[2]] <- c("STAAR-B")
	      }
	      
	      results_list_plof[[num]] <- results_plof
	    }
	  }
	}

	################################################
	#                  ptv
	################################################
	if(variant_type=="SNV")
	{
		lof.in.plof <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")
	}

	if(variant_type=="Indel")
	{
		lof.in.plof <- (GENCODE.EXONIC.Category=="frameshift deletion")|(GENCODE.EXONIC.Category=="frameshift insertion")
	}

	if(variant_type=="variant")
	{
		lof.in.plof.snv <- (GENCODE.EXONIC.Category=="stopgain")|(GENCODE.EXONIC.Category=="stoploss")|(GENCODE.Category=="splicing")|(GENCODE.Category=="exonic;splicing")
		lof.in.plof.indel <- (GENCODE.EXONIC.Category=="frameshift deletion")|(GENCODE.EXONIC.Category=="frameshift insertion")
		lof.in.plof <- lof.in.plof.snv|lof.in.plof.indel
	}

	variant.id.gene.category <- variant.id.gene[lof.in.plof]

	## Annotation
	Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.plof,]
	
	Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.id.gene.category,
	                                      sample.id=phenotype.id,
	                                      REF_AF=REF_AF[lof.in.plof],variant_maf_cutoff_filter=variant_maf_cutoff_filter,
	                                      Missing_rate=Missing_rate[lof.in.plof],
	                                      rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
	                                      annotation_phred=Anno.Int.PHRED.sub.category)
	Geno <- Genotype_sp$Geno
	Anno.Int.PHRED.sub.category <- Genotype_sp$annotation_phred
	rm(Genotype_sp)
	gc()

	results_list_ptv <- rep(list(c()), length(obj_nullmodel_list))
	phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))

	if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
	{
	  for(num in 1:length(obj_nullmodel_list))
	  {
	    id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
	    id.phenotype.num.merge <- data.frame(id.phenotype.num)
	    id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
	    phenotype.id.match <- id.phenotype.num.merge$index
	    samplesize.num <- length(phenotype.id.match)
	    
	    ## Extract the genotype matrix and variant information of the current trait
	    Geno.num <- Geno[phenotype.id.match, , drop = FALSE]
	    MAC.in <- Matrix::colSums(Geno.num,na.rm = TRUE)
	    Missing_num.in <- Missing_num.sp(Geno.num)
	    MAF.in <- MAC.in/(2*(samplesize.num-Missing_num.in))
	    
	    if (geno_missing_imputation == "mean")
	    {
	      Geno.num <- na.replace.sp(Geno.num,m=2*MAF.in)
	    }
	    if (geno_missing_imputation == "minor")
	    {
	      Geno.num <- na.replace.sp(Geno.num,is_NA_to_Zero=TRUE)
	      MAF.in <- MAC.in/(2*samplesize.num)
	    }
	    
	    
	    ## number of traits in analysis
	    n_pheno <- obj_nullmodel_list[[num]]$n.pheno
	    
	    ## SPA status
	    if(!is.null(obj_nullmodel_list[[num]]$use_SPA))
	    {
	      use_SPA <- obj_nullmodel_list[[num]]$use_SPA
	    }else
	    {
	      use_SPA <- FALSE
	    }
	    
	    pvalues <- 0
	    if(n_pheno == 1)
	    {
	      if(!use_SPA)
	      {
	        try(pvalues <- STAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
	      }else{
	        try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
	      }
	    }else
	    {
	      try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
	    }
	    
	    results_ptv <- c()
	    if(inherits(pvalues, "list"))
	    {
	      results_temp <- as.vector(genes[kk,])
	      results_temp[3] <- "ptv"
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
	      
	      results_ptv <- rbind(results_ptv,results_temp)
	    }
	    
	    if(!is.null(results_ptv))
	    {
	      if(!use_SPA)
	      {
	        colnames(results_ptv) <- colnames(results_ptv, do.NULL = FALSE, prefix = "col")
	        colnames(results_ptv)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
	        colnames(results_ptv)[(dim(results_ptv)[2]-1):dim(results_ptv)[2]] <- c("ACAT-O","STAAR-O")
	      }else
	      {
	        colnames(results_ptv) <- colnames(results_ptv, do.NULL = FALSE, prefix = "col")
	        colnames(results_ptv)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
	        colnames(results_ptv)[dim(results_ptv)[2]] <- c("STAAR-B")
	      }
	      
	      results_list_ptv[[num]] <- results_ptv
	    }
	  }
	}

	#############################################
	#             synonymous
	#############################################
	lof.in.synonymous <- (GENCODE.EXONIC.Category=="synonymous SNV")
	variant.id.gene.category <- variant.id.gene[lof.in.synonymous]

	## Annotation
	Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.synonymous,]
	
	Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.id.gene.category,
	                                      sample.id=phenotype.id,
	                                      REF_AF=REF_AF[lof.in.synonymous],variant_maf_cutoff_filter=variant_maf_cutoff_filter,
	                                      Missing_rate=Missing_rate[lof.in.synonymous],
	                                      rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
	                                      annotation_phred=Anno.Int.PHRED.sub.category)
	Geno <- Genotype_sp$Geno
	Anno.Int.PHRED.sub.category <- Genotype_sp$annotation_phred
	rm(Genotype_sp)
	gc()

	results_list_synonymous <- rep(list(c()), length(obj_nullmodel_list))
	phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))

	if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
	{
	  for(num in 1:length(obj_nullmodel_list))
	  {
	    id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
	    id.phenotype.num.merge <- data.frame(id.phenotype.num)
	    id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
	    phenotype.id.match <- id.phenotype.num.merge$index
	    samplesize.num <- length(phenotype.id.match)
	    
	    ## Extract the genotype matrix and variant information of the current trait
	    Geno.num <- Geno[phenotype.id.match, , drop = FALSE]
	    MAC.in <- Matrix::colSums(Geno.num,na.rm = TRUE)
	    Missing_num.in <- Missing_num.sp(Geno.num)
	    MAF.in <- MAC.in/(2*(samplesize.num-Missing_num.in))
	    
	    if (geno_missing_imputation == "mean")
	    {
	      Geno.num <- na.replace.sp(Geno.num,m=2*MAF.in)
	    }
	    if (geno_missing_imputation == "minor")
	    {
	      Geno.num <- na.replace.sp(Geno.num,is_NA_to_Zero=TRUE)
	      MAF.in <- MAC.in/(2*samplesize.num)
	    }
	    
	    ## number of traits in analysis
	    n_pheno <- obj_nullmodel_list[[num]]$n.pheno
	    
	    ## SPA status
	    if(!is.null(obj_nullmodel_list[[num]]$use_SPA))
	    {
	      use_SPA <- obj_nullmodel_list[[num]]$use_SPA
	    }else
	    {
	      use_SPA <- FALSE
	    }
	    
	    pvalues <- 0
	    if(n_pheno == 1)
	    {
	      if(!use_SPA)
	      {
	        try(pvalues <- STAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
	      }else{
	        try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
	      }
	    }else
	    {
	      try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
	    }
	    
	    results_synonymous <- c()
	    if(inherits(pvalues, "list"))
	    {
	      results_temp <- as.vector(genes[kk,])
	      results_temp[3] <- "synonymous"
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
	      
	      results_synonymous <- rbind(results_synonymous,results_temp)
	    }
	    
	    if(!is.null(results_synonymous))
	    {
	      if(!use_SPA)
	      {
	        colnames(results_synonymous) <- colnames(results_synonymous, do.NULL = FALSE, prefix = "col")
	        colnames(results_synonymous)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
	        colnames(results_synonymous)[(dim(results_synonymous)[2]-1):dim(results_synonymous)[2]] <- c("ACAT-O","STAAR-O")
	      }else
	      {
	        colnames(results_synonymous) <- colnames(results_synonymous, do.NULL = FALSE, prefix = "col")
	        colnames(results_synonymous)[1:5] <- c("Gene name","Chr","Category","#SNV","cMAC")
	        colnames(results_synonymous)[dim(results_synonymous)[2]] <- c("STAAR-B")
	      }
	      
	      results_list_synonymous[[num]] <- results_synonymous
	    }
	  }
	}

	#################################################
	#        missense
	#################################################
	lof.in.missense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")
	variant.id.gene.category <- variant.id.gene[lof.in.missense]

	## Annotation
	Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.missense,]
	
	Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.id.gene.category,
	                                      sample.id=phenotype.id,
	                                      REF_AF=REF_AF[lof.in.missense],variant_maf_cutoff_filter=variant_maf_cutoff_filter,
	                                      Missing_rate=Missing_rate[lof.in.missense],
	                                      rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
	                                      annotation_phred=Anno.Int.PHRED.sub.category)
	Geno <- Genotype_sp$Geno
	Anno.Int.PHRED.sub.category <- Genotype_sp$annotation_phred
	rm(Genotype_sp)
	gc()

	results_list_missense <- rep(list(c()), length(obj_nullmodel_list))
	phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))

	if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
	{
	  for(num in 1:length(obj_nullmodel_list))
	  {
	    id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
	    id.phenotype.num.merge <- data.frame(id.phenotype.num)
	    id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
	    phenotype.id.match <- id.phenotype.num.merge$index
	    samplesize.num <- length(phenotype.id.match)
	    
	    ## Extract the genotype matrix and variant information of the current trait
	    Geno.num <- Geno[phenotype.id.match, , drop = FALSE]
	    MAC.in <- Matrix::colSums(Geno.num,na.rm = TRUE)
	    Missing_num.in <- Missing_num.sp(Geno.num)
	    MAF.in <- MAC.in/(2*(samplesize.num-Missing_num.in))
	    
	    if (geno_missing_imputation == "mean")
	    {
	      Geno.num <- na.replace.sp(Geno.num,m=2*MAF.in)
	    }
	    if (geno_missing_imputation == "minor")
	    {
	      Geno.num <- na.replace.sp(Geno.num,is_NA_to_Zero=TRUE)
	      MAF.in <- MAC.in/(2*samplesize.num)
	    }
	    
	    ## number of traits in analysis
	    n_pheno <- obj_nullmodel_list[[num]]$n.pheno
	    
	    ## SPA status
	    if(!is.null(obj_nullmodel_list[[num]]$use_SPA))
	    {
	      use_SPA <- obj_nullmodel_list[[num]]$use_SPA
	    }else
	    {
	      use_SPA <- FALSE
	    }
	    
	    pvalues <- 0
	    if(n_pheno == 1)
	    {
	      if(!use_SPA)
	      {
	        try(pvalues <- STAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
	      }else{
	        try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
	      }
	    }else
	    {
	      try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
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
	      results_list_missense[[num]] <- results
	    }
	  }
	}

	#################################################
	#         disruptive missense
	#################################################
	lof.in.dmissense <- (GENCODE.EXONIC.Category=="nonsynonymous SNV")&(MetaSVM_pred=="D")
	variant.id.gene.category <- variant.id.gene[lof.in.dmissense]

	## Annotation
	Anno.Int.PHRED.sub.category <- Anno.Int.PHRED.sub[lof.in.dmissense,]
	
	Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.id.gene.category,
	                                      sample.id=phenotype.id,
	                                      REF_AF=REF_AF[lof.in.dmissense],variant_maf_cutoff_filter=variant_maf_cutoff_filter,
	                                      Missing_rate=Missing_rate[lof.in.dmissense],
	                                      rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
	                                      annotation_phred=Anno.Int.PHRED.sub.category)
	Geno <- Genotype_sp$Geno
	Anno.Int.PHRED.sub.category <- Genotype_sp$annotation_phred
	rm(Genotype_sp)
	gc()

	results_list_disruptive_missense <- rep(list(c()), length(obj_nullmodel_list))
	phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))

	for (num in 1:length(obj_nullmodel_list))
	{
	  id.phenotype.num <- as.character(obj_nullmodel_list[[num]]$id_include)
	  id.phenotype.num.merge <- data.frame(id.phenotype.num)
	  id.phenotype.num.merge <- dplyr::left_join(id.phenotype.num.merge,phenotype.id.merge,by = c("id.phenotype.num" = "phenotype.id"))
	  phenotype.id.match <- id.phenotype.num.merge$index
	  samplesize.num <- length(phenotype.id.match)
	  
	  if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
	  {
	    ## Extract the genotype matrix and variant information of the current trait
	    Geno.num <- Geno[phenotype.id.match, , drop = FALSE]
	    MAC.in <- Matrix::colSums(Geno.num,na.rm = TRUE)
	    Missing_num.in <- Missing_num.sp(Geno.num)
	    MAF.in <- MAC.in/(2*(samplesize.num-Missing_num.in))
	    
	    if (geno_missing_imputation == "mean")
	    {
	      Geno.num <- na.replace.sp(Geno.num,m=2*MAF.in)
	    }
	    if (geno_missing_imputation == "minor")
	    {
	      Geno.num <- na.replace.sp(Geno.num,is_NA_to_Zero=TRUE)
	      MAF.in <- MAC.in/(2*samplesize.num)
	    }
	  } else
	  {
	    Geno.num <- NULL
	    MAF.in <- numeric()
	  }
	  
	  ## number of traits in analysis
	  n_pheno <- obj_nullmodel_list[[num]]$n.pheno
	  
	  ## SPA status
	  if(!is.null(obj_nullmodel_list[[num]]$use_SPA))
	  {
	    use_SPA <- obj_nullmodel_list[[num]]$use_SPA
	  }else
	  {
	    use_SPA <- FALSE
	  }
	  
	  pvalues <- 0
	  if(n_pheno == 1)
	  {
	    if(!use_SPA)
	    {
	      try(pvalues <- STAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
	    }else{
	      try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
	    }
	  }else
	  {
	    try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub.category,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
	  }
	  
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
	    
	    results <- rbind(results_list_missense[[num]],results_temp)
	  }else
	  {
	    results <- results_list_missense[[num]]
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
	          results_missense <- results
	          results_ds <- c()
	        }else{
	          results <- cbind(results,matrix(1,1,2))
	          colnames(results)[(dim(results)[2]-1):dim(results)[2]] <- c("Burden(1,25)-Disruptive","Burden(1,1)-Disruptive")
	          results_missense <- results
	          results_ds <- c()
	        }
	        
	        results_list_missense[[num]] <- results
	      }else
	      {
	        results_missense <- c()
	        results_ds <- results
	        results <- c()
	        results_list_disruptive_missense[[num]] <- results_ds
	      }
	    }
	    
	    if(!is.null(results))
	    {
	      if(dim(results)[1]==2)
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
	          
	          results_ds <- c()
	          results_ds <- rbind(results_ds,results[2,])
	          
	          results <- c()
	          results <- rbind(results,results_m)
	          
	          results_list_missense[[num]] <- results
	          results_list_disruptive_missense[[num]] <- results_ds
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
	          pvalues_sub <- as.numeric(results_m[6:length(results_m)][p_seq])
	          if(sum(is.na(pvalues_sub))>0)
	          {
	            if(sum(is.na(pvalues_sub))==length(pvalues_sub))
	            {
	              results_m["STAAR-B"] <- 1
	            }else
	            {
	              ## not all NAs
	              pvalues_sub <- pvalues_sub[!is.na(pvalues_sub)]
	              if(sum(pvalues_sub[pvalues_sub<1])>0)
	              {
	                ## not all ones
	                results_m["STAAR-B"] <- CCT(pvalues_sub[pvalues_sub<1])
	                
	              }else
	              {
	                results_m["STAAR-B"] <- 1
	                
	              }
	            }
	          }else
	          {
	            if(sum(pvalues_sub[pvalues_sub<1])>0)
	            {
	              results_m["STAAR-B"] <- CCT(pvalues_sub[pvalues_sub<1])
	            }else
	            {
	              results_m["STAAR-B"] <- 1
	            }
	          }
	          
	          ## calculate STAAR-B(1,25)
	          pvalues_sub <- as.numeric(results_m[6:length(results_m)][c(1:apc_num,(length(results_m)-6))])
	          if(sum(is.na(pvalues_sub))>0)
	          {
	            if(sum(is.na(pvalues_sub))==length(pvalues_sub))
	            {
	              results_m["STAAR-B(1,25)"] <- 1
	            }else
	            {
	              ## not all NAs
	              pvalues_sub <- pvalues_sub[!is.na(pvalues_sub)]
	              if(sum(pvalues_sub[pvalues_sub<1])>0)
	              {
	                ## not all ones
	                results_m["STAAR-B(1,25)"] <- CCT(pvalues_sub[pvalues_sub<1])
	                
	              }else
	              {
	                results_m["STAAR-B(1,25)"] <- 1
	                
	              }
	            }
	          }else
	          {
	            if(sum(pvalues_sub[pvalues_sub<1])>0)
	            {
	              results_m["STAAR-B(1,25)"] <- CCT(pvalues_sub[pvalues_sub<1])
	            }else
	            {
	              results_m["STAAR-B(1,25)"] <- 1
	            }
	          }
	          
	          ## calculate STAAR-B(1,1)
	          pvalues_sub <- as.numeric(results_m[6:length(results_m)][c(1:apc_num+(apc_num+1),(length(results_m)-5))])
	          if(sum(is.na(pvalues_sub))>0)
	          {
	            if(sum(is.na(pvalues_sub))==length(pvalues_sub))
	            {
	              results_m["STAAR-B(1,1)"] <- 1
	            }else
	            {
	              ## not all NAs
	              pvalues_sub <- pvalues_sub[!is.na(pvalues_sub)]
	              if(sum(pvalues_sub[pvalues_sub<1])>0)
	              {
	                ## not all ones
	                results_m["STAAR-B(1,1)"] <- CCT(pvalues_sub[pvalues_sub<1])
	                
	              }else
	              {
	                results_m["STAAR-B(1,1)"] <- 1
	                
	              }
	            }
	          }else
	          {
	            if(sum(pvalues_sub[pvalues_sub<1])>0)
	            {
	              results_m["STAAR-B(1,1)"] <- CCT(pvalues_sub[pvalues_sub<1])
	            }else
	            {
	              results_m["STAAR-B(1,1)"] <- 1
	            }
	          }
	          
	          results_ds <- c()
	          results_ds <- rbind(results_ds,results[2,])
	          
	          results <- c()
	          results <- rbind(results,results_m)
	          
	          results_list_missense[[num]] <- results
	          results_list_disruptive_missense[[num]] <- results_ds
	        }
	      }
	    }
	  }else
	  {
	    results <- c()
	    results_ds <- c()
	  }
	}

	results_coding_list <- list(plof = results_list_plof,
	                            plof_ds = results_list_plof_ds,
	                            missense = results_list_missense,
	                            disruptive_missense = results_list_disruptive_missense,
	                            synonymous = results_list_synonymous,
	                            ptv = results_list_ptv,
	                            ptv_ds = results_list_ptv_ds)

	seqResetFilter(genofile)

	return(results_coding_list)
}

