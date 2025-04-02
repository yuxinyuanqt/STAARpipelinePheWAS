synonymous_PheWAS <- function(chr,gene_name,genofile,obj_nullmodel_list,genes,
                              rare_maf_cutoff=0.01,rv_num_cutoff=2,rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                              QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                              Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                              Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                              SPA_p_filter=FALSE,p_filter_cutoff=0.05,silent=FALSE){

	## evaluate choices
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	phenotype.id <- Reduce(union, lapply(obj_nullmodel_list, function(x) {
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

	seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

	### synonymous
	## Gencode_Exonic
	GENCODE.EXONIC.Category <- seqGetData(genofile, paste0(Annotation_dir,Annotation_name_catalog$dir[which(Annotation_name_catalog$name=="GENCODE.EXONIC.Category")]))

	variant.id.gene <- seqGetData(genofile, "variant.id")
	lof.in.synonymous <- (GENCODE.EXONIC.Category=="synonymous SNV")
	variant.id.gene <- variant.id.gene[lof.in.synonymous]
	
	if(length(variant.id.gene)==0)
	{
	  results_list <- rep(list(c()), length(obj_nullmodel_list))
	  seqResetFilter(genofile)
	  return(results_list)
	}

	seqSetFilter(genofile,variant.id=variant.id.gene,sample.id=phenotype.id)

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
	## variant id
	variant.id.gene <- seqGetData(genofile, "variant.id")
	variant_maf_cutoff_filter <- ifelse(rare_maf_cutoff<=0.01,0.05,1)
	seqResetFilter(genofile)
	
	Genotype_sp <- Genotype_sp_extraction(genofile,variant.id=variant.id.gene,
	                                      sample.id=phenotype.id,
	                                      REF_AF=REF_AF,variant_maf_cutoff_filter=variant_maf_cutoff_filter,
	                                      Missing_rate=Missing_rate,
	                                      rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
	                                      annotation_phred=Anno.Int.PHRED.sub)
	Geno <- Genotype_sp$Geno
	Anno.Int.PHRED.sub <- Genotype_sp$annotation_phred
	rm(Genotype_sp)
	gc()

	results_list <- rep(list(c()), length(obj_nullmodel_list))
	phenotype.id.merge <- data.frame(phenotype.id, index = seq(1, length(phenotype.id)))

	if(!is.null(Geno) & inherits(Geno, "dgCMatrix"))
	{
	  for (num in 1:length(obj_nullmodel_list))
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
	        try(pvalues <- STAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
	      }else{
	        try(pvalues <- STAAR_Binary_SPA_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max,SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff),silent=silent)
	      }
	    }else
	    {
	      try(pvalues <- MultiSTAAR_sp(Geno.num,MAF.in,obj_nullmodel_list[[num]],Anno.Int.PHRED.sub,rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,rv_num_cutoff_max=rv_num_cutoff_max),silent=silent)
	    }
	    
	    results <- c()
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
	      
	      results_list[[num]] <- results
	    }
	  }
	}

	seqResetFilter(genofile)

	return(results_list)
}

