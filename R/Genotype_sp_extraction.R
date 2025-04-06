#' Genotype extraction and filtering for association analysis
#'
#' The \code{Genotype_sp_extraction} function extracts genotype data for a given set of variants and samples from a GDS file.
#' It applies multiple filters based on allele frequency, missingness rate, and user-defined thresholds to categorize variants into different dosage groups.
#' The function returns a sparse genotype matrix (\code{"dgCMatrix"}), variant-level information, and optional annotation data.
#'
#' The extracted variants are processed based on the following three cases:
#' (1) Case 1: \code{ALT_AF > 0.5}
#'     - Use \code{"$dosage"} to extract Geno (dosages of the reference allele).
#'     - Convert the extracted data into the \code{"dgCMatrix"} format (sparse matrix).
#'
#' (2) Case 2: \code{ALT_AF ≤ 0.5} and (\code{MAF ≥ rare_maf_cutoff} or \code{Missing_rate ≥ Missing_cutoff})
#'     - Use \code{"$dosage_alt"} to extract Geno (dosages of the alternative allele).
#'     - Convert the extracted data into the \code{"dgCMatrix"} format (sparse matrix).
#'
#' (3) Case 3: \code{ALT_AF ≤ 0.5} and (\code{MAF < rare_maf_cutoff} and \code{Missing_rate < Missing_cutoff})
#'     - Use \code{"$dosage_ap"} to directly extract Geno in the \code{"dgCMatrix"} format.
#'
#' Note:
#' - \code{REF_AF} and \code{Missing_rate} can be efficiently computed using
#'   \code{SeqArray::seqGetAF_AC_Missing(genofile, minor=FALSE)}.
#'
#' @param genofile an object of opened annotated GDS (aGDS) file.
#' @param variant.id 	ID of selected variants.
#' @param sample.id ID of selected samples.
#' @param REF_AF a numeric vector of reference allele frequencies for each variant.
#' @param rare_maf_cutoff the cutoff of maximum minor allele frequency in defining rare variants (default = 0.01).
#' @param variant_maf_cutoff_filter a numeric value specifying the MAF threshold for excluding common variants in variant-set analysis. Default is 1 for individual-variant analysis.
#' @param Missing_rate a numeric vector of missing rates for each variant.
#' @param Missing_cutoff a numeric value specifying the maximum missing rate threshold for defining high-missing variants (default = 0.01).
#' @param subset_variants_num the number of variants to extract per subset for each time in Case 1 and Case 2 (default = 2e3).
#' @param rv_num_cutoff_max_prefilter the cutoff of maximum number of variants before extracting the genotype matrix (default = 1e+09).
#' @param annotation_phred a optional data frame or matrix of functional annotation data of dimension p*q (or a vector of a single annotation score with length p). See \code{STAAR::STAAR} for more details.
#'
#' @return A list containing:
#' \item{Geno}{A sparse matrix of genotypes.}
#' \item{results_information}{A data frame with variant-level details, including chromosome, position, reference and alternative alleles, MAF, ALT_AF, missing rate, and variant ID.}
#' \item{annotation_phred}{A data frame containing filtered functional annotation data if provided; otherwise, NULL.}
#'
#' @export

Genotype_sp_extraction <- function(genofile,variant.id,sample.id,
                                   REF_AF,rare_maf_cutoff=0.01,
                                   variant_maf_cutoff_filter=1,
                                   Missing_rate,Missing_cutoff=0.01,
                                   subset_variants_num=1e3,
                                   rv_num_cutoff_max_prefilter=1e9,
                                   annotation_phred=NULL)
{
  
  if((length(variant.id) != length(REF_AF)) | (length(Missing_rate) != length(REF_AF)))
  {
    stop(paste0("Lengths don't match for variant.id, REF_AF and Missing_rate!"))
  }
  
  annotation_phred <- as.data.frame(annotation_phred)
  if(dim(annotation_phred)[1] != 0 & length(variant.id) != dim(annotation_phred)[1])
  {
    stop(paste0("Dimensions don't match for genotype and annotation!"))
  }
  
  ## get AF
  ALT_AF <- 1-REF_AF
  MAF <- ifelse(REF_AF>=ALT_AF,ALT_AF,REF_AF)
  is.include <- !((MAF>=variant_maf_cutoff_filter) | (MAF==0) | is.na(MAF) | is.na(Missing_rate))
  
  variant.id.in <- variant.id[is.include]
  ALT_AF.in <- ALT_AF[is.include]
  MAF.in <- MAF[is.include]
  Missing_rate.in <- Missing_rate[is.include]
  
  annotation_phred.sub <- annotation_phred[is.include,,drop=FALSE]
  
  results <- list(Geno=NULL,results_information=NULL,annotation_phred=NULL)
  
  if(length(variant.id.in) == 0 | length(variant.id.in)>=rv_num_cutoff_max_prefilter)
  {
    return(results)
  }
  
  is.dosge <- (ALT_AF.in>0.5)
  is.dosge_alt <- (ALT_AF.in<=0.5) & (MAF.in>=rare_maf_cutoff | Missing_rate.in>=Missing_cutoff)
  is.dosge_sp <- (ALT_AF.in<=0.5) & (MAF.in<rare_maf_cutoff) & (Missing_rate.in<Missing_cutoff)
  
  ## dosage
  if(sum(is.dosge)>=1)
  {
    variant.id.dosage <- variant.id.in[is.dosge]
    subset.num.dosage <- ceiling(length(variant.id.dosage)/subset_variants_num)
    
    MAF.in.dosage <- MAF.in[is.dosge]
    ALT_AF.in.dosage <- ALT_AF.in[is.dosge]
    Missing_rate.in.dosage <- Missing_rate.in[is.dosge]
    
    Geno_dosage <- NULL
    results_dosage <- NULL
    
    for(kk in 1:subset.num.dosage)
    {
      if(kk < subset.num.dosage)
      {
        is.in.dosage <- ((kk-1)*subset_variants_num+1):(kk*subset_variants_num)
      }
      if(kk == subset.num.dosage)
      {
        is.in.dosage <- ((kk-1)*subset_variants_num+1):length(variant.id.dosage)
      }
      
      seqSetFilter(genofile,variant.id=variant.id.dosage[is.in.dosage],sample.id=sample.id)
      
      ## genotype id
      id.genotype <- seqGetData(genofile,"sample.id")
      
      id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
      phenotype.id.merge <- data.frame(sample.id)
      phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("sample.id"="id.genotype"))
      id.genotype.match <- phenotype.id.merge$index
      
      Geno_dosage_temp <- seqGetData(genofile, "$dosage")
      Geno_dosage_temp <- as(Geno_dosage_temp,"dgCMatrix")
      Geno_dosage_temp <- Geno_dosage_temp[id.genotype.match,,drop=FALSE]
      Geno_dosage <- cbind(Geno_dosage,Geno_dosage_temp)
      gc()
      
      CHR.dosage <- as.numeric(seqGetData(genofile, "chromosome"))
      position.dosage <- as.numeric(seqGetData(genofile, "position"))
      REF.dosage <- as.character(seqGetData(genofile, "$ref"))
      ALT.dosage <- as.character(seqGetData(genofile, "$alt"))
      variant.id.in.dosage <- seqGetData(genofile, "variant.id")
      
      results_dosage_temp <- data.frame(CHR=CHR.dosage,position=position.dosage,
                                        REF=REF.dosage,ALT=ALT.dosage,
                                        variant.id=variant.id.in.dosage)
      results_dosage <- rbind(results_dosage,results_dosage_temp)
      
      seqResetFilter(genofile)
    }
    
    annotation_phred.sub_dosage <- annotation_phred.sub[is.dosge,,drop=FALSE]
    ## Reorder
    variant_order.dosage <- match(results_dosage$variant.id, variant.id.dosage)
    annotation_phred.sub_dosage <- annotation_phred.sub_dosage[variant_order.dosage,,drop=FALSE]
    results_dosage <- dplyr::left_join(results_dosage,
                                       data.frame(MAF=MAF.in.dosage,
                                                  ALT_AF=ALT_AF.in.dosage,
                                                  Missing_rate=Missing_rate.in.dosage,
                                                  variant.id=variant.id.dosage),
                                       by="variant.id")
    Geno <- Geno_dosage
    results_info <- results_dosage
    annotation_phred.sub_sort <- annotation_phred.sub_dosage
  }
  
  ## dosage_alt
  if(sum(is.dosge_alt)>=1)
  {
    variant.id.dosage_alt <- variant.id.in[is.dosge_alt]
    subset.num.dosage_alt <- ceiling(length(variant.id.dosage_alt)/subset_variants_num)
    
    MAF.in.dosage_alt <- MAF.in[is.dosge_alt]
    ALT_AF.in.dosage_alt <- ALT_AF.in[is.dosge_alt]
    Missing_rate.in.dosage_alt <- Missing_rate.in[is.dosge_alt]
    
    Geno_dosage_alt <- NULL
    results_dosage_alt <- NULL
    
    for(kk in 1:subset.num.dosage_alt)
    {
      if(kk < subset.num.dosage_alt)
      {
        is.in.dosage_alt <- ((kk-1)*subset_variants_num+1):(kk*subset_variants_num)
      }
      if(kk == subset.num.dosage_alt)
      {
        is.in.dosage_alt <- ((kk-1)*subset_variants_num+1):length(variant.id.dosage_alt)
      }
      
      seqSetFilter(genofile,variant.id=variant.id.dosage_alt[is.in.dosage_alt],sample.id=sample.id)
      
      ## genotype id
      id.genotype <- seqGetData(genofile,"sample.id")
      
      id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
      phenotype.id.merge <- data.frame(sample.id)
      phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("sample.id"="id.genotype"))
      id.genotype.match <- phenotype.id.merge$index
      
      Geno_dosage_alt_temp <- seqGetData(genofile, "$dosage_alt")
      Geno_dosage_alt_temp <- as(Geno_dosage_alt_temp,"dgCMatrix")
      Geno_dosage_alt_temp <- Geno_dosage_alt_temp[id.genotype.match,,drop=FALSE]
      Geno_dosage_alt <- cbind(Geno_dosage_alt,Geno_dosage_alt_temp)
      gc()
      
      CHR.dosage_alt <- as.numeric(seqGetData(genofile, "chromosome"))
      position.dosage_alt <- as.numeric(seqGetData(genofile, "position"))
      REF.dosage_alt <- as.character(seqGetData(genofile, "$ref"))
      ALT.dosage_alt <- as.character(seqGetData(genofile, "$alt"))
      variant.id.in.dosage_alt <- seqGetData(genofile, "variant.id")
      
      results_dosage_alt_temp <- data.frame(CHR=CHR.dosage_alt,position=position.dosage_alt,
                                            REF=REF.dosage_alt,ALT=ALT.dosage_alt,
                                            variant.id=variant.id.in.dosage_alt)
      results_dosage_alt <- rbind(results_dosage_alt,results_dosage_alt_temp)
      
      seqResetFilter(genofile)
    }
    
    annotation_phred.sub_alt <- annotation_phred.sub[is.dosge_alt,,drop=FALSE]
    ## Reorder
    variant_order.dosage_alt <- match(results_dosage_alt$variant.id, variant.id.dosage_alt)
    annotation_phred.sub_alt <- annotation_phred.sub_alt[variant_order.dosage_alt,,drop=FALSE]
    results_dosage_alt <- dplyr::left_join(results_dosage_alt,
                                           data.frame(MAF=MAF.in.dosage_alt,
                                                      ALT_AF=ALT_AF.in.dosage_alt,
                                                      Missing_rate=Missing_rate.in.dosage_alt,
                                                      variant.id=variant.id.dosage_alt),
                                           by="variant.id")
    
    if(sum(is.dosge)>=1)
    {
      Geno <- cbind(Geno_dosage_alt,Geno_dosage)
      results_info <- rbind(results_dosage_alt,results_dosage)
      annotation_phred.sub_sort <- rbind(annotation_phred.sub_alt,annotation_phred.sub_dosage)
    } else
    {
      Geno <- Geno_dosage_alt
      results_info <- results_dosage_alt
      annotation_phred.sub_sort <- annotation_phred.sub_alt
    }
  }
  
  ## dosage_sp
  if(sum(is.dosge_sp)>=1)
  {
    variant.id.dosage_sp <- variant.id.in[is.dosge_sp]
    
    MAF.in.dosage_sp <- MAF.in[is.dosge_sp]
    ALT_AF.in.dosage_sp <- ALT_AF.in[is.dosge_sp]
    Missing_rate.in.dosage_sp <- Missing_rate.in[is.dosge_sp]
    
    Geno_dosage_sp <- NULL
    results_dosage_sp <- NULL
    
    seqSetFilter(genofile,variant.id=variant.id.dosage_sp,sample.id=sample.id)
    
    ## genotype id
    id.genotype <- seqGetData(genofile,"sample.id")
    
    id.genotype.merge <- data.frame(id.genotype,index=seq(1,length(id.genotype)))
    phenotype.id.merge <- data.frame(sample.id)
    phenotype.id.merge <- dplyr::left_join(phenotype.id.merge,id.genotype.merge,by=c("sample.id"="id.genotype"))
    id.genotype.match <- phenotype.id.merge$index
    
    Geno_dosage_sp <- seqGetData(genofile, "$dosage_sp")
    Geno_dosage_sp <- Geno_dosage_sp[id.genotype.match,,drop=FALSE]
    
    CHR.dosage_sp <- as.numeric(seqGetData(genofile, "chromosome"))
    position.dosage_sp <- as.numeric(seqGetData(genofile, "position"))
    REF.dosage_sp <- as.character(seqGetData(genofile, "$ref"))
    ALT.dosage_sp <- as.character(seqGetData(genofile, "$alt"))
    variant.id.in.dosage_sp <- seqGetData(genofile, "variant.id")
    
    annotation_phred.sub_sp <- annotation_phred.sub[is.dosge_sp,,drop=FALSE]
    results_dosage_sp <- data.frame(CHR=CHR.dosage_sp,position=position.dosage_sp,
                                    REF=REF.dosage_sp,ALT=ALT.dosage_sp,
                                    variant.id=variant.id.in.dosage_sp)
    
    ## Reorder
    variant_order.dosage_sp <- match(results_dosage_sp$variant.id, variant.id.dosage_sp)
    annotation_phred.sub_sp <- annotation_phred.sub_sp[variant_order.dosage_sp,,drop=FALSE]
    results_dosage_sp <- dplyr::left_join(results_dosage_sp,
                                           data.frame(MAF=MAF.in.dosage_sp,
                                                      ALT_AF=ALT_AF.in.dosage_sp,
                                                      Missing_rate=Missing_rate.in.dosage_sp,
                                                      variant.id=variant.id.dosage_sp),
                                           by="variant.id")
    
    seqResetFilter(genofile)
    
    if(sum(is.dosge)>=1 | sum(is.dosge_alt)>=1)
    {
      Geno <- cbind(Geno_dosage_sp,Geno)
      results_info <- rbind(results_dosage_sp,results_info)
      annotation_phred.sub_sort <- rbind(annotation_phred.sub_sp,annotation_phred.sub_sort)
    } else
    {
      Geno <- Geno_dosage_sp
      results_info <- results_dosage_sp
      annotation_phred.sub_sort <- annotation_phred.sub_sp
    }
  }
  
  if(dim(annotation_phred.sub_sort)[2]==0)
  {
    annotation_phred.sub_sort <- NULL
  }
  
  results <- list(Geno=Geno,results_information=results_info,annotation_phred=annotation_phred.sub_sort)
  return(results)
}