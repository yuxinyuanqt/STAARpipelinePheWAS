#' Gene-centric PheWAS analysis of coding functional categories using STAAR procedure
#'
#' The \code{Gene_Centric_Coding_PheWAS} function takes in chromosome, gene name, functional category,
#' the object of opened annotated GDS file, and the list of objects from fitting the null models to analyze the association between a series of
#' quantitative/dichotomous phenotypes (including imbalanced case-control design) and coding functional categories of a gene by using STAAR procedure.
#' For each coding functional category, the STAAR-O p-value is a p-value from an omnibus test
#' that aggregated SKAT(1,25), SKAT(1,1), Burden(1,25), Burden(1,1), ACAT-V(1,25),
#' and ACAT-V(1,1) together with p-values of each test weighted by each annotation
#' using Cauchy method. For imbalance case-control setting, the results correspond to the STAAR-B p-value, which is a p-value from
#' an omnibus test that aggregated Burden(1,25) and Burden(1,1) together with p-values of each test weighted by each annotation using Cauchy method.
#' For multiple phenotype analysis (\code{obj_nullmodel$n.pheno > 1}),
#' the results correspond to multi-trait association p-values (e.g. MultiSTAAR-O) by leveraging
#' the correlation structure between multiple phenotypes.
#' @param chr chromosome.
#' @param gene_name name of the gene to be analyzed using STAAR procedure.
#' @param category the coding functional category to be analyzed using STAAR procedure. Choices include
#' \code{all_categories}, \code{plof}, \code{plof_ds}, \code{missense}, \code{disruptive_missense}, \code{synonymous},
#' \code{ptv}, \code{ptv_ds}, \code{all_categories_incl_ptv}  (default = \code{all_categories}).
#' @param genofile an object of opened annotated GDS (aGDS) file.
#' @param obj_nullmodel_list a list of objects from fitting the null models, which are either the output from \code{fit_nullmodel} function in the \code{STAARpipeline} package,
#' or the output from \code{fitNullModel} function in the \code{GENESIS} package and transformed using the \cr \code{genesis2staar_nullmodel} function in the \code{STAARpipeline} package.
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
#' @param silent logical: should the report of error messages be suppressed (default = FALSE).
#' @return A list of data frames containing the STAAR p-values (including STAAR-O or STAAR-B in imbalanced case-control setting) corresponding to each coding functional category of the given gene for each phenotype.
#' @references Li, Z., Li, X., et al. (2022). A framework for detecting
#' noncoding rare-variant associations of large-scale whole-genome sequencing
#' studies. \emph{Nature Methods}, \emph{19}(12), 1599-1611.
#' (\href{https://doi.org/10.1038/s41592-022-01640-x}{pub})
#' @references Li, X., Li, Z., et al. (2020). Dynamic incorporation of multiple
#' in silico functional annotations empowers rare variant association analysis of
#' large whole-genome sequencing studies at scale. \emph{Nature Genetics}, \emph{52}(9), 969-983.
#' (\href{https://doi.org/10.1038/s41588-020-0676-4}{pub})
#' @export

Gene_Centric_Coding_PheWAS <- function(chr,gene_name,category=c("all_categories","plof","plof_ds","missense","disruptive_missense","synonymous","ptv","ptv_ds","all_categories_incl_ptv"),
                                       genofile,obj_nullmodel_list,rare_maf_cutoff=0.01,rv_num_cutoff=2,rv_num_cutoff_max=1e9,rv_num_cutoff_max_prefilter=1e9,
                                       QC_label="annotation/filter",variant_type=c("SNV","Indel","variant"),geno_missing_imputation=c("mean","minor"),
                                       Annotation_dir="annotation/info/FunctionalAnnotation",Annotation_name_catalog,
                                       Use_annotation_weights=c(TRUE,FALSE),Annotation_name=NULL,
                                       SPA_p_filter=TRUE,p_filter_cutoff=0.05,silent=FALSE){

	## evaluate choices
	category <- match.arg(category)
	variant_type <- match.arg(variant_type)
	geno_missing_imputation <- match.arg(geno_missing_imputation)

	genes <- genes_info[genes_info[,2]==chr,]

	if(category=="all_categories")
	{
		results_list <- coding_PheWAS(chr,gene_name,genofile,obj_nullmodel_list,genes,
		                              rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
		                              rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
		                              QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
		                              Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
		                              Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
		                              SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
	}

	if(category=="plof")
	{
		results_list <- plof_PheWAS(chr,gene_name,genofile,obj_nullmodel_list,genes,
		                            rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
		                            rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
		                            QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
		                            Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
		                            Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
		                            SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
	}

	if(category=="plof_ds")
	{
		results_list <- plof_ds_PheWAS(chr,gene_name,genofile,obj_nullmodel_list,genes,
		                               rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
		                               rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
		                               QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
		                               Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
		                               Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
		                               SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
	}

	if(category=="missense")
	{
		results_list <- missense_PheWAS(chr,gene_name,genofile,obj_nullmodel_list,genes,
		                                rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
		                                rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
		                                QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
		                                Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
		                                Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
		                                SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
	}

	if(category=="disruptive_missense")
	{
		results_list <- disruptive_missense_PheWAS(chr,gene_name,genofile,obj_nullmodel_list,genes,
		                                           rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
		                                           rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
		                                           QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
		                                           Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
		                                           Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
		                                           SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
	}

	if(category=="synonymous")
	{
		results_list <- synonymous_PheWAS(chr,gene_name,genofile,obj_nullmodel_list,genes,
		                                  rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
		                                  rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
		                                  QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
		                                  Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
		                                  Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
		                                  SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
	}


	if(category=="ptv")
	{
		results_list <- ptv_PheWAS(chr,gene_name,genofile,obj_nullmodel_list,genes,
		                           rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
		                           rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
		                           QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
		                           Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
		                           Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
		                           SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
	}


	if(category=="ptv_ds")
	{
		results_list <- ptv_ds_PheWAS(chr,gene_name,genofile,obj_nullmodel_list,genes,
		                              rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
		                              rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
		                              QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
		                              Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
		                              Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
		                              SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
	}


	if(category=="all_categories_incl_ptv")
	{
		results_list <- coding_incl_ptv_PheWAS(chr,gene_name,genofile,obj_nullmodel_list,genes,
		                                       rare_maf_cutoff=rare_maf_cutoff,rv_num_cutoff=rv_num_cutoff,
		                                       rv_num_cutoff_max=rv_num_cutoff_max,rv_num_cutoff_max_prefilter=rv_num_cutoff_max_prefilter,
		                                       QC_label=QC_label,variant_type=variant_type,geno_missing_imputation=geno_missing_imputation,
		                                       Annotation_dir=Annotation_dir,Annotation_name_catalog=Annotation_name_catalog,
		                                       Use_annotation_weights=Use_annotation_weights,Annotation_name=Annotation_name,
		                                       SPA_p_filter=SPA_p_filter,p_filter_cutoff=p_filter_cutoff,silent=silent)
	}

	return(results_list)
}

