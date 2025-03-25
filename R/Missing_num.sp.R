#' Count Missing Values in a Sparse Genotype Matrix
#'
#' The \code{Missing_num.sp} function counts the number of missing values (NA) 
#' in each column of a sparse genotype matrix (`dgCMatrix` format).
#'
#' @param genotype_sp A sparse genotype matrix of class `dgCMatrix` from the \code{\link{Matrix}} package.
#'
#' @return A numeric vector of length equal to the number of columns in \code{genotype_sp}, 
#' where each element represents the count of missing values in the corresponding column.
#'
#' @examples
#' library(Matrix)
#' set.seed(123)
#' # Create a sparse matrix with some NA values
#' mat <- Matrix(c(1, NA, 3, 0, NA, 2, 4, 5, NA), nrow = 3, sparse = TRUE)
#' print(mat)
#'
#' # Count missing values in each column
#' missing_counts <- Missing_num.sp(mat)
#' print(missing_counts)
#'
#' @export

Missing_num.sp <- function(genotype_sp) 
{
  if (!inherits(genotype_sp, "dgCMatrix")) stop("genotype_sp must be a dgCMatrix")
  
  x <- genotype_sp@x
  i <- genotype_sp@i
  p <- genotype_sp@p
  ncol <- ncol(genotype_sp)
  
  nas <- is.na(x)
  
  if (any(nas)) 
  {
    # Compute missing number for each column
    missing_num_col <- sapply(seq_len(ncol), function(j) 
    {
      if (p[j] < p[j + 1]) 
      {
        sum(nas[(p[j] + 1):p[j + 1]])  # Count NA values per column
      } else 
      {
        0
      }
    })
  } else
  {
    missing_num_col <- rep(0,ncol)
  }
  return(missing_num_col)
}
