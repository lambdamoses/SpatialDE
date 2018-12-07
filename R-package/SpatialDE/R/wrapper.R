spatialde <- NULL
.onLoad <- function(libname, pkgname) {
  if (!reticulate::py_module_available("SpatialDE")) {
    cat("Please install the Python module SpatialDE")
    #cat("Installing Python module SpatialDE.", "\n")
    #reticulate::py_install("spatialde")
  }
  spatialde <<- reticulate::import("SpatialDE", delay_load = TRUE)
}

#' Perform SpatialDE test
#'
#' Test which gene is dependent on spatial location in the tissue. For more information,
#' see the original paper
#' \href{https://doi.org/10.1038/nmeth.4636}{SpatialDE: identification of spatially variable genes}.
#'
#' @param X Data frame or matrix with spatial coordinates.
#' @param exp_mat Data frame or matrix with gene expression at each cell/loocation, with
#' genes in columns. The data is assumed to be properly normalized.
#' @param kernel_space The grid of covariance matrices to search over for the alternative
#' model can be specified using the kernel_space paramter. This should be a named list,
#' each name specifying a kernel, and each element a numeric vector specifying
#' different parameters for that kernel to be tried.
#'
#' @return A data frame with the following columns:
#' \describe{
#' \item{\code{FSV}}{Fraction of total variance explained by spatial variance}
#' \item{\code{M}}{Internal code for model used}
#' \item{\code{g}}{Gene names}
#' \item{\code{l}}{Length scale the gene changes expression over}
#' \item{\code{max_delta}}{Maximum likelihood estimate of ratio between variance of non-spatial variance
#' and spatial variance}
#' \item{\code{max_ll}}{Maximum log likelihood}
#' \item{\code{max_mu_hat}}{Mean expression level for the gene}
#' \item{\code{max_s2_t_hat}}{}
#' \item{\code{model}}{The model used to fit the data, which indicates spatial pattern
#' of gene expression}
#' \item{\code{n}}{}
#' \item{\code{s2_FSV}}{Standard error of FSV}
#' \item{\code{s2_logdelta}}{Standard error of log delta}
#' \item{\code{time}}{Time taken on this gene}
#' \item{\code{BIC}}{Akaike's An Information Criterion, with k = log(n); a smaller number
#' means the model fits the data better}
#' \item{\code{max_ll_null}}{Maximum likelihood of null model}
#' \item{\code{LLR}}{Log likelihood ratio}
#' \item{\code{pval}}{P-value of spatial differential expression}
#' \item{\code{qval}}{Significance after correcting for multiple testing}}
#' @export
run <- function(X, exp_mat, kernel_space = NULL) {
  if (is.matrix(X)) X <- as.data.frame(X)
  if (is.matrix(exp_mat)) exp_mat <- as.data.frame(exp_mat)
  if (nrow(X) != nrow(exp_mat)) {
    stop("X and exp_tab must have the same number of rows.")
  }
  spatialde$run(X, exp_mat, kernel_space)
}

#' Automatic Expression Histology (AEH)
#'
#' Classify the patterns of spatially differentially expressed genes.
#'
#' @param X Data frame or matrix with spatial coordinates.
#' @param exp_mat Data frame or matrix with gene expression at each cell/loocation, with
#' genes in columns. The data is assumed to be properly normalized.
#' @param DE_mll_results The results from SpatialDE, only with rows for statistically
#' significant genes.
#' @param C Number of expected patterns.
#' @param l Length scale for the patterns.
#' @param ... Extra arguments passed to the function \code{fit_patterns}.
#' See \href{https://github.com/Teichlab/SpatialDE/blob/23ec19131473bf87d5fb3b17508aae430f2e6310/Python-module/SpatialDE/aeh.py#L137}{the Python module this package wraps}
#' for the available arguments.
#'
#' @return A list of 2 data frames, one with pattern membership information
#' for each gene, and the other with the posterior mean underlying expression for genes in
#' given spatial patterns.
#' @export
aeh_spatial_patterns <- function(X, exp_mat, DE_mll_results, C, l, ...) {
  if (is.matrix(X)) X <- as.data.frame(X)
  if (is.matrix(exp_mat)) exp_mat <- as.data.frame(exp_mat)
  res <- spatialde$aeh$spatial_patterns(X, exp_mat, DE_mll_results, C, l, ...)
  # Convert pandas data frame to R data frame
  res <- lapply(res, reticulate::py_to_r)
  res
}
