spatialde <- NULL
.onLoad <- function(libname, pkgname) {
  if (!reticulate::py_module_available("SpatialDE")) {
    message("Python module SpatialDE not available.\n
          Please use the terminal to install Python module SpatialDE:\n
          reticulate::py_install(\"spatialde\", pip = TRUE)\n")
  } else {
    spatialde <<- reticulate::import("SpatialDE", delay_load = TRUE)
  }
}

#' Perform SpatialDE test
#'
#' Test which gene is dependent on spatial location in the tissue. For more information,
#' see the original paper
#' \href{https://doi.org/10.1038/nmeth.4636}{SpatialDE: identification of spatially variable genes}.
#'
#' @param exp_mat Data frame or matrix with gene expression at each cell/loocation.
#' The data is assumed to be properly normalized. Genes can be in either columns or
#' rows. If the genes are in columns (and samples/cells in rows), set \code{row_sample = TRUE}.
#' Alternatively, a \code{SingleCellExperiment} or a \code{seurat} object can be used.
#' If it's a \code{SingleCellExperiment} object, it must have spatial coordinates in
#' \code{colData}. If it's a \code{seurat} object, it must have spatial coordinates
#' in the \code{meta.data} slot.
#' @param X Data frame or matrix with spatial coordinates. All columns must be numeric,
#' referring to spatial coordinates. This argument is only needed when \code{exp_mat}
#' is a matrix or data frame.
#' @param covariates Covariates that are not genes to which the test is also run,
#' as a comparison to genes. For example, log total counts in cells. This should
#' be a character vector, as column names of \code{colData} in \code{SingleCellExperiment}
#' objects, or column names of \code{meta.data} in \code{seurat} objects.
#' @param kernel_space The grid of covariance matrices to search over for the alternative
#' model can be specified using the kernel_space paramter. This should be a named list,
#' each name specifying a kernel, and each element a numeric vector specifying
#' different parameters for that kernel to be tried.
#' @param row_sample Whether samples/cells are in rows and genes are in columns. Defaults
#' to \code{FALSE}. This argument is only needed when \code{exp_mat}
#' is a matrix or data frame.
#' @param location_names Column names of column data (for \code{SingleCellExperiment}) or names of
#' columns of \code{meta.data} (for \code{seurat}) that specify the spatial coordinates
#' of cells.
#'
#' @return A data frame with the following columns:
#' \describe{
#' \item{\code{FSV}}{Fraction of total variance explained by spatial variance}
#' \item{\code{M}}{Internal code for model used}
#' \item{\code{g}}{Gene names}
#' \item{\code{l}}{Length scale the gene changes expression over}
#' \item{\code{max_delta}}{Maximum likelihood estimate of ratio between variance
#'  of non-spatial variance
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
#' @importFrom zeallot %<-%
#' @export
#'
RunSpatialDE <- function(exp_mat, ...) {
  UseMethod("RunSpatialDE", exp_mat)
}

#' @describeIn RunSpatialDE Perform SpatialDE test on matrices
#' @export RunSpatialDE.default
#' @method RunSpatialDE default
RunSpatialDE.default <- function(exp_mat, X, kernel_space = NULL, row_sample = FALSE) {
  ensure_load()
  c(X, exp_mat) %<-% convert2df(X, exp_mat, row_sample)
  if (nrow(X) != nrow(exp_mat)) {
    stop("X and exp_tab must have the same number of rows.")
  }
  spatialde$run(X, exp_mat, kernel_space)
}

#' @describeIn RunSpatialDE Perform SpatialDE test on \code{SingleCellExperiment} objects
#' @export RunSpatialDE.SingleCellExperiment
#' @method RunSpatialDE SingleCellExperiment
RunSpatialDE.SingleCellExperiment <- function(exp_mat, covariates = NULL,
                                              location_names = c("x", "y"),
                                              kernel_space = NULL) {
  c(m, X) %<-% sce2mat(exp_mat, covariates, location_names)
  RunSpatialDE.default(m, X, kernel_space)
}

#' @describeIn RunSpatialDE Perform SpatialDE test on \code{seurat} objects
#' @export RunSpatialDE.seurat
#' @method RunSpatialDE seurat
RunSpatialDE.seurat <- function(exp_mat, covariates = NULL,
                                location_names = c("x", "y"),
                                kernel_space = NULL) {
  c(m, X) %<-% seu2mat(exp_mat, covariates, location_names)
  RunSpatialDE.default(m, X, kernel_space)
}

#' Compare Model Fits with Different Models
#'
#' This way DE genes can be classified to interpretable function classes.
#' The strategy is based on ABCD in the Automatic Statistician, but
#' using precomputed covariance matrices for all models in the search space.
#' By default searches a grid of periodic covariance matrices and a linear
#' covariance matrix.
#'
#' @inheritParams RunSpatialDE
#' @param DE_mll_results The results from SpatialDE, only with rows for statistically
#' significant genes.
#'
#' @return A data frame
#' @export
ModelSearch <- function(exp_mat, ...) {
  UseMethod("ModelSearch", exp_mat)
}

#' @describeIn ModelSearch Compare models with matrix or data frame input
#' @export ModelSearch.default
#' @method ModelSearch default
ModelSearch.default <- function(exp_mat, X, DE_mll_results, kernel_space = NULL,
                        row_sample = FALSE) {
  ensure_load()
  c(X, exp_mat) %<-% convert2df(X, exp_mat, row_sample)
  spatialde$model_search(X, exp_mat, DE_mll_results, kernel_space)
}

#' @describeIn ModelSearch Compare models with \code{SingleCellExperiment} input
#' @export ModelSearch.SingleCellExperiment
#' @method ModelSearch SingleCellExperiment
ModelSearch.SingleCellExperiment <- function(exp_mat, DE_mll_results,
                                             covariates = NULL,
                                             location_names = c("x", "y"),
                                             kernel_space = NULL) {
  c(m, X) %<-% sce2mat(exp_mat, covariates, location_names)
  ModelSearch.default(m, X, DE_mll_results, kernel_space)
}

#' @describeIn ModelSearch Compare models with \code{seurat} input
#' @export ModelSearch.seurat
#' @method ModelSearch seurat
ModelSearch.seurat <- function(exp_mat, DE_mll_results, covariates = NULL,
                               location_names = c("x", "y"),
                               kernel_space = NULL) {
  c(m, X) %<-% seu2mat(exp_mat, covariates, location_names)
  ModelSearch.default(m, X, DE_mll_results, kernel_space)
}

#' Automatic Expression Histology (AEH)
#'
#' Classify the patterns of spatially differentially expressed genes.
#'
#' @inheritParams ModelSearch
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
AEHSpatialPatterns <- function(exp_mat, ...) {
  UseMethod("AEHSpatialPatterns", exp_mat)
}

#' @describeIn AEHSpatialPatterns AEH with matrix or data frame input
#' @export AEHSpatialPatterns.default
#' @method AEHSpatialPatterns default
AEHSpatialPatterns.default <- function(exp_mat, X, DE_mll_results, C, l,
                               row_sample = FALSE, ...) {
  ensure_load()
  c(X, exp_mat) %<-% convert2df(X, exp_mat, row_sample)
  res <- spatialde$aeh$spatial_patterns(X, exp_mat, DE_mll_results, C, l, ...)
  # Convert pandas data frame to R data frame
  res <- lapply(res, reticulate::py_to_r)
  res
}

#' @describeIn AEHSpatialPatterns AEH with \code{SingleCellExperiment} input
#' @export AEHSpatialPatterns.SingleCellExperiment
#' @method AEHSpatialPatterns SingleCellExperiment
AEHSpatialPatterns.SingleCellExperiment <- function(exp_mat, DE_mll_results, C, l,
                                                    covariates = NULL,
                                                    location_names = c("x", "y"),
                                                    kernel_space = NULL, ...) {
  c(m, X) %<-% sce2mat(exp_mat, covariates, location_names)
  AEHSpatialPatterns.default(m, X, DE_mll_results, C, l, ...)
}

#' @describeIn AEHSpatialPatterns AEH with \code{seurat} input
#' @export AEHSpatialPatterns.seurat
#' @method AEHSpatialPatterns seurat
AEHSpatialPatterns.seurat <- function(exp_mat, DE_mll_results, C, l,
                                      covariates = NULL,
                                      location_names = c("x", "y"),
                                      kernel_space = NULL, ...) {
  c(m, X) %<-% seu2mat(exp_mat, covariates, location_names)
  AEHSpatialPatterns.default(m, X, DE_mll_results, C, l, ...)
}
