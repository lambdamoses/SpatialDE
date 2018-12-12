#' Ensusre that the \code{SpatialDE} Python module is loaded
#'
#' By default, if the Python module \code{SpatialDE} is instatlled and \code{reticulate}
#' is calling the correct version of Python, then the Python module is loaded when
#' this R package is loaded. But the user might not have the Python module installed.
#' This function is to ensure that the Python module is loaded before it's called from
#' this R package.
#'
#' @return An error if the Python module is not installed. Otherwise loads the
#' Python module.
ensure_load <- function() {
  if (!reticulate::py_module_available("SpatialDE")) {
    stop("Python module SpatialDE not available.\n
          Please use the terminal to install Python module SpatialDE:\n
          reticulate::py_install(\"spatialde\", pip = TRUE)\n")
  } else if (!exists("spatialde", mode = "environment")) {
    spatialde <<- reticulate::import("SpatialDE", delay_load = TRUE)
  }
}

#' Convert inputs into the correct format
#'
#' The Python module takes in data frames as input, but you may pass a matrix to the
#' wrapper functions. This function converts whatever passed to the wrapper functions
#' into data frames.
#'
#' @param X Data frame or matrix with spatial coordinates. All columns must be numeric,
#' referring to spatial coordinates.
#' @param exp_mat Data frame or matrix with gene expression at each cell/loocation.
#' The data is assumed to be properly normalized. Genes can be in either columns or
#' rows. If the genes are in columns (and samples/cells in rows), set \code{row_sample = TRUE}.
#' @param row_sample Whether samples/cells are in rows and genes are in columns. Defaults
#' to \code{FALSE}.
#' @return A list of two data frames
convert2df <- function(X, exp_mat, row_sample = FALSE) {
  if (any(grepl("Matrix", class(X)))) X <- as.matrix(X)
  if (any(grepl("Matrix", class(exp_mat)))) exp_mat <- as.matrix(exp_mat)
  if (is.matrix(X)) X <- as.data.frame(X)
  if (!row_sample) {
    exp_mat <- t(exp_mat)
  }
  if (is.matrix(exp_mat)) exp_mat <- as.data.frame(exp_mat)
  list(X, exp_mat)
}

#' Convert \code{SingleCellExperiment} to matrix
#'
#' This function pulls out relevant parts of the \code{SingleCellExperiment} object
#' and converts it into matrix before calling the Python module `SpatialDE`.
#'
#' @param object The \code{SingleCellExperiment} object.
#' @param covariates Covariates that are not genes to which the test is also run,
#' as a comparison to genes. For example, log total counts in cells. This should
#' be a character vector, as column names of \code{colData} in \code{SingleCellExperiment}
#' objects, or column names of \code{meta.data} in \code{Seurat} objects.
#' @param location_names Column names of column data (for \code{SingleCellExperiment}) or names of
#' columns of \code{meta.data} (for \code{Seurat}) that specify the spatial coordinates
#' of cells.
#' @return A list with a matrix and a data frame
#' @importFrom SingleCellExperiment normcounts colData
sce2mat <- function(object, covariates, location_names) {
  m <- normcounts(object)
  if (!is.null(covariates)) {
    covar <- t(as.matrix(colData(object)[,covariates]))
    rns <- c(rownames(m), covariates)
    m <- rbind(m, covar)
    rownames(m) <- rns
  }
  X <- as.data.frame(colData(object))[,location_names]
  list(m, X)
}

#' Convert \code{Seurat} to matrix
#'
#' This function pulls out relevant parts of the \code{Seurat} object
#' and converts it into matrix before calling the Python module `SpatialDE`.
#'
#' @param object The \code{Seurat} object.
#' @param covariates Covariates that are not genes to which the test is also run,
#' as a comparison to genes. For example, log total counts in cells. This should
#' be a character vector, as column names of \code{colData} in \code{SingleCellExperiment}
#' objects, or column names of \code{meta.data} in \code{Seurat} objects.
#' @param location_names Column names of column data (for \code{SingleCellExperiment}) or names of
#' columns of \code{meta.data} (for \code{Seurat}) that specify the spatial coordinates
#' of cells.
#' @return A list with a matrix and a data frame
seu2mat <- function(object, covariates, location_names) {
  m <- object@data
  if (!is.null(covariates)) {
    covar <- t(as.matrix(object@meta.data[,covariates]))
    rns <- c(rownames(m), covariates)
    m <- rbind(m, covar)
    rownames(m) <- rns
  }
  X <- object@meta.data[,location_names]
  list(m, X)
}
