#' @include wrapper.R
NULL

#' Inverse log scale
#'
#' Custom transform for inverted log transformed axis in ggplot2.
#' @importFrom scales trans_new log_breaks
reverse_log10 <- function() {
  trans <- function(x) -log10(x)
  inv <- function(x) 10^(-x)
  scales::trans_new("reverse_log10", trans, inv, scales::log_breaks(base = 10),
                    domain = c(1e-100, Inf))
}

#' Plot Fraction Spatial Variance vs Q-value
#'
#' Optionally provide model selection results to the function will color points by model.
#' Point size corresponds to certinety of the FSV value.
#'
#' @param results Results from SpatialDE.
#' @param ms_results Model selection results, should be a data frame with columns
#' \code{g} for gene names and \code{model} for the model selected.
#' @param certain_only Only plot results with narrow 95\% confidence interval.
#' @param log_x Whether to display x axis in log scale.
#' @param do_label Display gene names for statistically significant genes, default to \code{TRUE}.
#' @param covariate_names Names of covariates as a reference, default to \code{NULL}.
#'
#' @return A \code{ggplot2} object.
#' @import ggplot2
#' @importFrom dplyr mutate filter select left_join case_when
#' @importFrom magrittr %>%
#' @export
FSV_sig <- function(results, ms_results = NULL, certain_only = FALSE, log_x = FALSE,
                    do_label = TRUE, covariate_names = NULL) {
  if (!is.null(ms_results)) {
    results <- results %>%
      full_join(ms_results[,c("g", "model")], by = "g", suffix = c("", "_bic"))
  } else {
    results <- results %>%
      rename(model_bic = model)
  }
  results <- results %>%
    mutate(FSV95conf = 2 * sqrt(s2_FSV),
           conf_categories = fct_rev(cut(FSV95conf, c(0, 0.1, 1, Inf))),
           is_covariate = FALSE,
           # More user friendly model labels
           color_categories = case_when(model_bic == "SE" ~ "general",
                                        model_bic == "PER" ~ "periodic",
                                        model_bic == "linear" ~ "linear"))
  if (!is.null(covariate_names)) {
    results <- results %>%
      mutate(is_covariate = g %in% covariate_names)
  }
  if (certain_only) {
    results <- results %>%
      filter(conf_categories == "(0,0.1]")
  }
  colors_use <- scales::hue_pal()(length(unique(results$model_bic)))
  p <- ggplot(results, aes(FSV, qval)) +
    geom_hline(yintercept = 0.05, linetype = 2) +
    scale_color_manual(values = colors_use, na.translate = TRUE, na.value = "black",
                       guide = guide_legend(title = "model")) +
    scale_y_continuous(trans = reverse_log10()) +
    annotate(geom = "text", x = 0, y = 0.05, label = "0.05")
  if (!is.null(covariate_names)) {
    p <- p +
      scale_shape_manual(values = c(16,4), guide = guide_legend(title = "covariate"))
  }
  if (!certain_only) {
    p <- p +
      geom_point(aes(color = color_categories, size = conf_categories, shape = is_covariate),
                 alpha = 0.5) +
      scale_size_manual(values = c(0.7, 1.5, 3),
                        guide = guide_legend(title = "confidence \n category"))
  } else {
    p <- p +
      geom_point(aes(color = color_categories, shape = is_covariate),
                 alpha = 0.5)
  }
  if (log_x) {
    p <- p +
      scale_x_log10()
  }
  if (do_label) {
    gene_label <- results %>%
      filter(qval < 0.05) %>%
      select(FSV, qval, g)
    p <- p + ggrepel::geom_label_repel(aes(label = g), data = gene_label)
  }
  p
}

#' Plot Spatial Patterns of Multiple Genes
#'
#' While \code{ggplot2} allows us to easily generate multi-faceted plots, we can't
#' use a separate color scale for each facet. This function allows you to make
#' multi-faceted plots with a separate color scale for each facet, so high expression
#' of one gene in the facets will not affect the dynamic range of colors for the other,
#' not so highly expressed genes, making it easier to discern spatial patterns.
#'
#' @inheritParams RunSpatialDE
#' @param genes_plot Character vector specifying which genes are to be plotted.
#' Covariate names are permitted in this argument.
#' @param facet_titles Title for each facet. This must be a character vector with the
#' same length as \code{genes_plot}.
#' @param viridis_option This function uses the \code{viridis} palette to color cells
#' for gene expression. Four options are available: "magma" (or "A"), "inferno" (or "B"),
#' "plasma" (or "C"), "viridis" (or "D", the default option) and "cividis" (or "E").
#' @param ncol Number of columns to arrange the plots.
#' @param dark_theme Whether dark background should be used; this is helpful to highlight
#' cells with high expression when using the \code{viridis} palette.
#' @return This function draws on the current divice without returning an object.
#' @importFrom gridExtra grid.arrange
#' @export
MultiGenePlot <- function(exp_mat, ...) {
  UseMethod("MultiGenePlot", exp_mat)
}

#' @describeIn MultiGenePlot Plot multiple genes for matrix or data frame input
#' @export MultiGenePlot.default
#' @method MultiGenePlot default
MultiGenePlot.default <- function(exp_mat, X, genes_plot, row_sample = FALSE,
                                  facet_titles = NULL, viridis_option = "D",
                                  ncol = 2, dark_theme = TRUE) {
  if (!is.null(facet_titles) && length(genes_plot) != length(facet_titles)) {
    stop("Arguments genes_plot and facet_titles must have the same length.")
  }
  if (is.matrix(X)) {
    X <- as.data.frame(X)
  }
  if (row_sample) {
    exp_mat <- t(exp_mat)
  }
  df <- cbind(X, t(exp_mat[genes_plot,]))
  pls <- lapply(seq_along(genes_plot),
                function(i) {
                  col <- paste0("`", genes_plot[i], "`")
                  p <- ggplot(df, aes_string("x", "y", color = col)) +
                    geom_point() +
                    coord_equal() +
                    scale_color_viridis_c(option = viridis_option)
                  if (!is.null(facet_titles)) {
                    p <- p +
                      ggtitle(label = facet_titles[i])
                  }
                  if (dark_theme) {
                    p <- p +
                      theme_dark()
                  }
                  p
                })
  grid.arrange(grobs = pls, ncol = ncol)
}

#' @describeIn MultiGenePlot Plot multiple genes for \code{SingleCellExperiment}
#' objects
#' @export MultiGenePlot.SingleCellExperiment
#' @method MultiGenePlot SingleCellExperiment
MultiGenePlot.SingleCellExperiment <- function(exp_mat, genes_plot,
                                               location_names = c("x", "y"),
                                               facet_titles = NULL,
                                               viridis_option = "D", ncol = 2,
                                               dark_theme = TRUE) {
  covariates <- genes_plot[!genes_plot %in% rownames(sce)]
  c(m, X) %<-% sce2mat(exp_mat, covariates, location_names)
  MultiGenePlot.default(m, X, genes_plot, row_sample = FALSE,facet_titles,
                        viridis_option, ncol, dark_theme)
}

#' @describeIn MultiGenePlot Plot multiple genes for \code{seurat} objects
#' @export MultiGenePlot.seurat
#' @method MultiGenePlot seurat
MultiGenePlot.seurat <- function(exp_mat, genes_plot,
                                 location_names = c("x", "y"),
                                 facet_titles = NULL,
                                 viridis_option = "D", ncol = 2,
                                 dark_theme = TRUE) {
  covariates <- genes_plot[!genes_plot %in% rownames(exp_mat@data)]
  c(m, X) %<-% seu2mat(exp_mat, covariates, location_names)
  MultiGenePlot.default(m, X, genes_plot, row_sample = FALSE,facet_titles,
                        viridis_option, ncol, dark_theme)
}
