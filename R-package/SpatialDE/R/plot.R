#' Inverse log scale
#'
#' Custom transform for inverted log transformed axis in ggplot2.
#' @importFrom scales trans_new
reverse_log10 <- function() {
  trans <- function(x) -log10(x)
  inv <- function(x) 10^(-x)
  scales::trans_new("reverse_log10", trans, inv, log_breaks(base = 10), domain = c(1e-100, Inf))
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
#' @param do_label Display gene names for statistically significant genes, default to \code{TRUE}.
#' @param covariate_names Names of covariates as a reference.
#'
#' @return A \code{ggplot2} object.
#' @import ggplot2
#' @importFrom dplyr mutate filter select left_join case_when
#' @importFrom magrittr %>%
#' @export
FSV_sig <- function(results, ms_results = NULL, certain_only = FALSE, do_label = TRUE,
                    covariate_names = "log_total_count") {
  if (!is.null(ms_results)) {
    results <- results %>%
      select(-model) %>%
      left_join(ms_results, by = "g")
  }
  results <- results %>%
    mutate(FSV95conf = 2 * sqrt(s2_FSV),
           conf_categories = fct_rev(cut(FSV95conf, c(0, 0.1, 1, Inf))),
           is_covariate = g %in% covariate_names,
           # More user friendly model labels
           color_categories = case_when(model == "SE" ~ "general",
                                        model == "PER" ~ "periodic"),
           color_categories = ifelse(qval < 0.05 & !is_covariate, color_categories, NA))
  if (certain_only) {
    results <- results %>%
      filter(conf_categories == "(0,0.1]")
  }
  colors_use <- scales::hue_pal()(length(unique(results$model)))
  p <- ggplot(results, aes(FSV, qval)) +

    geom_hline(yintercept = 0.05, linetype = 2) +
    scale_color_manual(values = colors_use, na.translate = TRUE, na.value = "black",
                       guide = guide_legend(title = "model")) +
    scale_y_continuous(trans = reverse_log10()) +
    scale_shape_manual(values = c(16,4),
                       guide = guide_legend(title = "covariate")) +
    annotate(geom = "text", x = 0, y = 0.05, label = "0.05")
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
  if (do_label) {
    gene_label <- results %>%
      filter(qval < 0.05) %>%
      select(FSV, qval, g)
    p <- p + ggrepel::geom_label_repel(aes(label = g), data = gene_label)
  }
  p
}
