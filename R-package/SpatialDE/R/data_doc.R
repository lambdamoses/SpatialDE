#' Mouse Olfactory Bulb Data from Stahl et al.
#'
#' This is a sample of the mouse olfactory bulb data from
#' \href{https://doi.org/10.1126/science.aaf2403}{Stahl et al}.
#' Out of the over 16000 genes in the original data set, 1000 genes were chosen at
#' random. We used the Rep11 data. See
#' \href{http://www.spatialtranscriptomicsresearch.org/datasets/doi-10-1126science-aaf2403/}{this site}
#' for the raw data.
#'
#' @format A matrix with 1000 columns (genes) and 260 rows (spots), with the gene
#' symbols in the column name and spot names in the row name.
#' @source \url{https://doi.org/10.1126/science.aaf2403}
"mouseob"

#' Location of Spots in the Mouse Olfactory Bulb Data
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{X1}{Spots names}
#'   \item{x}{x coordinate of the location of the spots}
#'   \item{y}{y coordinate of the location of the spots}
#' }
#' @source \url{https://doi.org/10.1126/science.aaf2403}
"locations"

#' Smaller Subset of Mouse Olfactory Bulb Data
#'
#' This is a smaller subset of the mouse olfactory bulb data from
#' \href{https://doi.org/10.1126/science.aaf2403}{Stahl et al}, with only 100 randomly
#' selected genes. This dataset is only used for examples and unit testing.
#'
#' @format A matrix with 100 columns (genes) and 260 rows (spots), with the gene
#' symbols in the column name and spot names in the row name.
#' @source \url{https://doi.org/10.1126/science.aaf2403}
"mouseob_small"
