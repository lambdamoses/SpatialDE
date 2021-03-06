% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{ensure_load}
\alias{ensure_load}
\title{Ensusre that the \code{SpatialDE} Python module is loaded}
\usage{
ensure_load()
}
\value{
An error if the Python module is not installed. Otherwise loads the
Python module.
}
\description{
By default, if the Python module \code{SpatialDE} is instatlled and \code{reticulate}
is calling the correct version of Python, then the Python module is loaded when
this R package is loaded. But the user might not have the Python module installed.
This function is to ensure that the Python module is loaded before it's called from
this R package.
}
