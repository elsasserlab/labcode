#' elsasserlib: Yet another utilities package.
#'
#' \code{elsasserlib} package is meant to contain all functions and utilities
#' that we use in the context of the lab, but are not large enough to make a
#' standalone package. Everybody is welcome to contribute!
#'
#' At this moment, the package provides two types of functions:
#'
#' @section Plotting themes:
#'
#' This contains some functions that play nice with \code{ggplot2}. You can use them
#' on top of any plot you make, just by using the + operator.
#'
#' There are two themes, one for screen and another one for printed version of
#' plots. This approach helps in two ways: 1) centralizes the general look of
#' the plots, so they look consistent with each other. 2) Sets up a font size
#' that can be read in large publication panels.
#'
#' @section BigWig manipulation:
#'
#' Here you can find functions to handle bigwig data files, such as binning,
#' averaging and intersecting with bed files.
#'
#' @docType package
#' @name elsasserlib
NULL
