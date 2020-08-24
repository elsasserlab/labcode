
#' Quantile Normalization of Gene Count Matrix
#'
#' This function replaces original values with estimated values according
#' to the rank within each sample.
#'
#' @param gene.counts Gene count table for normalization
#' @return Normalised Counts
#' @export
quantile_norm <- function(gene.counts) {
  rank_mat <- apply(gene.counts, 2, rank, ties.method = "random")
  means <- rowMeans(gene.counts)
  rank_mean <- rank(means, ties.method = "random")
  norm_mat <- apply(rank_mat, 2, function(x) means[x])
  return(norm_mat)
}

#' Size Factor Normalization of Gene Count Matrix
#'
#' This function devides original values with size factors of each sample according
#' to DESeq method.
#'
#' @param gene.counts Gene count table for normalization
#' @return Normalised Counts
#' @export
size_factor_norm <- function(gene.counts) {
  geoMeans = apply(gene.counts, 1, function(x) exp(sum(log(x[x > 0]), na.rm=T) / length(x)) )
  sf = apply(sweep(gene.counts, 1, geoMeans,'/'), 2, median)
  return(sweep(gene.counts, 2, sf,'/'))
}

#' Median Ratio Normalization of Gene Count Matrix
#'
#' This function devides original values with median ratios comparing to the first sample
#' according to edgeR method.
#'
#' @param gene.counts Gene count table for normalization
#' @return Normalised Counts
#' @export
mean_ratio_norm <- function(gene.counts) {
  prenorm_mat <- t(t(gene.counts) / colSums(gene.counts))
  tpm_sf <- colMedians(prenorm_mat / prenorm_mat[, 1])
  norm_mat <- sweep(gene.counts, 2, tpm_sf / exp(mean(log(tpm_sf))), "/")
  return(norm_mat)
}

#' Limit Outliers of Gene Count Matrix
#'
#' This function trims extreme values above or below the indicated quantile, 
#' processes by column and returns as a matrix.
#'
#' @param x Gene count table or a vector
#' @param q Limit of quantile
#' @return Normalised Counts
#' @export
trim_quantile <- function(x, q = 0.995) {
  if (!is.null(dim(x))) {
    cbind(trim_quantile(x[, 1], q), trim_quantile(x[, -1], q))
  } else {
    x[x > quantile(x, q, na.rm = T)] = quantile(x, q, na.rm = T)
    x[x < quantile(x, 1-q, na.rm = T)] = quantile(x, 1-q, na.rm = T)
    return(x)
  }
}