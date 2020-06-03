
#' Quantile Normalization of Gene Count Matrix
#'
#' This function replaces original values with estimated values according
#' to the rank within each sample.
#'
#' @param gene.counts Gene count table for normalization
#' @return Normalised Counts
#' @export
quantile_norm <- function(gene.counts)
{
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
size_factor_norm <- function(gene.counts)
{
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
mrn_norm <- function(gene.counts)
{
  prenorm_mat <- t(t(gene.counts) / colSums(gene.counts))
  tpm_sf <- colMedians(prenorm_mat / prenorm_mat[, 1])
  norm_mat <- sweep(gene.counts, 2, tpm_sf / exp(mean(log(tpm_sf))), "/")
  return(norm_mat)
}
