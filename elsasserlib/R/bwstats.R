#' Run DESeq2 analysis on genome-wide bins
#'
#' Runs a DESeq2 analysis on genome-wide bins of a specified bin size.
#' The particularity of this analysis is that it skips the estimateSizeFactors
#' step by default, because this is accounted for in the scaling step of
#' MINUTE-ChIP samples.
#'
#' @param bwfiles_c1 Path or array of paths to the bigWig files for first condition.
#' @param bwfiles_c2 Path or array of paths to the bigWig files for second condition.
#' @param label_c1 Condition name for first condition
#' @param label_c2 Condition name for second condition
#' @param genome Genome. Available choices are mm9, hg38.
#' @param bin_size Bin size.
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateDispersions nbinomWaldTest `sizeFactors<-` results
#' @return a DESeqResults object as returned by DESeq2::results function
#' @export
bw_bins_diff_analysis <- function(bwfiles_c1,
                                  bwfiles_c2,
                                  label_c1,
                                  label_c2,
                                  bin_size = 10000,
                                  genome = "mm9") {

  bins_c1 <- bw_bins(bwfiles_c1, genome = genome, bin_size = bin_size)
  bins_c2 <- bw_bins(bwfiles_c2, genome = genome, bin_size = bin_size)

  cts_c1 <- get_nreads_columns(bins_c1)
  cts_c2 <- get_nreads_columns(bins_c2)

  cts <- cbind(cts_c1, cts_c2)

  condition_labels <- c(rep(label_c1, ncol(cts_c1)),
                        rep(label_c2, ncol(cts_c2)))

  coldata <- data.frame(colnames(cts), condition = condition_labels)

  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ condition)

  # Since files are scaled, we do not want to estimate size factors, so give it
  # an array of ones
  ncolumns <- length(bwfiles_c1) + length(bwfiles_c2)
  sizeFactors(dds) <- c(rep(1, ncolumns))

  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)

  results(dds)
}


#' Get values in a granges object as round numeric values in a matrix.
#'
#' This is an auxiliary function for stats. It drops NAs or NaN values, only
#' complete cases are used.
#'
#' @param granges Target granges object
#' @param fraglen Estimated fragment length
#'
#' @return An integer matrix
#' @importFrom stats complete.cases
get_nreads_columns <- function(granges, fraglen = 150) {
  # TODO: Consider whether to multiply by locus length. For bins analysis
  # this should not affect results, but for genes or loci of different length
  # it might. Since we skip the size factor step, we may bias the results?
  # So right now it's only fragment length
  length_factor <- fraglen

  bins.df <- data.frame(granges)
  cts <- as.matrix(mcols(granges))
  cts <- as.matrix(cts[complete.cases(cts),])
  cts <- round(cts*length_factor)
  cts
}


