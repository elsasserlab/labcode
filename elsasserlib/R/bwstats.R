#' Run DESeq2 analysis on genome-wide bins
#'
#' Runs a DESeq2 analysis on genome-wide bins of a specified bin size.
#' The particularity of this analysis is that it skips the estimateSizeFactors
#' step by default, because this is accounted for in the scaling step of
#' MINUTE-ChIP samples.
#'
#' @param bwfiles_c1 Path or array of paths to the bigWig files for first condition.
#' @param bwfiles_c2 Path or array of paths to the bigWig files for second condition.
#' @param genome Genome. Available choices are mm9, hg38.
#' @param bin_size Bin size.
#' @inheritParams bw_granges_diff_analysis
#' @return a DESeqResults object as returned by DESeq2::results function
#' @export
bw_bins_diff_analysis <- function(bwfiles_c1,
                                  bwfiles_c2,
                                  label_c1,
                                  label_c2,
                                  bin_size = 10000,
                                  genome = "mm9",
                                  estimate_size_factors = FALSE,
                                  as_granges = FALSE) {

  bins_c1 <- bw_bins(bwfiles_c1, genome = genome, bin_size = bin_size)
  bins_c2 <- bw_bins(bwfiles_c2, genome = genome, bin_size = bin_size)

  bw_granges_diff_analysis(bins_c1, bins_c2, label_c1, label_c2,
                           estimate_size_factors = estimate_size_factors,
                           as_granges = as_granges)
}

#' Run DESeq2 analysis on bed file
#'
#' Runs a DESeq2 analysis on a set of loci specified in a BED file.
#' The particularity of this analysis is that it skips the estimateSizeFactors
#' step by default, because this is accounted for in the scaling step of
#' MINUTE-ChIP samples.
#'
#' @param bwfiles_c1 Path or array of paths to the bigWig files for first condition.
#' @param bwfiles_c2 Path or array of paths to the bigWig files for second condition.
#' @param bedfile BED file for locus specific analysis.
#' @inheritParams bw_granges_diff_analysis
#' @return a DESeqResults object as returned by DESeq2::results function
#' @importFrom wigglescout bw_loci
#' @export
bw_bed_diff_analysis <- function(bwfiles_c1,
                                 bwfiles_c2,
                                 bedfile,
                                 label_c1,
                                 label_c2,
                                 estimate_size_factors = FALSE,
                                 as_granges = FALSE) {

  loci_c1 <- bw_loci(bwfiles_c1, loci = bedfile)
  loci_c2 <- bw_loci(bwfiles_c2, loci = bedfile)

  bw_granges_diff_analysis(loci_c1, loci_c2, label_c1, label_c2,
                           estimate_size_factors = estimate_size_factors,
                           as_granges = as_granges)
}


#' Compute DESeq2 differential analysis on GRanges objects
#'
#' Runs a DESeq2 analysis on loci specified on GRanges objects.
#' The particularity of this analysis is that it skips the estimateSizeFactors
#' step by default, because this is accounted for in the scaling step of
#' MINUTE-ChIP samples.
#'
#' @param granges_c1 GRanges object containing the values for condition 1.
#' @param granges_c2 GRanges object containing the values for condition 2.
#'     Note that these objects MUST CORRESPOND to the same loci.
#' @param label_c1 Condition name for condition 1.
#' @param label_c2 Condition name for condition 2.
#' @param estimate_size_factors If TRUE, normal DESeq2 procedure is done. Set it
#'     to true to analyze non-MINUTE data.
#' @param as_granges If TRUE, returns a GRanges object. This is useful when
#'     loci coordinates are relevant.
#' @importFrom DESeq2 DESeqDataSetFromMatrix estimateDispersions nbinomWaldTest `sizeFactors<-` results estimateSizeFactors
#' @return a DESeqResults object as returned by DESeq2::results function. A GRanges object
#'     if as_granges parameter is TRUE.
#' @export
bw_granges_diff_analysis <- function(granges_c1,
                                     granges_c2,
                                     label_c1,
                                     label_c2,
                                     estimate_size_factors = FALSE,
                                     as_granges = FALSE) {

  # Bind first, get numbers after (drop complete cases separately could cause error)
  granges_c1 <- sortSeqlevels(granges_c1)
  granges_c1 <- sort(granges_c1)

  granges_c2 <- sortSeqlevels(granges_c2)
  granges_c2 <- sort(granges_c2)

  names_values <- NULL
  if ("name" %in% names(mcols(granges_c1))) {
    names_values <- granges_c1$name

    granges_c1$name <- NULL
  }

  if ("name" %in% names(mcols(granges_c2))) {
    granges_c2$name <- NULL
  }

  cts_df <- cbind(data.frame(granges_c1), mcols(granges_c2))

  # Needs to drop non-complete cases and match rows
  complete <- complete.cases(cts_df)
  cts_df <- cts_df[complete, ]

  cts <- get_nreads_columns(cts_df[, 6:ncol(cts_df)])

  condition_labels <- c(rep(label_c1, length(mcols(granges_c1))),
                        rep(label_c2, length(mcols(granges_c2))))

  coldata <- data.frame(colnames(cts), condition = as.factor(condition_labels))

  dds <- DESeqDataSetFromMatrix(countData = cts,
                                colData = coldata,
                                design = ~ condition,
                                rowRanges = granges_c1[complete, ])


  if (estimate_size_factors == TRUE) {
    dds <- estimateSizeFactors(dds)
  }
  else {
    # Since files are scaled, we do not want to estimate size factors, so give it
    # an array of ones
    sizeFactors(dds) <- c(rep(1, ncol(cts)))
  }

  dds <- estimateDispersions(dds)
  dds <- nbinomWaldTest(dds)

  if (as_granges) {
    result <- results(dds, format = "GRanges")
    if (!is.null(names_values)) {
      result$name <- names_values[complete]
    }

  }
  else {
   result <- results(dds, format="DataFrame")
  }

  result
}


#' Get values in a data frame object as round numeric values in a matrix.
#'
#' This is an auxiliary function for stats. It drops NAs or NaN values, only
#' complete cases are used.
#'
#' @param df Target data frame
#' @param length_factor Scaling factor to multiply coverage values by.
#'
#' @return An integer matrix
#' @importFrom stats complete.cases
get_nreads_columns <- function(df, length_factor = 1000) {
  # TODO: Consider whether to multiply by locus length. For bins analysis
  # this should not affect results, but for genes or loci of different length
  # it might. Since we skip the size factor step, we may bias the results?
  # So right now it's only fragment length
  cts <- as.matrix(df)
  cts <- as.matrix(cts[complete.cases(cts),])
  cts <- round(cts*length_factor)
  cts
}


