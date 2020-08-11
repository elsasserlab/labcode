#' Build a bed-scored GRanges object from a bigWig file list and a BED file.
#'
#' Build a scored GRanges object from a bigWig file list and a BED file.
#' The aggregating function (per locus) can be min, max, sd, mean.
#'
#' bwfiles and bg_bwfiles must have the same length. If you are using same
#' background for several files, then file paths must be repeated accordingly.
#'
#' Values can be normalized using background bigWig files (usually input
#' files). By default, the value obtained will be bigwig / bg_bigwig per bin,
#' per bigWig.
#'
#' If norm_func is specified, this can be changed to any given function, for
#' instance, if norm_func = log2, values will represent log2(bigwig / bg_bigwig)
#' per bin.
#'
#' @param bedfile BED file to intersect with the BigWig file.
#' @param aggregate_by Statistic to aggregate per group. If NULL, values are
#'    not aggregated. This is the behavior by default.
#' @export
#' @inheritParams bw_bins
#' @importFrom rtracklayer import BigWigFile
#' @importFrom GenomeInfoDb sortSeqlevels
bw_bed <- function(bwfiles,
                   bedfile,
                   bg_bwfiles=NULL,
                   labels=NULL,
                   per_locus_stat="mean",
                   aggregate_by=NULL,
                   norm_func=identity) {

  validate_filelist(bwfiles)
  validate_filelist(bedfile)

  if (is.null(labels)) {
    labels <- basename(bwfiles)
  }

  bed <- import(bedfile)
  bed <- sortSeqlevels(bed)
  bed <- sort(bed, ignore.strand=FALSE)

  result <- NULL
  if (is.null(bg_bwfiles)) {
    result <- multi_bw_ranges(bwfiles,
                              labels,
                              granges=bed,
                              per_locus_stat=per_locus_stat)

  } else {
    if (is.null(aggregate_by)) {
      # Only want to normalize per-locus if not aggregating
      result <- multi_bw_ranges_norm(bwfiles,
                                     bg_bwfilelist=bg_bwfiles,
                                     labels=labels,
                                     granges=bed,
                                     per_locus_stat=per_locus_stat,
                                     norm_func=norm_func)
    }
  }

  if ( 'name' %in% names(mcols(bed)) ) {
    result$name <- bed$name
  }

  if (!is.null(aggregate_by)) {
    df <- aggregate_scores(result,
                           group_col="name",
                           aggregate_by=aggregate_by)

    result <- natural_sort_by_field(df, "name")

    if (!is.null(bg_bwfiles)) {
      bg <- multi_bw_ranges(bg_bwfiles,
                            labels,
                            granges=bed,
                            per_locus_stat=per_locus_stat)

      if ( 'name' %in% names(mcols(bed)) ) {
        bg$name <- bed$name
      }

      bg_df <- aggregate_scores(bg,
                                group_col="name",
                                aggregate_by=aggregate_by)

      values <- cbind(norm_func(df[, labels]/ bg_df[, labels]), df$name)
      colnames(values) <- c(labels, "name")
      result <- natural_sort_by_field(values, "name")
    }
  }

  result
}


#' Build a binned-scored GRanges object from a bigWig file
#'
#' Build a binned-scored GRanges object from a bigWig file. The aggregating
#' function per bin can be min, max, sd, mean.
#'
#' bwfiles and bg_bwfiles must have the same length. If you are using same
#' background for several files, then file paths must be repeated accordingly.
#'
#' Values can be normalized using background bigWig files (usually input
#' files). By default, the value obtained will be bigwig / bg_bigwig per bin,
#' per bigWig.
#'
#' If norm_func is specified, this can be changed to any given function, for
#' instance, if norm_func = log2, values will represent log2(bigwig / bg_bigwig)
#' per bin.
#'
#' @param bwfiles Path or array of paths to the bigWig files to be summarized.
#' @param bg_bwfiles Path or array of paths to the bigWig files to be used as
#'   background.
#' @param labels List of names to give to the mcols of the returned GRanges
#'     object. If NULL, file names are used.
#' @param per_locus_stat Aggregating function (per locus). Mean by default.
#'     Choices: min, max, sd, mean.
#' @param bin_size Bin size.
#' @param genome Genome. Available choices are mm9, hg38.
#' @param selection A GRanges object to restrict binning to a certain set of
#'     intervals. This is useful for debugging and improving performance of locus
#'     specific analyses.
#' @param norm_func Function to apply to normalize bin values f(bw / bw_bg).
#' @return A GRanges object with each bwfile as a metadata column named
#'     after labels, if provided, or after filenames otherwise.
#' @export
bw_bins <- function(bwfiles,
                    bg_bwfiles=NULL,
                    labels=NULL,
                    per_locus_stat="mean",
                    bin_size=10000,
                    genome="mm9",
                    selection=NULL,
                    norm_func=identity) {

  validate_filelist(bwfiles)

  if (is.null(labels)) {
    labels <- make.names(basename(bwfiles))
  }

  tiles <- build_bins(bin_size=bin_size, genome=genome)

  if (is.null(bg_bwfiles)) {
    result <- multi_bw_ranges(bwfiles,
                              labels,
                              tiles,
                              per_locus_stat=per_locus_stat,
                              selection=selection)
  } else {
    # FIXME: mcols of result may end up being <matrix> instead of <numeric>.
    result <- multi_bw_ranges_norm(bwfiles,
                                   bg_bwfiles,
                                   labels,
                                   tiles,
                                   per_locus_stat=per_locus_stat,
                                   selection=selection,
                                   norm_func=norm_func)
  }

  result
}


#' Calculate profile values of a bigWig file over a BED file.
#'
#' Calculates profile values of a set of tracks over the loci speficied in a
#' BED file. Data points are taken each bin_size base pairs.
#'
#' Loci are aligned depending on mode parameter:
#'
#' - stretch. Aligns all starts and all ends, sort of stretching the loci.
#' The median of these lenghts is taken as the pseudo-length in order to show
#' a realistic plot when displayed.
#'
#' - start. All loci are aligned by start.
#'
#' - end. All loci are aligned by end.
#'
#' - center. All loci are aligned by center.
#'
#' @param bedfile BED file to summarize
#' @param mode How to handle differences in lengths across loci:
#'
#'   stretch: Anchor each locus on both sides.
#'
#'   start: Anchor all loci on start.
#'
#'   end: Anchor all loci on end.
#'
#'   center: Center all loci.
#'
#' @param bin_size Bin size. Length of bin in base pairs. The lower, the higher the resolution.
#' @param upstream Number of base pairs to include upstream of loci.
#' @param downstream Number of base pairs to include downstream of loci.
#' @param ignore_strand Whether to use strand information in BED file.
#' @inheritParams bw_bins
#' @return a data frame in long format
#' @export
bw_profile <- function(bwfiles,
                       bg_bwfiles=NULL,
                       bedfile=NULL,
                       labels=NULL,
                       mode="stretch",
                       bin_size=100,
                       upstream=2500,
                       downstream=2500,
                       ignore_strand=F,
                       norm_func=identity) {

  validate_filelist(bwfiles)
  validate_filelist(bedfile)
  granges <- rtracklayer::import(bedfile)

  if (bin_size <= 0) {
    stop(paste("bin size must be a positive value:", bin_size))
  }

  if (upstream <= 0) {
    stop(paste("upstream size must be a positive value:", upstream))
  }

  if (downstream <= 0) {
    stop(paste("downstream size must be a positive value:", downstream))
  }

  if (bin_size > upstream || bin_size > downstream) {
    stop("bin size must be smaller than flanking regions")
  }

  if (is.null(labels)) {
    labels <- basename(bwfiles)
  }

  if (length(bwfiles) != length(labels)) {
    stop("labels and bwfiles must have the same length")
  }

  calculate_bw_profile_fixed <- purrr::partial(calculate_bw_profile,
                                               granges=granges,
                                               mode=mode,
                                               bin_size=bin_size,
                                               upstream=upstream,
                                               downstream=downstream,
                                               ignore_strand=ignore_strand,
                                               norm_func=norm_func)
  if (is.null(bg_bwfiles)) {
    values_list <- purrr::map2(bwfiles,
                               labels,
                               calculate_bw_profile_fixed,
                               bg_bw=NULL)
  } else {
    values_list <- purrr::pmap(list(bwfiles, bg_bwfiles, labels),
                               calculate_bw_profile_fixed)
  }

  values <- do.call(rbind, values_list)
  values
}


#' Build a unscored bins GRanges object.
#'
#' Build a GRanges of bins of a given size, for a specific genome. Supported
#' genomes (required for the package): mm9, hg38.
#'
#' @param bin_size Bin size.
#' @param genome Genome.
#' @importFrom GenomicRanges tileGenome
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom BSgenome.Mmusculus.UCSC.mm9 BSgenome.Mmusculus.UCSC.mm9
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @return A GRanges object
#' @export
build_bins <- function(bin_size=10000, genome='mm9') {
  seq_lengths <- NULL

  if (genome == "mm9") {
    seq_lengths <- seqinfo(BSgenome.Mmusculus.UCSC.mm9)
  } else {
    if (genome == "hg38") {
      seq_lengths <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
    } else {
      stop("Supported genomes: mm9, hg38")
    }
  }

  bins <- tileGenome(seq_lengths,
                     tilewidth=bin_size,
                     cut.last.tile.in.chrom=T)

  bins
}


#' Intersect a list of bigWig files with a GRanges object
#'
#' Build a binned-scored GRanges object from a list of bigWig files. The
#' aggregating function per locus can be min, max, sd, mean.
#'
#' @param granges GRanges object to intersect
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomeInfoDb sortSeqlevels
#' @inheritParams bw_bins
#' @return A sorted GRanges object.
multi_bw_ranges <- function(bwfiles,
                            labels,
                            granges,
                            per_locus_stat="mean",
                            selection=NULL) {

  if (length(bwfiles) != length(labels)) {
    stop("BigWig file list and column names must have the same length.")
  }

  summaries <- purrr::map(bwfiles,
                          bw_ranges,
                          granges=granges,
                          per_locus_stat=per_locus_stat,
                          selection=selection)

  # granges_cbind sorts each element so it's safer to merge and no need to
  # sort after
  result <- granges_cbind(summaries, labels)
  result
}


#' Intersect a list of bigWig files with a GRanges object (with background)
#'
#' Intersect a list of bigWig files with a GRanges object and normalize values
#' with a set of background bigWig files.
#'
#' @param bwfilelist BigWig file list to be summarized.
#' @param bg_bwfilelist Background BigWig files.
#' @param granges GRanges object to intersect
#' @param selection A GRanges object to restrict analysis to.
#' @param norm_func Function to apply after bw/bg_bw.
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @inheritParams bw_bins
#' @return a sorted GRanges object
multi_bw_ranges_norm <- function(bwfilelist,
                                 bg_bwfilelist,
                                 labels,
                                 granges,
                                 per_locus_stat="mean",
                                 selection=NULL,
                                 norm_func=identity) {

  if (length(bwfilelist) != length(bg_bwfilelist)) {
    stop("Background and signal bwfile lists must have the same length.")
  }

  result <- multi_bw_ranges(bwfilelist,
                            labels,
                            granges,
                            per_locus_stat=per_locus_stat,
                            selection=selection)

  bg <- multi_bw_ranges(bg_bwfilelist,
                        labels,
                        granges,
                        per_locus_stat=per_locus_stat,
                        selection=selection)

  result_df <- data.frame(result)
  result_df[, labels] <- norm_func(as.matrix(mcols(result)) / as.matrix(mcols(bg)))

  makeGRangesFromDataFrame(result_df, keep.extra.columns = T)
}

#' Score a GRanges object against a BigWig file
#'
#' Build a scored GRanges object from a GRanges object and a single BigWig file.
#' The aggregating function can be min, max, sd, mean.
#'
#' @param bwfile Path to a single BigWig file to be summarized.
#' @param granges GRanges object file to be summarized.
#' @importFrom rtracklayer BigWigFile
#' @importFrom IRanges subsetByOverlaps
#' @importFrom methods getMethod
#' @inheritParams bw_bins
#' @return GRanges with column score.
bw_ranges <- function (bwfile, granges, per_locus_stat='mean', selection=NULL) {
  bw <- BigWigFile(bwfile)
  explicit_summary <- getMethod("summary", "BigWigFile")

  if (! is.null(selection)) {
    granges <- subsetByOverlaps(granges, selection)
  }
  result <- unlist(explicit_summary(bw, granges, type=per_locus_stat))
  result
}


#' Aggregate scores of a GRanges object on a field
#'
#' Aggregates scores of a GRanges object on a specific field. Used for summary
#' functions.
#'
#' @param scored_granges A GRanges object with numerical metadata columns
#' @param group_col A column among the mcols that can be seen as a factor.
#' @param aggregate_by Function used to aggregate: mean, median, true_mean.
#'
#'     true_mean: Mean coverage taking all elements in a class as one large bin.
#'
#'     mean: mean-of-distribution approach. The mean of the aggregated value per
#'       locus is reported.
#'
#'     median: median-of-distribution. The median of the aggregated value per
#'       locus is reported.
#'
#' @importFrom dplyr group_by_at summarise across `%>%`
#' @importFrom rtracklayer mcols
#' @return A data frame with aggregated scores.
aggregate_scores <- function(scored_granges, group_col, aggregate_by) {
  if ( !is.null(group_col) && group_col %in% names(mcols(scored_granges))) {

    # GRanges objects are 1-based and inclusive [start, end]
    scored_granges$length <- GenomicRanges::end(scored_granges) - GenomicRanges::start(scored_granges) + 1

    df <- data.frame(mcols(scored_granges))
    validate_categories(df[, group_col])
    # Make sure special characters are taken into account
    score_cols <- colnames(mcols(scored_granges))
    score_cols <- make.names(score_cols)
    score_cols <- score_cols[!score_cols %in% c(group_col)]

    if (aggregate_by == "true_mean") {
      sum_vals <- df[, score_cols]*df$length
      colnames(sum_vals) <- score_cols
      sum_vals[, group_col] <- df[, group_col]
      sum_vals$length <- df$length

      # Summarize SUM only
      sum_vals <- sum_vals %>%
        group_by_at(group_col) %>%
        summarise(across(where(is.numeric), sum))

      # Divide sum(scores) by sum(length) and keep only scores
      df <- sum_vals[, score_cols]/sum_vals$length
      df[, group_col] <- sum_vals[, group_col]

    } else if (aggregate_by %in% c("mean", "median")) {
      f <- get(aggregate_by)
      df <- df %>%
        group_by_at(group_col) %>%
        summarise(across(where(is.numeric), f))

    } else {
      stop(paste("Function not implemented as aggregate_by:", aggregate_by))
    }

    score_cols <- score_cols[! score_cols %in% c('length')]
    df <- df[, c(score_cols, group_col), drop=FALSE]
    data.frame(df)

  } else {
    stop("Grouping column not provided or not present in GRanges object.")
  }
}



#' Calculate profile values of a bigWig file over GRanges object.
#'
#' Calculates profile values of a bigWig file over the loci speficied in a
#' GRanges object. Data points are taken each bin_size base pairs.
#'
#' Loci are aligned depending on mode parameter:
#'
#' - stretch. Aligns all starts and all ends, sort of stretching the loci.
#' The median of these lenghts is taken as the pseudo-length in order to show
#' a realistic plot when displayed.
#'
#' - start. All loci are aligned by start.
#'
#' - end. All loci are aligned by end.
#'
#' - center. All loci are aligned by center.
#'
#' Adapted from seqplots rtracklayer wrapping functions to compute coverage values:
#' For more on seqplots: https://www.bioconductor.org/packages/release/bioc/html/seqplots.html
#'
#' @param bw BigWig file to be summarized.
#' @param granges GRanges object
#' @param bg_bw BigWig file to be used as background.
#' @param label Name to give to the values
#' @importFrom rtracklayer BigWigFile import
#' @inheritParams bw_profile
#' @return A DataFrame with the aggregated scores
calculate_bw_profile <- function(bw,
                                 granges,
                                 bg_bw=NULL,
                                 label=NULL,
                                 mode="stretch",
                                 bin_size=100,
                                 upstream=2500,
                                 downstream=2500,
                                 ignore_strand=F,
                                 norm_func=identity) {


  bwfile <- BigWigFile(path=bw)

  if (is.null(label)) {
    label <- basename(bw)
  }

  if (mode == 'stretch') {

    full <- calculate_stretch_matrix(bwfile,
                                     granges,
                                     bin_size=bin_size,
                                     upstream=upstream,
                                     downstream=downstream,
                                     ignore_strand=ignore_strand)

    if (!is.null(bg_bw)) {
      bg_bwfile <- BigWigFile(path=bg_bw)
      bg <- calculate_stretch_matrix(bg_bwfile,
                                   granges,
                                   bin_size=bin_size,
                                   upstream=upstream,
                                   downstream=downstream,
                                   ignore_strand=ignore_strand)

      full <- norm_func(full/bg)
    }

  } else {
    granges <- GenomicRanges::promoters(GenomicRanges::resize(granges, 1, fix=mode),
                                   upstream,
                                   downstream)

    npoints <- floor((upstream + downstream) / bin_size)

    full <- intersect_bw_and_granges(bwfile, granges, npoints=npoints, ignore_strand=F)

    if (!is.null(bg_bw)) {
      bg_bwfile <- BigWigFile(path=bg_bw)
      bg <- intersect_bw_and_granges(bg_bwfile, granges, npoints=npoints, ignore_strand=F)

      full <- norm_func(full/bg)
    }
  }

  result_df <- summarize_matrix(full, label)
}

#' Calculate matrix for stretch mode
#'
#' @param bw BigWigFile object
#' @param granges GRanges object
#' @param bin_size Bin size
#' @param upstream Number of basepairs upstream
#' @param downstream Number of basepairs downstream
#' @param ignore_strand Ignore strand (bool)
#'
#' @return Summary matrix
calculate_stretch_matrix <- function(bw,
                                   granges,
                                   bin_size=100,
                                   upstream=2500,
                                   downstream=2500,
                                   ignore_strand=F) {

  left_npoints <- floor(upstream/bin_size)
  right_npoints <- floor(downstream/bin_size)
  # Stretch to the median value of the GR object
  middle_npoints <- floor(median(GenomicRanges::width(granges))/bin_size )

  left <- intersect_bw_and_granges(bw,
                         GenomicRanges::flank(granges, upstream, start=TRUE),
                         npoints=left_npoints,
                         ignore_strand=ignore_strand)

  right <- intersect_bw_and_granges(bw,
                          GenomicRanges::flank(granges, downstream, start=FALSE),
                          npoints=right_npoints,
                          ignore_strand=ignore_strand)

  middle <- intersect_bw_and_granges(bw,
                           granges,
                           npoints=middle_npoints,
                           ignore_strand=ignore_strand)

  cbind(left, middle, right)

}



#' Intersect a BigWig file over loci on a GRanges object
#'
#' Intersect a BigWig file over loci on a GRanges object, taking npoints points
#' per locus.
#'
#' @param bw BigWigFile object (rtracklayer).
#' @param granges GRanges object.
#' @param npoints How many points to take (different to bin size!).
#' @param ignore_strand Ignore strand information in granges.
#' @return A value matrix of dimensions len(granges) x npoints.
intersect_bw_and_granges <- function(bw, granges, npoints, ignore_strand=F) {
  values <- rtracklayer::summary(bw,
                    which=granges,
                    as='matrix',
                    size=npoints)

  # Reverse minus strand rows
  if (!ignore_strand) {
    values[as.character(GenomicRanges::strand(granges))=="-",] <- values[
      as.character(GenomicRanges::strand(granges))=="-", ncol(values):1]
  }

  values
}


#' Summarize a intersect_bw_and_granges matrix
#'
#' Compute averages and standard error values for a matrix returned by
#' intersect_bw_and_granges.
#'
#' @param matrix A matrix returned by intersect_bw_and_granges.
#' @param label Label for the sample.
#' @importFrom stats median sd
#' @return A dataframe with summarized values, sderror and medians, plus a label.
summarize_matrix <- function(matrix, label) {
  # Ignore Inf and NaN in the computation of means/SD
  matrix[is.infinite(matrix)] <- NA

  omitted_vals <- sum(is.na(matrix))
  if(omitted_vals > 100) {
    mean_per_locus <- omitted_vals / nrow(matrix)
    warning(paste("Profile plot:", omitted_vals, "generated (", mean_per_locus,"per locus)"))
  }


  df <- data.frame(mean=colMeans(matrix, na.rm=TRUE),
                   sderror=apply(matrix, 2,
                                 function (n) { sd(n, na.rm=TRUE) / sqrt( sum(!is.na(n))) } ),
                   median=apply(matrix, 2, median, na.rm=TRUE))

  df$index <- as.integer(rownames(df))
  df$sample <- label
  df
}


#' GRanges cbind-like operation
#'
#' Perform a cbind operation on a GRanges list, appending scores to mcols.
#' It will sort the GRanges elements in order to ensure the match is proper.
#'
#' This assumes that you are trying to cbind things that match (bins generated
#' from the same parameters, BED intersections from the same BED.)
#'
#' @param grlist A list of GRanges objects that have all the same fields.
#' @param labels Vector of names for the score columns.
#' @importFrom GenomeInfoDb sortSeqlevels
granges_cbind <- function(grlist, labels) {
  fixed_fields <- c("seqnames", "start", "end", "width", "strand")

  grlist[[1]] <- sortSeqlevels(grlist[[1]])
  grlist[[1]] <- sort(grlist[[1]])

  result <- data.frame(grlist[[1]])[, fixed_fields]
  for (i in seq(1, length(grlist))) {
    grlist[[i]] <- sortSeqlevels(grlist[[i]])
    grlist[[i]] <- sort(grlist[[i]])

    result[, labels[[i]]] <- grlist[[i]]$score
  }

  result <- makeGRangesFromDataFrame(result, keep.extra.columns=T)
  result
}


#' Take a column of a data frame to use as row names, drop such column and
#' reorder the data frame so row names follow natural order.
#'
#' @param df A dataframe
#' @param col The column (usually name)
#' @importFrom stringr str_sort
#' @return A sorted df
natural_sort_by_field <- function(df, col) {
  rownames(df) <- df[, col]
  order <- str_sort(df[,col], numeric=TRUE)
  df[, col] <- NULL
  df[order, , drop=FALSE]
}


#' Validate a category array
#'
#' Checks whether the number of categories is reasonable for an aggregation.
#' Throws a warning if it finds more than 50 different values.
#'
#' @param cat_values An array of values
validate_categories <- function(cat_values) {
  MAX_CATEGORIES <- 50
  # Test number of values in group_col
  ncat <- length(levels(as.factor(cat_values)))
  if (ncat > MAX_CATEGORIES) {
    warning(paste("Number of values in group column field very large:",
                  ncat,
                  "(does BED file have unique IDs instead of categories?)"))
  }
}


#' Validate an array of paths
#'
#' Check that a list of files is valid: not empty and contents exist.
#'
#' @param filelist An array of files
#' @return NULL
validate_filelist <- function(filelist) {
  if (length(filelist) == 0) {
    stop("File list provided is empty.")
  }

  if (!all(file.exists(filelist))) {
    msg <- paste("Files not found:", filelist[!file.exists(filelist)])
    stop(msg)
  }
}

