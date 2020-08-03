
#' Build a binned-scored GRanges object from a bigWig file
#'
#' Build a binned-scored GRanges object from a bigWig file. The aggregating
#' function can be min, max, sd, mean.
#'
#' @param bwfiles BigWig file to be summarized (or list).
#' @param colnames List of names to give to the mcols of the returned GRanges
#'     object. If NULL, filenames are used (default).
#' @param stat Aggregating function (per locus). Mean by default.
#'     Choices: min, max, sd, mean.
#' @param bsize Bin size. Default 10000.
#' @param genome Genome. Available choices are mm9, hg38.
#' @param selection A GRanges object to restrict binning to a certain set of
#'     intervals. It is useful for debugging and improving performance of locus
#'     specific analyses.
#' @return A GenomicRanges object with each bwfile as a metadata column named
#'     after colnames.
#' @export
bw_bins <- function(bwfiles,
                    colnames=NULL,
                    stat='mean',
                    bsize=10000,
                    genome='mm9',
                    selection=NULL) {

  check_filelist(bwfiles)

  if (is.null(colnames)) {
    colnames <- basename(bwfiles)
  }

  tiles <- build_bins(bsize=bsize, genome=genome)
  result <- multi_bw_ranges(bwfiles,
                            colnames,
                            tiles,
                            per.locus.stat=stat,
                            selection=selection)
  result
}

#' Intersect a list of bw files with a GRanges object
#'
#' Build a binned-scored GRanges object from a list of bigWig files. The
#' aggregating function per locus can be min, max, sd, mean.
#'
#' @param bwfilelist BigWig file to be summarized.
#' @param colnames Names to be assigned to the columns
#' @param gr GRanges object to intersect
#' @param per.locus.stat Aggregating function per stat
#' @param selection A GRanges object to restrict analysis to.
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomeInfoDb sortSeqlevels
#' @return a sorted GRanges object
multi_bw_ranges <- function(bwfilelist,
                            colnames,
                            gr,
                            per.locus.stat='mean',
                            selection=NULL) {

  if (length(bwfilelist) != length(colnames)) {
    stop("BigWig file list and column names must have the same length.")
  }

  summaries <- purrr::map(bwfilelist,
                          bw_ranges,
                          gr=gr,
                          per.locus.stat=per.locus.stat,
                          selection=selection)

  # granges_cbind sorts each element so it's safer to merge and no need to
  # sort after
  result <- granges_cbind(summaries, colnames)
  result
}

#' Performs a cbind operation on a GRanges list, appending scores
#'
#' It will sort the GRanges elements in order to ensure the match is proper.
#'
#' @param grlist A list of GRanges objects that have all the same fields.
#' @param colnames Vector of names for the score columns.
#' @importFrom GenomeInfoDb sortSeqlevels
granges_cbind <- function(grlist, colnames) {
  fixed.fields <- c('seqnames', 'start', 'end', 'width', 'strand')

  grlist[[1]] <- sortSeqlevels(grlist[[1]])
  grlist[[1]] <- sort(grlist[[1]])

  result <- data.frame(grlist[[1]])[, fixed.fields]
  for (i in seq(1, length(grlist))) {
    grlist[[i]] <- sortSeqlevels(grlist[[i]])
    grlist[[i]] <- sort(grlist[[i]])

    result[, colnames[[i]]] <- grlist[[i]]$score
  }

  result <- makeGRangesFromDataFrame(result, keep.extra.columns=T)
  result
}

#' Build a scored GRanges object from a BED file.
#'
#' Build a scored GRanges object from a bigWig file and a BED file.
#' The aggregating function (per locus) can be min, max, sd, mean.
#'
#' @param bwfiles BigWig file (or list) to be summarized.
#' @param bedfile BED file to intersect with the BigWig file.
#' @param colnames Column names of the score fields. Must have the same
#'    length as bigwig file list. If not provided, colnames are the names of
#'    the files in bwfiles.
#' @param per.locus.stat Aggregate per locus function.
#' @param aggregate.by Statistic to aggregate per group. If NULL, values are
#'    not aggregated. This is the behavior by default.
#' @export
#' @importFrom rtracklayer import BigWigFile
#' @importFrom GenomeInfoDb sortSeqlevels
bw_bed <- function(bwfiles,
                   bedfile,
                   colnames=NULL,
                   per.locus.stat='mean',
                   aggregate.by=NULL) {

  check_filelist(bwfiles)
  check_filelist(bedfile)

  if (is.null(colnames)) {
    colnames <- basename(bwfiles)
  }

  bed <- import(bedfile)
  bed <- sortSeqlevels(bed)
  bed <- sort(bed, ignore.strand=FALSE)

  result <- multi_bw_ranges(bwfiles,
                            colnames,
                            gr=bed,
                            per.locus.stat=per.locus.stat)

  if ( 'name' %in% names(mcols(bed)) ) {
    result$name <- bed$name
  }
  if (! is.null(aggregate.by)) {
    df <- aggregate_scores(result,
                           group.col='name',
                           aggregate.by=aggregate.by)

    result <- natural_sort_by_field(df, 'name')
  }
  result
}

check_filelist <- function(filelist) {
  if (length(filelist) == 0) {
    stop("File list provided is empty.")
  }

  if (!all(file.exists(filelist))) {
    msg <- paste("Files not found:", filelist[!file.exists(filelist)])
    stop(msg)
  }
}

#' Moves a column of a dataframe to rownames and naturally sorts the rows
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

#' Build a bins GRanges object.
#'
#' Build a GRanges of bins of a given size, for a specific genome. Supported
#' genomes (required for the package): mm9, hg38.
#'
#' @param bsize Bin size. 10000 by default.
#' @param genome Genome used. mm9 by default. Valid values: mm9, hg38
#' @return A GRanges object
#' @export
#' @importFrom GenomicRanges tileGenome
#' @importFrom GenomeInfoDb seqinfo
#' @importFrom BSgenome.Mmusculus.UCSC.mm9 BSgenome.Mmusculus.UCSC.mm9
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
build_bins <- function(bsize=10000, genome='mm9') {
  seq_lengths <- NULL

  if (genome == 'mm9') {
    seq_lengths <- seqinfo(BSgenome.Mmusculus.UCSC.mm9)
  } else {
    if (genome == 'hg38') {
      seq_lengths <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
    } else {
      stop("Supported genomes: mm9, hg38")
    }
  }

  bins <- tileGenome(seq_lengths,
                     tilewidth=bsize,
                     cut.last.tile.in.chrom=T)

  bins
}

#' Score a GRanges object against a BigWig file
#'
#' Build a scored GRanges object from a GRanges object and a BigWig file.
#' The aggregating function can be min, max, sd, mean.
#'
#' @param bwfile BigWig file to be summarized.
#' @param gr GRanges object file to be summarized.
#' @param per.locus.stat Aggregating function (per locus). Mean by default.
#'    Choices: min, max, sd, mean. These choices depend on rtracklayer library.
#' @param selection A GRanges object to restrict analysis to.
#' @importFrom rtracklayer BigWigFile
#' @importFrom IRanges subsetByOverlaps
#' @importFrom methods getMethod
#' @return GRanges with column score.`
bw_ranges <- function (bwfile, gr, per.locus.stat='mean', selection=NULL) {
  bw <- BigWigFile(bwfile)
  explicit_summary <- getMethod("summary", "BigWigFile")

  if (! is.null(selection)) {
    gr <- subsetByOverlaps(gr, selection)
  }
  result <- unlist(explicit_summary(bw, gr, type=per.locus.stat))
  result
}


#' Checks that the number of categories is reasonable for an aggregation.
#'
#' Throws a warning if found more than  50 values.
#'
#' @param cat.values An array of values
validate_categories <- function(cat.values) {
  MAX_CATEGORIES <- 50
  # Test number of values in group.col
  ncat <- length(levels(as.factor(cat.values)))
  if (ncat > MAX_CATEGORIES) {
    warning(paste("Number of values in group column field very large:",
                  ncat,
                  "(does BED file have unique IDs instead of categories?)"))
  }
}

#'
#' Aggregate scores of a GRanges object on a specific field
#'
#' Aggregates scores of a GRanges object on a specific field.
#' @param scored.gr A GRanges object with numerical metadata columns
#' @param group.col A column among the mcols that can be seen as a factor.
#' @param aggregate.by Function used to aggregate: mean, median, true_mean.
#'     true_mean: Mean coverage taking all elements in a class as one large bin
#'     mean: mean-of-distribution approach. The mean of the aggregated value per
#'       locus is reported.
#'     median: median-of-distribution. The median of the aggregated value per
#'       locus is reported
#' @return A DataFrame with the aggregated scores (any numerical column will be
#'     aggregated).
#' @importFrom dplyr group_by_at summarise across `%>%`
#' @importFrom rtracklayer mcols
aggregate_scores <- function(scored.gr, group.col, aggregate.by) {
  if ( !is.null(group.col) && group.col %in% names(mcols(scored.gr))) {

    # GRanges objects are 1-based and inclusive [start, end]
    scored.gr$length <- GenomicRanges::end(scored.gr) - GenomicRanges::start(scored.gr) + 1

    df <- data.frame(mcols(scored.gr))
    validate_categories(df[, group.col])
    # Make sure special characters are taken into account
    score.cols <- colnames(mcols(scored.gr))
    score.cols <- make.names(score.cols)
    score.cols <- score.cols[!score.cols %in% c(group.col)]

    if (aggregate.by == 'true_mean') {
       sum.vals <- df[, score.cols]*df$length
       colnames(sum.vals) <- score.cols
       sum.vals[, group.col] <- df[, group.col]
       sum.vals$length <- df$length

       # Summarize SUM only
       sum.vals <- sum.vals %>%
         group_by_at(group.col) %>%
         summarise(across(where(is.numeric), sum))

       # Divide sum(scores) by sum(length) and keep only scores
       df <- sum.vals[, score.cols]/sum.vals$length
       df[, group.col] <- sum.vals[, group.col]

    } else if (aggregate.by %in% c('mean', 'median')) {
      f <- get(aggregate.by)
      df <- df %>%
        group_by_at(group.col) %>%
        summarise(across(where(is.numeric), f))

    } else {
      stop(paste("Function not implemented as aggregate.by:", aggregate.by))
    }

    score.cols <- score.cols[! score.cols %in% c('length')]
    df <- df[, c(score.cols, group.col), drop=FALSE]
    data.frame(df)

  } else {
    stop("Grouping column not provided or not present in GRanges object.")
  }
}


#'
#' Mean coverage of a BigWig list of files on a set of loci specified as BED file.
#'
#' Adapted from seqplots rtracklayer wrapping functions to compute coverage values:
#' For more on seqplots: https://www.bioconductor.org/packages/release/bioc/html/seqplots.html
#'
#' @param bwfile BigWig file to be summarized.
#' @param bed BED file file to be summarized.
#' @param mode How to handle differences in lengths across loci:
#'   stretch: Anchor each locus on both sides.
#'   start: Anchor all loci on start.
#'   end: Anchor all loci on end.
#'   center: Center all loci.
#' @param bin Bin size. Length of bin in base pairs. The lower, the higher the resolution.
#' @param upstream Number of base pairs to include upstream of loci.
#' @param downstream Number of base pairs to include downstream of loci.
#' @param ignore_strand Whether to use strand information in BED file.
#' @importFrom rtracklayer BigWigFile import
#' @return A DataFrame with the aggregated scores
calculate_bw_profile <- function(bw,
                                 gr,
                                 label=NULL,
                                 mode='stretch',
                                 bin=100,
                                 upstream=2500,
                                 downstream=2500,
                                 ignore_strand=F) {

  bwfile <- BigWigFile(path=bw)
  if (is.null(label)) {
    label <- basename(bw)
  }

  if (mode == 'stretch') {
    left_npoints <- floor(upstream/bin)
    right_npoints <- floor(downstream/bin)
    # Stretch to the median value of the GR object
    middle_npoints <- floor(median(GenomicRanges::width(gr))/bin )

    # So beautiful
    left <- summary_matrix(bwfile,
                           GenomicRanges::flank(gr, upstream, start=TRUE),
                           npoints=left_npoints,
                           ignore_strand=ignore_strand)

    right <- summary_matrix(bwfile,
                            GenomicRanges::flank(gr, downstream, start=FALSE),
                            npoints=right_npoints,
                            ignore_strand=ignore_strand)

    middle <- summary_matrix(bwfile,
                             gr,
                             npoints=middle_npoints,
                             ignore_strand=ignore_strand)

    full <- cbind(left, middle, right)

  } else {
    gr <- GenomicRanges::promoters(GenomicRanges::resize(gr, 1, fix=mode),
                                   upstream,
                                   downstream)

    npoints <- floor((upstream + downstream) / bin)

    full <- summary_matrix(bwfile, gr, npoints=npoints, ignore_strand=F)
  }

  make_averages_df(full, label)

}


#' Compute a coverage profile on a list of bw files with a given GRanges object.
#'
#'
#' @param bwfiles BigWig file list to be summarized.
#' @param bedfile BED file to summarize
#' @param colnames Names to be assigned to the columns
#' @return a data frame in long format
#' @export
bw_profile <- function(bwfiles,
                       bedfile,
                       colnames,
                       mode='stretch',
                       bin=100,
                       upstream=2500,
                       downstream=2500,
                       ignore_strand=F) {

  check_filelist(bwfiles)
  check_filelist(bedfile)
  gr <- rtracklayer::import(bedfile)

  values_list <- purrr::map2(bwfiles,
                             colnames,
                             calculate_bw_profile,
                             gr=gr,
                             mode=mode,
                             bin=bin,
                             upstream=upstream,
                             downstream=downstream,
                             ignore_strand=ignore_strand)

  values <- do.call(rbind, values_list)
  values
}


#' Summarize a BigWig file over loci on a GRanges object.
#'
#' @param track BigWigFile object (rtracklayer)
#' @param gr GRanges object
#' @param npoints How many points to take
#' @param ignore_strand Ignore strand information in gr. Default false.
#' @return A value matrix dimensions len(gr) x npoints
summary_matrix <- function(track, gr, npoints, ignore_strand=F) {
  values <- rtracklayer::summary(track,
                    which=gr,
                    as='matrix',
                    size=npoints)

  # Reverse minus strand rows
  if (!ignore_strand) {
    values[as.character(GenomicRanges::strand(gr))=='-',] <- values[
      as.character(GenomicRanges::strand(gr))=='-', ncol(values):1]
  }

  values
}


#' Compute averages and standard error values for a matrix returned by
#' summary_matrix.
#'
#' @param matrix A matrix returned by summary_matrix
#' @param label Label for the sample
#' @return A dataframe with summarized values, sderror and medians, plus a label.
make_averages_df <- function(matrix, label) {
  df <- data.frame(mean=colMeans(matrix, na.rm=TRUE),
                   sderror=apply(matrix, 2,
                                 function (n) { sd(n, na.rm=TRUE) / sqrt( sum(!is.na(n))) } ),
                   median=apply(matrix, 2, median, na.rm=TRUE))

  df$index <- as.integer(rownames(df))
  df$sample <- label
  df
}
