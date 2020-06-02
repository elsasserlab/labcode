
#' Build a binned-scored GRanges object from a bigWig file
#'
#' Build a binned-scored GRanges object from a bigWig file. The aggregating
#' function can be min, max, sd, mean.
#'
#' @param bwfile BigWig file to be summarized (or list)
#' @param colnames List of names to give to the mcols of the returned GRanges object
#' @param stat Aggregating function (per locus). Mean by default. Choices: min, max, sd, mean.
#' @param bsize Bin size. Default 10000.
#' @param genome Genome. Available choices are mm9, hg38.
#' @return A GenomicRanges object with each bwfile as a metadata column named after colnames.
#' @export
bw_bins <- function(bwfiles, colnames, stat='mean', bsize=10000, genome='mm9') {
  tiles <- build_bins(bsize=bsize, genome=genome)
  if (length(bwfiles) != length(colnames)) {
    stop("BigWig file list and column names must have the same length.")
  }
  multi_bw_ranges(bwfiles, colnames, tiles, per.locus.stat=stat)
}

#' Intersect a list of bw files with a GRanges object
#'
#' Build a binned-scored GRanges object from a list of bigWig files. The aggregating
#' function can be min, max, sd, mean.
#'
#' @param bwfilelist BigWig file to be summarized.
#' @param col.names Names to be assigned to the columns
#' @param gr GRanges object to intersect
#' @param stat Aggregating function per stat
#' @importFrom GenomicRanges makeGRangesFromDataFrame
multi_bw_ranges <- function(bwfilelist, colnames, gr, per.locus.stat='mean') {
  summaries <- purrr::map(bwfilelist, bw_ranges, gr=gr, per.locus.stat=per.locus.stat)
  with.names <- purrr::map2(summaries, colnames, rename_score)
  df.fg <- Reduce(function(...) merge(..., all=TRUE), with.names)
  if (length(bwfilelist) > 1) {
    df.fg <- df.fg %>% dplyr::arrange(seqnames, start)
  }

  makeGRangesFromDataFrame(df.fg, keep.extra.columns = T)
}


multi_bw_bed <- function(bwfiles, colnames, bedfile, per.locus.stat='mean', aggregate.by=NULL) {
  bed <- import(bedfile)
  result <- multi_bw_ranges(bwfiles, colnames, gr=bed, per.locus.stat=per.locus.stat)
  if ( 'name' %in% names(mcols(bed)) ) {
    result$name <- bed$name
  }
  if (! is.null(aggregate.by)) {
    df <- aggregate_scores(result, group.col='name', aggregate.by=aggregate.by)
    result <- df
  }

  result
}


#' Rename score function of a GRanges object
#'
#' @param gr GRanges object
#' @param new.name Name to give to the column.
rename_score <- function(gr, new.name) {
  colnames(mcols(gr)) <- replace(colnames(mcols(gr)), colnames(mcols(gr))=='score', new.name)
  gr
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
#' @param per.locus.stat Aggregating function (per locus). Mean by default. Choices: min, max, sd, mean. These choices depend on rtracklayer library.
#' @importFrom rtracklayer BigWigFile
#' @export
#' @return Data frame with columns score and group.col (if provided).
bw_ranges <- function (bwfile, gr, per.locus.stat='mean') {
  bw <- BigWigFile(bwfile)
  explicit_summary <- getMethod("summary", "BigWigFile")
  unlist(explicit_summary(bw, gr, type=per.locus.stat))
}

#'
#' #' Build a scored GRanges object from a BED file.
#' #'
#' #' Build a scored GRanges object from a bigWig file and a BED file.
#' #' The aggregating function (per locus) can be min, max, sd, mean.
#' #'
#' #' @param bwfile BigWig file to be summarized.
#' #' @param bedfile BED file to intersect with the BigWig file.
#' #' @param per.locus.stat Aggregate per locus
#' #' @param aggregate.by Statistic to aggregate per group. If NULL, values are not aggregated. This is the behavior by default.
#' #' @param keep.name Keep the name of specific locus (adds one mdata col to the GRanges object). True by default.
#' #' @export
#' #' @importFrom rtracklayer import BigWigFile
#' bw_bed <- function(bwfile, bedfile, per.locus.stat='mean', aggregate.by=NULL) {
#'   bed <- import(bedfile)
#'   result <- bw_ranges(bwfile, bed, per.locus.stat=per.locus.stat)
#'   if ( 'name' %in% names(mcols(bed)) ) {
#'     result$name <- bed$name
#'   }
#'   if (! is.null(aggregate.by)) {
#'     df <- aggregate_scores(result[,c('score','name')], group.col='name', aggregate.by=aggregate.by)
#'     result <- df
#'   }
#'
#'   result
#' }

#'
#' Aggregate scores of a GRanges object on a specific field
#'
#' Aggregates scores of a GRanges object on a specific field.
#' @param scored.gr A GRanges object with numerical metadata columns
#' @param group.col A column among the mcols that can be seen as a factor (usually name)
#' @param aggregate.by Function used to aggregate (mean, median, any valid function)
#' @return A DataFrame with the aggregated scores (any numerical column will be aggregated).
#' @importFrom dplyr group_by_ summarise across `%>%`
aggregate_scores <- function(scored.gr, group.col, aggregate.by) {
  df <- data.frame(mcols(scored.gr))
  if ( !is.null(group.col) && group.col %in% names(mcols(scored.gr))) {
    df <- df %>%
      group_by_(group.col) %>%
      summarise(across(where(is.numeric), aggregate.by))
  }
  as.data.frame(df)
}



