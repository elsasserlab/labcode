#' Print a basic scatterplot from two bigwig files and an optional set of BED
#' files as highlighted annotations.
#'
#'
#' @param xtrack BigWig file to be used for the x axis
#' @param ytrack BigWig file to be used for the y axis
#' @param bsize Bin size. Default 10000.
#' @param stat Statistics used in the calculation of bins.
#' @param genome Genome. Available choices are mm9, hg38.
#' @param highlight List of bed files to use as highlight for subgroups.
#' @param minoverlap Minimum overlap required for a bin to be highlighted
#' @param highlight_label Labels for the highlight groups.
#' @param xlimits Limits used for x axis
#' @param ylimits Limits used for y axis
#' @param xlabel Label for x axis
#' @param ylabel Label used for y axis
#'  If not provided, filenames are used.
#' @return A ggplot object
#' @import ggplot2
#' @export

bw_bins_scatterplot <- function(xtrack,
                                ytrack,
                                bsize=10000,
                                stat='mean',
                                genome='mm9',
                                highlight=NULL,
                                minoverlap=0L,
                                highlight_label=NULL,
                                xlimits=NULL,
                                ylimits=NULL,
                                xlabel=NULL,
                                ylabel=NULL) {

  label_df <- function(df, name) {
    data.frame(df, group=name)
  }

  names <- make.names(basename(c(xtrack, ytrack)))
  bins_values <- bw_bins(c(xtrack, ytrack),
                         bsize=bsize,
                         genome=genome,
                         stat=stat)

  xcol <- names[[1]]
  ycol <- names[[2]]

  bins_df <- data.frame(bins_values)
  p <- ggplot(bins_df, aes_string(x=xcol, y=ycol)) +
    geom_point(color='#cccccc', alpha=0.8) +
        theme_elsasserlab_screen(base_size = 18) +
    xlab(xcol) +
    ylab(ycol) +
    ggtitle(paste('Bin coverage (bsize = ', bsize, ')', sep=''))

  if (!is.null(highlight)) {
    gr_list <- lapply(highlight, rtracklayer::import, format='BED')

    subset_bins <- purrr::partial(IRanges::subsetByOverlaps,
                                  x=bins_values,
                                  minoverlap=minoverlap)

    gr_values <- lapply(gr_list, subset_bins)
    df_values <- lapply(gr_values, data.frame)

    if (is.null(highlight_label)) {
      highlight_label <- basename(highlight)
    }

    df_values_labeled <- purrr::map2(df_values, highlight_label, label_df)
    highlight_values <- do.call(rbind, df_values_labeled)
    p <- p + geom_point(data=highlight_values, aes_string(x=xcol,
                                                          y=ycol,
                                                          color='group'), alpha=0.8)

  }

  if (!is.null(ylimits)) {
    p <- p + ylim(ylimits)
  }

  if (!is.null(xlimits)) {
    p <- p + xlim(xlimits)
  }

  if (!is.null(xlabel)) {
    p <- p + xlab(xlabel)
  }

  if (!is.null(ylabel)) {
    p <- p + ylab(ylabel)
  }
  p
}
