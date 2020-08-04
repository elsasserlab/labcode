#' Normalize a column of a dataframe, replacing it.
#'
#' @param df Data frame to be used.
#' @param fg_col Column that corresponds to main signal.
#' @param bg_col Column that corresponds to background.
#' @param norm_func Function to apply after dividing.
norm_col <- function(df, fg_col, bg_col, norm_func) {
  cols <- colnames(df)
  cols <- cols[! cols %in% c(fg_col)]
  new_col <- norm_func(df[, fg_col] / df[, bg_col])
  df[, fg_col] <- new_col
  data.frame(df)
}


#' Basic scatterplot from two bigwig files and an optional set of BED
#' files as highlighted annotations.
#'
#'
#' @param x BigWig file to be used for the x axis
#' @param y BigWig file to be used for the y axis
#' @param x_bg BigWig file to be used for the x axis background (us. input)
#' @param y_bg BigWig file to be used for the y axis background (us. input)
#' @param bsize Bin size. Default 10000.
#' @param stat Statistics used in the calculation of bins.
#' @param genome Genome. Available choices are mm9, hg38.
#' @param highlight List of bed files to use as highlight for subgroups.
#' @param minoverlap Minimum overlap required for a bin to be highlighted
#' @param highlight_label Labels for the highlight groups.
#'  If not provided, filenames are used.
#' @param x_func Function to use after x / x_bg
#' @param y_func Function to use after y / y_bg
#' @return A ggplot object
#' @import ggplot2
#' @export
bw_bins_scatterplot <- function(x,
                                y,
                                x_bg=NULL,
                                y_bg=NULL,
                                bsize=10000,
                                stat='mean',
                                genome='mm9',
                                highlight=NULL,
                                minoverlap=0L,
                                highlight_label=NULL,
                                x_func=identity,
                                y_func=identity) {

  label_df <- function(df, name) {
    data.frame(df, group=name)
  }

  names <- make.names(basename(c(x, y, x_bg, y_bg)))
  bins_values <- bw_bins(c(x, y, x_bg, y_bg),
                         bsize=bsize,
                         genome=genome,
                         stat=stat)

  xcol <- names[[1]]
  ycol <- names[[2]]

  bins_df <- data.frame(bins_values)

  if (!is.null(x_bg)) {
    x_bgcol <- make.names(basename(x_bg))
    bins_df <- norm_col(bins_df, xcol, x_bgcol, x_func)
  }

  if (!is.null(y_bg)) {
    y_bgcol <- make.names(basename(y_bg))
    bins_df <- norm_col(bins_df, ycol, y_bgcol, y_func)
  }

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

    subset_values <- lapply(gr_list, subset_bins)
    subset_df <- lapply(subset_values, data.frame)

    if (!is.null(x_bg)) {
      x_bgcol <- make.names(basename(x_bg))
      subset_df <- lapply(subset_df, norm_col, fg_col=xcol, bg_col=x_bgcol, norm_func=x_func)
    }

    if (!is.null(y_bg)) {
      y_bgcol <- make.names(basename(y_bg))
      subset_df <- lapply(subset_df, norm_col, fg_col=ycol, bg_col=y_bgcol, norm_func=y_func)
    }

    if (is.null(highlight_label)) {
      highlight_label <- basename(highlight)
    }

    df_values_labeled <- purrr::map2(subset_df, highlight_label, label_df)
    highlight_values <- do.call(rbind, df_values_labeled)
    p <- p + geom_point(data=highlight_values, aes_string(x=xcol,
                                                          y=ycol,
                                                          color='group'), alpha=0.8)

  }

  p
}


#' Violin plot of bin distribution of a set of bigWig files overlayed with
#' annotated bins (i.e. bins overlapping a given BED file)
#'
#'
#' @param bw BigWig files
#' @param bg_bw BigWig files used as background (us. input)
#' @param bw_label Labels to use for in the plot for the bw files.
#' @param bsize Bin size. Default 10000.
#' @param stat Statistics used in the calculation of bins.
#' @param genome Genome. Available choices are mm9, hg38.
#' @param highlight Bed file to use as highlight for subgroups.
#' @param minoverlap Minimum overlap required for a bin to be highlighted
#' @param norm_func Function to use on top of dividing bw/bg_bw (usually identity or log2)
#' @return A ggplot object
#' @import ggplot2
#' @importFrom reshape2 melt
#' @export
bw_bins_violinplot <- function(bw,
                               bg_bw=NULL,
                               bw_label=NULL,
                               bsize=10000,
                               stat='mean',
                               genome='mm9',
                               highlight=NULL,
                               minoverlap=0L,
                               norm_func=identity) {

  if (!is.null(bg_bw) && length(bg_bw) != length(bw)) {
    stop('BigWig file list and background bw file list must have the same length')
  }

  if (!is.null(bw_label) && length(bw) != length(bw_label)) {
    stop('colnames and bw lists must have the same length')
  }

  bins_values <- bw_bins(c(bw, bg_bw),
                         bsize=bsize,
                         genome=genome,
                         stat=stat)

  bins_df <- data.frame(bins_values)
  bwnames <- make.names(basename(bw))

  if (!is.null(bg_bw)) {
    bgnames <- make.names(basename(bg_bw))
    bins_df[, bwnames] <- norm_func(bins_df[, bwnames] / bins_df[, bgnames])
  }

  bin_id <- c('seqnames', 'start', 'end')
  bins_df <- bins_df[, c(bin_id, bwnames)]

  if (!is.null(bw_label)) {
    colnames(bins_df) <- c(bin_id, bw_label)
  }

  melted_bins <- melt(bins_df, id.vars=bin_id)

  title <- paste('Global bin distribution (binsize=', bsize, ')', sep='')
  if (!is.null(highlight)) {
    title <- paste(title, 'vs', basename(highlight))
  }

  p <- ggplot(melted_bins, aes(x=variable, y=value)) +
    geom_violin(fill='#cccccc') +
    theme_elsasserlab_screen(base_size = 18) +
    xlab('') +
    ylab('Mean coverage') +
    ggtitle(title) +
    theme(legend.position='none', axis.text.x=element_text(angle=45, hjust=1))

  if (!is.null(highlight)) {
    gr_highlight <- rtracklayer::import(highlight, format='BED')

    overlapping_values <- IRanges::subsetByOverlaps(bins_values,
                                                    gr_highlight,
                                                    minoverlap=minoverlap)

    overlapping_df <- data.frame(overlapping_values)
    overlapping_df_norm <- cbind(overlapping_df[, bin_id],
                                 norm_func(overlapping_df[, bwnames] / overlapping_df[, bgnames]))

    if (!is.null(bw_label)) {
      colnames(overlapping_df_norm) <- c(bin_id, bw_label)
    }

    melted_highlight <- melt(overlapping_df_norm[, c(bin_id, bw_label)], id.vars=bin_id)

    p <- p + geom_jitter(data=melted_highlight,
                         aes(x=variable, y=value, color=variable),
                         alpha=0.7)
  }

  p
}