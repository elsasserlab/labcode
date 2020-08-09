#' Basic scatterplot from two bigwig files and an optional set of BED
#' files as highlighted annotations.
#'
#'
#' @param x BigWig file to be used for the x axis
#' @param y BigWig file to be used for the y axis
#' @param bg_x BigWig file to be used for the x axis background (us. input)
#' @param bg_y BigWig file to be used for the y axis background (us. input)
#' @param bsize Bin size. Default 10000.
#' @param stat Statistics used in the calculation of bins.
#' @param genome Genome. Available choices are mm9, hg38.
#' @param highlight List of bed files to use as highlight for subgroups.
#' @param minoverlap Minimum overlap required for a bin to be highlighted
#' @param highlight_label Labels for the highlight groups.
#'  If not provided, filenames are used.
#' @param norm_func_x Function to use after x / x_bg
#' @param norm_func_y Function to use after y / y_bg
#' @return A ggplot object
#' @import ggplot2
#' @export
bw_bins_scatterplot <- function(x,
                                y,
                                bg_x=NULL,
                                bg_y=NULL,
                                bsize=10000,
                                stat='mean',
                                genome='mm9',
                                highlight=NULL,
                                minoverlap=0L,
                                highlight_label=NULL,
                                norm_func_x=identity,
                                norm_func_y=identity) {

  label_df <- function(df, name) {
    data.frame(df, group=name)
  }

  bins_values_x <- bw_bins(x,
                      bg_bwfiles=bg_x,
                      bsize=bsize,
                      genome=genome,
                      stat=stat,
                      norm_func=norm_func_x,
                      colnames='x')

  bins_values_y <- bw_bins(y,
                      bg_bwfiles=bg_y,
                      bsize=bsize,
                      genome=genome,
                      stat=stat,
                      norm_func=norm_func_y,
                      colnames='y')


  bins_df <- cbind(data.frame(bins_values_x), data.frame(bins_values_y)[,'y'])
  colnames(bins_df) <- c(colnames(bins_df)[1:ncol(bins_df)-1], 'y')
  bins_values <- makeGRangesFromDataFrame(bins_df, keep.extra.columns = T)

  extra_plot <- NULL

  if (!is.null(highlight)) {
    gr_list <- lapply(highlight, rtracklayer::import, format='BED')

    subset_func <- purrr::partial(IRanges::subsetByOverlaps,
                                  x=bins_values,
                                  minoverlap=minoverlap)

    bins_subset <- lapply(gr_list, subset_func)
    subset_df <- lapply(bins_subset, data.frame)

    if (is.null(highlight_label)) {
      highlight_label <- basename(highlight)
    }

    df_values_labeled <- purrr::map2(subset_df, highlight_label, label_df)
    highlight_values <- do.call(rbind, df_values_labeled)

    extra_plot <- geom_point(data=highlight_values,
                             aes_string(x='x', y='y', color='group'),
                             alpha=0.8)
  }

  p <- ggplot(bins_df, aes_string(x='x', y='y')) +
    geom_point(color='#cccccc', alpha=0.8) +
    theme_elsasserlab_screen(base_size = 18) +
    xlab(make.names(basename(x))) +
    ylab(make.names(basename(y))) +
    ggtitle(paste('Bin coverage (bsize = ', bsize, ')', sep='')) +
    extra_plot

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

  bins_values <- bw_bins(bw,
                         bg_bwfiles=bg_bw,
                         colnames=bw_label,
                         bsize=bsize,
                         genome=genome,
                         stat=stat)

  bins_df <- data.frame(bins_values)
  bwnames <- colnames(mcols(bins_values))
  bin_id <- c('seqnames', 'start', 'end')

  melted_bins <- melt(bins_df[, c(bin_id, bwnames)], id.vars=bin_id)
  title <- paste('Global bin distribution (binsize=', bsize, ')', sep='')
  extra_plot <- NULL

  if (!is.null(highlight)) {
    title <- paste(title, 'vs', basename(highlight))
    gr_highlight <- rtracklayer::import(highlight, format='BED')

    overlapping_values <- IRanges::subsetByOverlaps(bins_values,
                                                    gr_highlight,
                                                    minoverlap=minoverlap)

    overlapping_df <- data.frame(overlapping_values)

    melted_highlight <- melt(overlapping_df[, c(bin_id, bwnames)], id.vars=bin_id)

    extra_plot <- geom_jitter(data=melted_highlight,
                                  aes_string(x='variable', y='value', color='variable'),
                                  alpha=0.7)
  }

  ggplot(melted_bins, aes_string(x='variable', y='value')) +
    geom_violin(fill='#cccccc') +
    theme_elsasserlab_screen(base_size = 18) +
    xlab('') +
    ylab('Mean coverage') +
    ggtitle(title) +
    theme(legend.position='none', axis.text.x=element_text(angle=45, hjust=1)) +
    extra_plot
}


#' Violin plot of bin distribution of a set of bigWig files overlayed with
#' annotated bins (i.e. bins overlapping a given BED file)
#'
#' @param bw BigWig files.
#' @param bed BED file to use to summarize. It needs to have an adequate `name` field
#'   (where names correspond to categories that can be grouped).
#' @param bg_bw BigWig files used as background (us. input).
#' @param bw_label Labels to use for in the plot for the bw files.
#' @param aggregate_by Can be true_mean, mean (mean of means), median (median of means).
#' @param norm_func Function to use on top of dividing bw/bg_bw (usually identity or log2).
#' @param file_out Output the plot to a file.
#' @return A plot object
#' @export
bw_bed_summary_heatmap <- function(bw,
                                   bed,
                                   bg_bw=NULL,
                                   bw_label=NULL,
                                   aggregate_by='true_mean',
                                   norm_func=identity,
                                   file_out=NA) {

  title <- paste('Coverage per region (', aggregate_by, ')')
  summary_values <- bw_bed(bw,
                           bed,
                           bg_bwfiles = bg_bw,
                           aggregate.by=aggregate_by,
                           norm_func=norm_func)

  if (!is.null(bg_bw)) {
    title <- paste(title, '- Norm to background')
  }
  print(summary_values)
  summary_heatmap(t(summary_values), title=title, file_out=file_out)
}


#' Plot a pretty heatmap using pheatmap library
#'
#' This function ignores NA values to calculate min and max values.
#'
#' @param values Value matrix
#' @param title Title of the plot
#' @param size Size of the square
#' @param file_out Optionally save directly to a file. This will take care of
#'   cropping the heatmap to a right size so it fits.
#' @importFrom RColorBrewer brewer.pal
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @return Heatmap plot
summary_heatmap <- function(values, title, size=35, file_out=NA) {
  bcolor <- "white"

  breakslist <- compute_breakslist(values)
  palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(length(breakslist))

  cellsize_inches <- size / 72

  margin <- 4
  plot_width_inches <- cellsize_inches*ncol(values) + margin
  plot_height_inches <- cellsize_inches*nrow(values) + margin

  plot <- pheatmap(values,
                   main=title,
                   cellwidth=size,
                   cluster_rows=F,
                   cluster_cols=F,
                   cellheight=size,
                   border_color=bcolor,
                   breaks=breakslist,
                   color=palette,
                   width=plot_width_inches,
                   height=plot_height_inches,
                   display_numbers=TRUE,
                   filename = file_out)

  plot
}

#' Given a value matrix, compute breakslist for the summarized heatmap.
#'
#' This is done so log scale and norm scale heatmaps have the same color scale
#' on the positive values, and zero values are shown as white always.
#'
#' This function ignores NA values to calculate min and max values.
#' @param mat Value matrix
#' @return Sequence of breaks used by heatmap function
compute_breakslist <- function(mat) {
  # Compute the largest deviation from zero ignoring NAs
  maxmat <- mat
  maxmat[is.na(maxmat)] <- -Inf
  maxval <- max(maxmat)

  minmat <- mat
  minmat[is.na(minmat)] <- Inf
  minval <- min(minmat)

  breaklim <- ceiling(abs(max(abs(c(minval, maxval)))))

  # Compute a reasonable amount of steps
  nsteps <- 21.0
  stepsize <- 2*breaklim / nsteps

  breakslist <- seq(-breaklim, +breaklim, by=stepsize)
  breakslist
}


#' Profile plot of a set of bw files
#' annotated bins (i.e. bins overlapping a given BED file)
#'
#' @param bw BigWig files
#' @param bed BED file to use to summarize. It needs to have an adequate `name` field
#'   (where names correspond to categories that can be grouped).
#' @param bg_bw BigWig files used as background (us. input)
#' @param mode How to align BED loci: Can be stretch, start, end, center.
#' @param bsize Bin size.
#' @param upstream Number of base pairs to include upstream of loci.
#' @param downstream Number of base pairs to include downstream of loci.
#' @param ignore_strand Ignore strand information in BED file (default false).
#' @param show_error Show standard error.
#' @param norm_func Function to apply after bw / bg (default identity)
#' @import ggplot2
#' @return A plot object
#' @export
bw_profile_plot <- function(bw,
                            bed,
                            bg_bw=NULL,
                            mode='stretch',
                            bsize=100,
                            upstream=2500,
                            downstream=2500,
                            ignore_strand=F,
                            show_error=F,
                            norm_func=identity) {

  values <- bw_profile(bw,
                       bed,
                       bg_bwfiles=bg_bw,
                       mode=mode,
                       bin=bsize,
                       upstream=upstream,
                       downstream=downstream,
                       ignore_strand=ignore_strand,
                       norm_func=identity)

  ylabel <- 'RPGC'

  if (!is.null(bg_bw)) {
    ylabel <- paste('Enrichment over background')
  }

  nrows <- max(values$index)
  left_flank_size <- floor(upstream/bsize)
  right_flank_size <- floor(downstream/bsize)

  # index value starts at 1
  axis_breaks <- c(1, left_flank_size, nrows)

  axis_labels <- c(paste('-', upstream/1000, 'kb', sep=''),
                   mode,
                   paste('+', downstream/1000, 'kb', sep=''))

  lines <- axis_breaks[2]

  if (mode == 'stretch') {
    axis_breaks <- c(1,
                     left_flank_size,
                     nrows-right_flank_size,
                     nrows)

    axis_labels <- c(paste('-', upstream/1000, 'kb', sep=''),
                     'start',
                     'end',
                     paste('+', downstream/1000, 'kb', sep=''))

    lines <- axis_breaks[2:3]

  }

  loci <- length(import(bed, format='BED'))
  xtitle <- paste(basename(bed), '-', loci, 'loci', sep=' ')

  p <- ggplot(values, aes_string(x='index', y='mean', color='sample', fill='sample')) +
    geom_line(size=0.8) +
    geom_vline(xintercept=lines, linetype='dashed', color='#cccccc', alpha=0.8) +
    scale_x_continuous(breaks=axis_breaks, labels=axis_labels) +
    xlab(xtitle) +
    ylab(ylabel) +
    ggtitle('Profile plot') +
    theme_elsasserlab_screen(base_size=18) +
    theme(legend.position=c(0.80,0.90), legend.direction = 'vertical', legend.title = element_blank())
  if (show_error) {
    p <- p + geom_ribbon(aes(x=index, ymin=mean-sderror, ymax=mean+sderror), color=NA, alpha=0.3)
  }

  p
}


