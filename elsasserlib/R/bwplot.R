#' Summary heatmap of a categorized BED file
#'
#' Make a summary heatmap where each cell contains an aggregated value of a
#' bigWig file from bwfiles and a category of a BED file (bedfile). The
#' provided BED file must have a name field that is valid (i.e. can be grouped,
#' representing some type of category, not a per-locus unique ID).
#'
#' @param labels Labels to use for in the plot for the bw files.
#' @param file_out Output the plot to a file.
#' @inheritParams bw_bed
#' @return A pheatmap object
#' @export
plot_bw_bed_summary_heatmap <- function(bwfiles,
                                        bedfile,
                                        bg_bwfiles = NULL,
                                        labels = NULL,
                                        aggregate_by = "true_mean",
                                        norm_func = identity,
                                        file_out = NA) {

  summary_values <- bw_bed(bwfiles, bedfile,
                      bg_bwfiles = bg_bwfiles,
                      aggregate_by = aggregate_by,
                      norm_func = norm_func,
                      labels = labels
                    )


  if (sum(summary_values) == 0) {
    warning("All zero-values matrix. Using same background as bw input?")
  }

  title <- paste("Coverage per region (", aggregate_by, ")")
  title <- paste(title, "-", make_norm_label(substitute(norm_func), bg_bwfiles))
  summary_heatmap(t(summary_values), title = title, file_out = file_out)
}

#' Bin-based scatterplot of a pair of bigWig files
#'
#' Plots a scatter plot from two given bigWig files and an optional set of BED
#' files as highlighted annotations. Bins are highlighted if there is at least
#' minoverlap base pairs overlap with any loci in BED file.
#'
#' If specifying minoverlap, you must take into account the bin_size parameter
#' and the size of the loci you are providing as BED file.
#'
#' Values in x and y axis can be normalized using background bigWig files
#' (usually input files). By default, the value shown will be x / bg_x per bin.
#' If norm_func_x or norm_func_y are provided, this can be changed to any given
#' function, for instance, if norm_func_x = log2, values on the x axis will
#' represent log2(x / bg_x) for each bin.
#'
#' Values that are invalid (NaN, Inf, -Inf) in doing such normalization will
#' be ignored and shown as warnings, as this is ggplot default behavior.
#'
#' @param x BigWig file corresponding to the x axis.
#' @param y BigWig file corresponding to the y axis.
#' @param bg_x BigWig file to be used as x axis background (us. input).
#' @param bg_y BigWig file to be used as y axis background (us. input).
#' @param highlight List of bed files to use as highlight for subgroups.
#' @param minoverlap Minimum overlap required for a bin to be highlighted
#' @param highlight_label Labels for the highlight groups.
#'  If not provided, filenames are used.
#' @param norm_func_x Function to use after x / x_bg.
#' @param norm_func_y Function to use after y / y_bg.
#' @param highlight_colors Array of color values for the highlighting groups
#' @import ggplot2
#' @inheritParams bw_bins
#' @return A ggplot object.
#' @export
plot_bw_bins_scatter <- function(x,
                                y,
                                bg_x = NULL,
                                bg_y = NULL,
                                bin_size = 10000,
                                per_locus_stat = "mean",
                                genome = "mm9",
                                highlight = NULL,
                                minoverlap = 0L,
                                highlight_label = NULL,
                                norm_func_x = identity,
                                norm_func_y = identity,
                                highlight_colors = NULL) {

  label_df <- function(df, name) {
    data.frame(df, group = name)
  }

  bins_values_x <- bw_bins(x,
                      bg_bwfiles = bg_x,
                      bin_size = bin_size,
                      genome = genome,
                      per_locus_stat = per_locus_stat,
                      norm_func = norm_func_x,
                      labels = "x"
                    )

  bins_values_y <- bw_bins(y,
                      bg_bwfiles = bg_y,
                      bin_size = bin_size,
                      genome = genome,
                      per_locus_stat = per_locus_stat,
                      norm_func = norm_func_y,
                      labels = "y"
                    )


  bins_df <- cbind(data.frame(bins_values_x),
                   data.frame(bins_values_y)[, "y"]
             )

  colnames(bins_df) <- c(colnames(bins_df)[seq_len(ncol(bins_df) - 1)], "y")
  bins_values <- makeGRangesFromDataFrame(bins_df, keep.extra.columns = TRUE)

  extra_plot <- NULL
  extra_colors <- NULL

  if (!is.null(highlight)) {
    gr_list <- lapply(highlight, rtracklayer::import, format = "BED")

    subset_func <- purrr::partial(
                     IRanges::subsetByOverlaps,
                     x = bins_values,
                     minoverlap = minoverlap
                   )

    bins_subset <- lapply(gr_list, subset_func)
    subset_df <- lapply(bins_subset, data.frame)

    if (is.null(highlight_label)) {
      highlight_label <- basename(highlight)
    }

    df_values_labeled <- purrr::map2(subset_df, highlight_label, label_df)
    highlight_values <- do.call(rbind, df_values_labeled)

    # Order of factors need to match to assign properly colors to points
    highlight_values$group <- factor(highlight_values$group,
                                     levels = highlight_label)

    extra_plot <- geom_point(
                    data = highlight_values,
                    aes_string(x = "x", y = "y", color = "group"),
                    alpha = 0.8
                  )

    if (!is.null(colors)) {
      extra_colors <- scale_color_manual(values=highlight_colors)
    }
  }

  x_label <- paste(make_label_from_filename(x), "-",
                   make_norm_label(substitute(norm_func_x), bg_x))

  y_label <- paste(make_label_from_filename(y), "-",
                   make_norm_label(substitute(norm_func_y), bg_y))

  ggplot(bins_df, aes_string(x = "x", y = "y")) +
    geom_point(color = "#cccccc", alpha = 0.8) +
    theme_elsasserlab_screen(base_size = 18) +
    xlab(x_label) +
    ylab(y_label) +
    ggtitle(paste("Bin coverage (bin_size = ", bin_size, ")", sep = "")) +
    extra_plot +
    extra_colors
}


#' Bin-based violin plot of a set of bigWig files
#'
#' Plots a violin plot of bin distribution of a set of bigWig files optionally
#' overlaid with annotated bins. Bins overlapping loci of the provided BED
#' file will be shown as a jitter plot on top of the violin plot.
#'
#' @param highlight BED file to use as highlight for subgroups.
#' @param minoverlap Minimum overlap required for a bin to be highlighted.
#' @param highlight_colors Array of color values for the highlighted groups.
#' @inheritParams bw_bins
#' @import ggplot2
#' @importFrom reshape2 melt
#' @return A ggplot object.
#' @export
plot_bw_bins_violin <- function(bwfiles,
                                bg_bwfiles = NULL,
                                labels = NULL,
                                bin_size = 10000,
                                per_locus_stat = "mean",
                                genome = "mm9",
                                highlight = NULL,
                                minoverlap = 0L,
                                norm_func = identity,
                                highlight_colors = NULL) {

  bins_values <- bw_bins(bwfiles,
                         bg_bwfiles = bg_bwfiles,
                         labels = labels,
                         bin_size = bin_size,
                         genome = genome,
                         per_locus_stat = per_locus_stat,
                         norm_func = norm_func)

  bins_df <- data.frame(bins_values)
  bwnames <- colnames(mcols(bins_values))
  bin_id <- c("seqnames", "start", "end")

  melted_bins <- melt(bins_df[, c(bin_id, bwnames)], id.vars = bin_id)
  title <- paste("Global bin distribution (binsize=", bin_size, ")", sep = "")
  extra_plot <- NULL
  extra_colors <- NULL

  if (!is.null(highlight)) {
    title <- paste(title, "vs", basename(highlight))
    gr_highlight <- rtracklayer::import(highlight, format = "BED")

    overlapping_values <- IRanges::subsetByOverlaps(
                            bins_values,
                            gr_highlight,
                            minoverlap = minoverlap
                          )

    overlap_df <- data.frame(overlapping_values)

    melted_highlight <- melt(overlap_df[, c(bin_id, bwnames)], id.vars = bin_id)

    extra_plot <- geom_jitter(data = melted_highlight,
                    aes_string(x = "variable", y = "value", color = "variable"),
                    alpha = 0.7
                  )

    if (!is.null(colors)) {
      extra_colors <- scale_color_manual(values = highlight_colors)
    }

  }

  y_label <- make_norm_label(substitute(norm_func), bg_bwfiles)

  # Avoid R_CMD_CHECK note on variable / value with aes
  ggplot(melted_bins, aes_string(x = "variable", y = "value")) +
    geom_violin(fill = "#cccccc") +
    theme_elsasserlab_screen(base_size = 18) +
    xlab("") +
    ylab(y_label) +
    ggtitle(title) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    extra_plot +
    extra_colors
}

#' Plot a heatmap of a given bigWig file over a set of loci
#'
#' @inheritParams plot_bw_profile
#' @param bwfile BigWig file to plot
#' @param bg_bwfile Background bw file. Use this with care. Depending on bin
#'   size and actual values, this may result in a very noisy plot.
#' @param zmin Minimum of the color scale. Majority of tools set
#'   this to 0.01 percentile of the data distribution.
#' @param zmax Maximum of the color scale. Majority of tools set
#'   this to 0.99.
#' @param cmap Color map. Any RColorBrewer palette name is accepted here.
#' @param max_rows_allowed Maximum number of loci that will be allowed in the
#'   plot. If the amount of loci exceeds this value, the plot will be binned
#'   on the y axis until it fits max_rows_allowed. This speeds up plotting of
#'   very large matrices, where higher resolution would not be perceivable by eye.
#' @importFrom dplyr group_by summarise
#' @export
plot_bw_heatmap <- function(bwfile,
                            bedfile,
                            bg_bwfile = NULL,
                            mode = "stretch",
                            bin_size = 100,
                            upstream = 2500,
                            downstream = 2500,
                            middle = NULL,
                            ignore_strand = FALSE,
                            norm_func = identity,
                            cmap = "Reds",
                            zmin = NULL,
                            zmax = NULL,
                            max_rows_allowed = 10000) {

  values <- bw_heatmap(bwfile, bedfile,
                       bg_bwfiles = bg_bwfile,
                       mode = mode,
                       bin_size = bin_size,
                       upstream = upstream,
                       downstream = downstream,
                       middle = middle,
                       ignore_strand = ignore_strand,
                       norm_func = norm_func
  )

  # Order matrix by mean and transpose it (image works flipped)
  m <- t(values[[1]][order(rowMeans(values[[1]]), decreasing=F),])

  zlim <- define_color_limits(m, zmin, zmax)

  zmin <- zlim[[1]]
  zmax <- zlim[[2]]

  # Cap values out of zlim
  m[m<zmin] <- zmin
  m[m>zmax] <- zmax

  df <- melt(m)
  colnames(df) <- c('x', 'y', 'value')

  df2 <- df

  if (ncol(m) > max_rows_allowed) {
    # Downsample rows only and downsample only enough to fit max_rows. So
    # we make sure we do not extremely downsample a value that only slightly
    # exceeds our max resolution.
    warning(paste0("Large matrix of ", ncol(m), ". Downscaling to ", max_rows_allowed))
    downsample_factor <- round(ncol(m) / max_rows_allowed)
    # downsample <- 10
    df2 <- df %>%
      dplyr::group_by(x = x,
                      y = downsample_factor * round(y / downsample_factor)) %>%
      dplyr::summarise(value = mean(value))
  }

  axis_breaks <- calculate_profile_breaks(nrow(m), upstream, downstream, bin_size, mode)
  axis_labels <- calculate_profile_labels(upstream, downstream, mode)

  lines <- axis_breaks[2]
  if (mode == "stretch") {
    lines <- axis_breaks[2:3]
  }

  loci <- ncol(m)
  y_label <- paste(basename(bedfile), "-", loci, "loci", sep = " ")
  x_title <- make_label_from_filename(bwfile)

  p <- ggplot(df2, aes(x=x, y=y, fill=value)) + geom_raster()

  gcol <- colorRampPalette(brewer.pal(n=8, name=cmap))
  p <- p +
    scale_x_continuous(breaks = axis_breaks, labels = axis_labels, expand=c(0,0)) +
    # scale_y_continuous(expand=c(0,0)) +
    scale_y_continuous(breaks = c(1, loci), labels = c(loci, "0"), expand=c(0,0)) +
    scale_fill_gradientn(
      colours = gcol(100), limits = c(zmin, zmax),
      breaks=c(zmin, zmax), labels=format(c(zmin, zmax), digits=2)
    ) +
    xlab(x_title) +
    ylab(y_label) +
    geom_vline(
      xintercept = lines,
      linetype = "dashed",
      color = "#111111",
      size = 0.2
    ) +
    ggtitle("Heatmap plot") +
    theme_elsasserlab_screen(base_size = 14) +
    theme(axis.line = element_blank(),
          panel.border = element_rect(color = "black", fill=NA, size = 0.1)) +
    labs(fill = make_norm_label(substitute(norm_func), bg_bwfile))

  p
}

#' Profile plot of a set of bigWig files
#'
#' Plots a profile of a set of bigWig files over a set of loci in a BED file.
#'
#' @param show_error Show standard error.
#' @param colors Array of colors that will  be assigned to labels or files
#'    (in that order)
#' @inheritParams bw_profile
#' @import ggplot2
#' @return A ggplot object.
#' @export
plot_bw_profile <- function(bwfiles,
                            bedfile,
                            bg_bwfiles = NULL,
                            mode = "stretch",
                            bin_size = 100,
                            upstream = 2500,
                            downstream = 2500,
                            middle = NULL,
                            ignore_strand = FALSE,
                            show_error = FALSE,
                            norm_func = identity,
                            labels = NULL,
                            colors = NULL) {

  values <- bw_profile(bwfiles, bedfile,
              bg_bwfiles = bg_bwfiles,
              mode = mode,
              bin_size = bin_size,
              upstream = upstream,
              downstream = downstream,
              middle = middle,
              ignore_strand = ignore_strand,
              norm_func = norm_func,
              labels = labels
            )

  y_label <- make_norm_label(substitute(norm_func), bg_bwfiles)

  nrows <- max(values$index)

  axis_breaks <- calculate_profile_breaks(nrows, upstream, downstream, bin_size, mode)
  axis_labels <- calculate_profile_labels(upstream, downstream, mode)

  lines <- axis_breaks[2]
  if (mode == "stretch") {
    lines <- axis_breaks[2:3]
  }

  loci <- length(import(bedfile, format = "BED"))
  x_title <- paste(basename(bedfile), "-", loci, "loci", sep = " ")

  p <- ggplot(values,
              aes_string(
                x = "index",
                y = "mean",
                color = "sample",
                fill = "sample"
              )) +
    geom_line(size = 0.8) +
    geom_vline(
      xintercept = lines,
      linetype = "dashed",
      color = "#cccccc",
      alpha = 0.8
    ) +
    scale_x_continuous(breaks = axis_breaks,
                       labels = axis_labels,
                       limits = c(0.5, nrows + 0.5)) +
    xlab(x_title) +
    ylab(y_label) +
    ggtitle("Profile plot") +
    theme_elsasserlab_screen(base_size = 18) +
    theme(
      legend.position = c(0.80, 0.90),
      legend.direction = "vertical",
      legend.title = element_blank()
    )
  if (show_error) {
    p <-
      p + geom_ribbon(
        aes(
          x = index,
          ymin = mean - sderror,
          ymax = mean + sderror
        ),
        color = NA,
        alpha = 0.3
      )
  }

  if (!is.null(colors)) {
    p <- p + scale_color_manual(values=colors)
  }

  p
}


#' Plot a pretty heatmap using pheatmap library
#'
#' This function ignores NA values to calculate min and max values.
#'
#' @param values Value matrix.
#' @param title Plot title.
#' @param size Square size in points.
#' @param file_out Optional file output.
#' @importFrom RColorBrewer brewer.pal
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette
#' @return pheatmap object
summary_heatmap <- function(values, title, size = 35, file_out = NA) {
  bcolor <- "white"

  breakslist <- calculate_breakslist(values)
  html_colors <- rev(RColorBrewer::brewer.pal(n = 7, name = "RdBu"))
  palette <- colorRampPalette(html_colors)(length(breakslist))

  cellsize_inches <- size / 72

  margin <- 4
  plot_width_inches <- cellsize_inches * ncol(values) + margin
  plot_height_inches <- cellsize_inches * nrow(values) + margin

  pheatmap(
    values,
    main = title,
    cellwidth = size,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    cellheight = size,
    border_color = bcolor,
    breaks = breakslist,
    color = palette,
    width = plot_width_inches,
    height = plot_height_inches,
    display_numbers = TRUE,
    filename = file_out
  )
}


#' Calculate breaks list for the values in a matrix
#'
#' Given a numerical matrix, calculate breakslist to use in the summarized
#' heatmap plot. This is done in order to unify color scale regardless of
#' normalization used (mostly log vs linear), keeping white as the zero value.
#'
#' This function ignores NA values to calculate min and max values.
#'
#' @param values Value matrix.
#' @return Sequence of breaks used by heatmap function
calculate_breakslist <- function(values) {
  # Compute the largest deviation from zero ignoring NAs
  maxmat <- values
  maxmat[is.na(maxmat)] <- -Inf
  maxval <- max(maxmat)

  minmat <- values
  minmat[is.na(minmat)] <- Inf
  minval <- min(minmat)

  breaklim <- ceiling(abs(max(abs(c(minval, maxval)))))

  # Compute a reasonable amount of steps
  nsteps <- 21.0
  stepsize <- 2 * breaklim / nsteps

  breakslist <- seq(-breaklim, +breaklim, by = stepsize)
  breakslist
}


calculate_profile_breaks <- function(nrows, upstream, downstream, bin_size, mode) {
  upstream_nbins <- floor(upstream / bin_size)
  downstream_nbins <- floor(downstream / bin_size)

  # index value starts at 1
  axis_breaks <- c(1, upstream_nbins + 1, nrows + 1)

  # Put ticks on the edges
  axis_breaks <- axis_breaks - 0.5

  if (mode == "stretch") {
    axis_breaks <- c(1, upstream_nbins + 1, nrows - downstream_nbins + 1, nrows + 1)
    # center ticks on the middle of the bins
    axis_breaks <- axis_breaks - 0.5
  }

  axis_breaks
}

calculate_profile_labels <- function(upstream, downstream, mode) {
  if (mode == "stretch") {
    c(paste("-", upstream / 1000, "kb", sep = ""),
      "start", "end",
      paste("+", downstream / 1000, "kb", sep = ""))

  } else {
    c(paste("-", upstream / 1000, "kb", sep = ""),
      mode,
      paste("+", downstream / 1000, "kb", sep = ""))
  }
}


#' Define color limits from a value matrix and provided parameters
#'
#' @param m Value matrix.
#' @param zmin Min value. Overrides percentile.
#' @param zmax Max value. Overrides percentile.
#'
#' @return A pair of c(min, max)
define_color_limits <- function(m, zmin, zmax) {
  # colorscale limits percentiles: 0.01 - 0.99
  zlim <- quantile(unlist(m), c(0.01, 0.99), na.rm=TRUE)

  if (!is.null(zmin)) {
    zlim[[1]] <- zmin
  }

  if (!is.null(zmax)) {
    zlim[[2]] <- zmax
  }
  zlim
}


#' Generate a human-readable normalization function string
#'
#' @param f String representing normalization function.
#' @param bg Background file.
#'
#' @return A string describing normalization.
make_norm_label <- function(f, bg) {
  label <- "RPGC"
  if (!is.null(bg)) {
    if (f != "identity") {
      label <- paste(f, "(", label, " / background)", sep = "")
    } else {
      label <- paste(label, " / background", sep = "")
    }
  }
  label
}
