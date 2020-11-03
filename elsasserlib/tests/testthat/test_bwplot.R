context("Test functions for bigWig plots")
library(GenomicRanges)
library(testthat)
library(mockery)

get_file_path <- function(filename) {
  system.file("extdata", filename, package = "elsasserlib")
}

bw1 <- get_file_path("ES_H33_00h_rep1_hoxc.bw")
bw2 <- get_file_path("ES_H33_03h_rep1_hoxc.bw")
bg_bw <- get_file_path("ES_H33_inp_rep1_hoxc.bw")
bed <- get_file_path("chromhmm_hoxc.bed")

bw_limits <- GRanges(seqnames = c("chr15"),
                     ranges = IRanges(c(102723600, 102959000)))

reduced_bins <- bw_bins(bw1, selection = bw_limits, labels = "x")
reduced_bins_2 <- bw_bins(bw2, selection = bw_limits, labels = "y")
summary_values <- bw_bed(c(bw1, bw2), bed, aggregate_by = "mean")
profile_values <- bw_profile(bw1,
                    bedfile = bed,
                    upstream = 1000,
                    downstream = 1000
                  )

test_that("Setup files exist", {
  expect_true(file_test("-f", bw1))
  expect_true(file_test("-f", bw2))
  expect_true(file_test("-f", bg_bw))
  expect_true(file_test("-f", bed))
})

test_that("plot_bw_bins_scatter with defaults returns a plot", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_bins_scatter(bw1, bw2)
    expect_is(p, "ggplot")
  })
})

test_that("plot_bw_bins_scatter with highlight set returns a plot", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_bins_scatter(bw1, bw2, highlight = bed)
    expect_is(p, "ggplot")
  })
})

test_that("plot_bw_bins_scatter with bg files passes on parameters", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_bins_scatter(bw1, bw2, bg_x = bg_bw, bg_y = bg_bw)
  })

  expect_call(m, 1,
    bw_bins(
      x,
      bg_bwfiles = bg_x,
      bin_size = bin_size,
      genome = genome,
      per_locus_stat = per_locus_stat,
      norm_func = norm_func_x,
      labels = "x"
    )
  )

  expect_args(m, 1,
    x = bw1,
    bg_bwfiles = bg_bw,
    bin_size = 10000,
    genome = "mm9",
    per_locus_stat = "mean",
    norm_func = identity,
    labels = "x"
  )

  expect_args(m, 2,
    x = bw2,
    bg_bwfiles = bg_bw,
    bin_size = 10000,
    genome = "mm9",
    per_locus_stat = "mean",
    norm_func = identity,
    labels = "y"
  )
})

test_that("plot_bw_bins_violin with defaults returns a plot", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_bins_violin(c(bw1, bw2))
    expect_is(p, "ggplot")
  })
})

test_that("plot_bw_bins_violin with bg files passes on parameters", {
  m <- mock(reduced_bins, reduced_bins_2)
  with_mock(bw_bins = m, {
    p <- plot_bw_bins_violin(c(bw1, bw2),
      bg_bwfiles = c(bg_bw, bg_bw),
      labels = c("A", "B"),
      bin_size = 5000,
      norm_func = log2,
      genome = "hg38",
    )
  })

  expect_call(m, 1,
    bw_bins(bwfiles,
      bg_bwfiles = bg_bwfiles,
      labels = labels,
      bin_size = bin_size,
      genome = genome,
      per_locus_stat = per_locus_stat,
      norm_func = norm_func
    )
  )

  expect_args(m, 1,
    bwfiles = c(bw1, bw2),
    bg_bwfiles = c(bg_bw, bg_bw),
    labels = c("A", "B"),
    bin_size = 5000,
    genome = "hg38",
    per_locus_stat = "mean",
    norm_func = log2
  )

})

test_that("plot_bw_bins_violin with highlight returns a plot with jitter", {
    m <- mock(reduced_bins, reduced_bins_2)
    with_mock(bw_bins = m, {
      p <- plot_bw_bins_violin(c(bw1, bw2),
        bg_bwfiles = c(bg_bw, bg_bw),
        labels = c("A", "B"),
        highlight = bed,
        bin_size = 5000,
        norm_func = log2,
        genome = "hg38",
      )
    })

    expect_call(m, 1,
      bw_bins(bwfiles,
        bg_bwfiles = bg_bwfiles,
        labels = labels,
        bin_size = bin_size,
        genome = genome,
        per_locus_stat = per_locus_stat,
        norm_func = norm_func
      )
  )

  expect_args(m, 1,
    bwfiles = c(bw1, bw2),
    bg_bwfiles = c(bg_bw, bg_bw),
    labels = c("A", "B"),
    bin_size = 5000,
    genome = "hg38",
    per_locus_stat = "mean",
    norm_func = log2
  )

})

test_that(
  "plot_bw_bed_summary_heatmap with defaults returns a pheatmap object", {
  m <- mock(summary_values)
  with_mock(bw_bed = m, {
    p <- plot_bw_bed_summary_heatmap(c(bw1, bw2), bedfile = bed)
    expect_is(p, "pheatmap")
  })
})

test_that(
  "plot_bw_bed_summary_heatmap passes parameters on", {
    m <- mock(summary_values)
    with_mock(bw_bed = m, {
      p <- plot_bw_bed_summary_heatmap(c(bw1, bw2),
            bedfile = bed,
            bg_bwfiles <- c(bg_bw, bg_bw),
            aggregate_by = "median",
            norm_func = log2,
            labels = c("bw1", "bw2")
        )
    })

    expect_call(m, 1,
      bw_bed(bwfiles,
        bedfile,
        bg_bwfiles = bg_bwfiles,
        aggregate_by = aggregate_by,
        norm_func = norm_func,
        labels = labels
      )
    )

    expect_args(m, 1,
      bwfiles = c(bw1, bw2),
      bedfile = bed,
      bg_bwfiles = c(bg_bw, bg_bw),
      aggregate_by = "median",
      norm_func = log2,
      labels = c("bw1", "bw2")
    )
})

test_that(
  "plot_bw_profile with defaults returns a ggplot object", {
    m <- mock(profile_values)
    with_mock(bw_profile = m, {
      p <- plot_bw_profile(c(bw1, bw2), bedfile = bed)
      expect_is(p, "ggplot")
    })
})

test_that(
  "plot_bw_profile passes on parameters to bw_profile call", {
    m <- mock(profile_values)
    with_mock(bw_profile = m, {
      p <- plot_bw_profile(c(bw1, bw2),
             bedfile = bed,
             bg_bwfiles = c(bg_bw, bg_bw),
             mode = "start",
             bin_size = 1000,
             upstream = 1500,
             downstream = 1500,
             middle = 1000,
             ignore_strand = TRUE,
             show_error = TRUE,
             norm_func = log2,
             labels = c("bw1", "bw2"))
    })

    expect_call(m, 1,
      bw_profile(bwfiles, bedfile,
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
    )

    expect_args(m, 1,
      c(bw1, bw2),
      bedfile = bed,
      bg_bwfiles = c(bg_bw, bg_bw),
      mode = "start",
      bin_size = 1000,
      upstream = 1500,
      downstream = 1500,
      middle = 1000,
      ignore_strand = TRUE,
      norm_func = log2,
      labels = c("bw1", "bw2")
    )
})
