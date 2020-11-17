context("Test functions for bigWig stats")
library(GenomicRanges)
library(testthat)
library(mockery)

get_file_path <- function(filename) {
  system.file("extdata", filename, package = "elsasserlib")
}

bw1 <- get_file_path("sample_H33_ChIP.bw")
bw2 <- get_file_path("sample_H3K9me3_ChIP.bw")
bg_bw <- get_file_path("sample_Input.bw")
bed <- get_file_path("sample_chromhmm.bed")

bw_limits <- GRanges(seqnames = c("chr15"),
                     ranges = IRanges(c(102723600, 102959000)))

reduced_bins <- bw_bins(c(bw1, bw2), selection = bw_limits)
reduced_bg_bins <- bw_bins(c(bg_bw), selection = bw_limits)

test_that("Setup files exist", {
  expect_true(file_test("-f", bw1))
  expect_true(file_test("-f", bw2))
  expect_true(file_test("-f", bg_bw))
  expect_true(file_test("-f", bed))
})

test_that("bw_bins_diff_analysis passes on parameters", {
  m_func <- mock()
  m_bins <- mock(reduced_bins, reduced_bg_bins)
  with_mock(
    bw_granges_diff_analysis = m_func,
    bw_bins = m_bins,
    bw_bins_diff_analysis(c(bw1, bw2), bg_bw, "treated", "untreated")
  )
  expect_call(m_func, 1,
              bw_granges_diff_analysis(bins_c1,
                                       bins_c2,
                                       label_c1,
                                       label_c2,
                                       estimate_size_factors = estimate_size_factors,
                                       as_granges = as_granges)
  )

  expect_call(m_bins, 1,
              bw_bins(bwfiles_c1, genome = genome, bin_size = bin_size)
  )

  expect_call(m_bins, 2,
              bw_bins(bwfiles_c2, genome = genome, bin_size = bin_size)
  )

  expect_args(m_func, 1,
              granges_c1 = reduced_bins,
              granges_c2 = reduced_bg_bins,
              label_c1 = "treated",
              label_c2 = "untreated",
              estimate_size_factors = FALSE,
              as_granges = FALSE)

  expect_args(m_bins, 1,
              c(bw1, bw2),
              genome = "mm9",
              bin_size = 10000)

  expect_args(m_bins, 2,
              bg_bw,
              genome = "mm9",
              bin_size = 10000)
})


test_that("bw_bed_diff_analysis passes on parameters", {
  m_func <- mock()
  m_bed <- mock(reduced_bins, reduced_bg_bins)
  with_mock(
    bw_granges_diff_analysis = m_func,
    bw_bed = m_bed,
    bw_bed_diff_analysis(c(bw1, bw2), bg_bw, bed, "treated", "untreated", as_granges = TRUE)
  )
  expect_call(m_func, 1,
              bw_granges_diff_analysis(loci_c1,
                                       loci_c2,
                                       label_c1,
                                       label_c2,
                                       estimate_size_factors = estimate_size_factors,
                                       as_granges = as_granges)
  )

  expect_call(m_bed, 1,
              bw_bed(bwfiles_c1, bedfile = bedfile)
  )

  expect_call(m_bed, 2,
              bw_bed(bwfiles_c2, bedfile = bedfile)
  )

  expect_args(m_func, 1,
              granges_c1 = reduced_bins,
              granges_c2 = reduced_bg_bins,
              label_c1 = "treated",
              label_c2 = "untreated",
              estimate_size_factors = FALSE,
              as_granges = TRUE)

  expect_args(m_bed, 1,
              c(bw1, bw2),
              bedfile = bed)

  expect_args(m_bed, 2,
              bg_bw,
              bedfile = bed)
})

