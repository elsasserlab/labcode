context("Test functions for bigWig handling")
library(GenomicRanges)
library(rtracklayer)


toy_example <- function(bw1,
                        bw2,
                        bw_special,
                        bed_with_names,
                        bg1,
                        bg2) {
  granges <- GRanges(
    seqnames = Rle(c("chr1", "chr2"), c(10, 10)),
    ranges = IRanges(c(seq(1, 181, by = 20), seq(1, 181, by = 20)),
                     c(seq(20, 200, by = 20), seq(20, 200, by = 20))),
    score = 1:20
  )

  gr2 <- GRanges(
    seqnames = Rle(c("chr1", "chr2"), c(10, 10)),
    ranges = IRanges(c(seq(1, 181, by = 20), seq(1, 181, by = 20)),
                     c(seq(20, 200, by = 20), seq(20, 200, by = 20))),
    score = 20:1
  )

  gr_bg1 <- GRanges(
    seqnames = Rle(c("chr1", "chr2"), c(10, 10)),
    ranges = IRanges(c(seq(1, 181, by = 20), seq(1, 181, by = 20)),
                     c(seq(20, 200, by = 20), seq(20, 200, by = 20))),
    score = rep(1, 20)
  )

  gr_bg2 <- GRanges(
    seqnames = Rle(c("chr1", "chr2"), c(10, 10)),
    ranges = IRanges(c(seq(1, 181, by = 20), seq(1, 181, by = 20)),
                     c(seq(20, 200, by = 20), seq(20, 200, by = 20))),
    score = rep(2, 20)
  )

  labeled_gr <- GRanges(
    seqnames = c("chr1", "chr1", "chr2", "chr2", "chr2"),
    ranges = IRanges(c(21,  61, 21, 111, 161),
                     c(40, 100, 40, 130, 180))
  )

  labeled_gr$name <- c("typeA", "typeB", "typeA", "typeB", "typeB")

  chromsizes <- c(200, 200)

  seqlengths(granges) <- chromsizes
  seqlengths(gr2) <- chromsizes
  seqlengths(labeled_gr) <- chromsizes
  seqlengths(gr_bg1) <- chromsizes
  seqlengths(gr_bg2) <- chromsizes

  export(granges, bw1)
  export(gr2, bw2)
  export(gr2, bw_special)
  export(labeled_gr, bed_with_names)
  export(gr_bg1, bg1)
  export(gr_bg2, bg2)
}

bw1 <- tempfile("bigwig", fileext = ".bw")
bw2 <- tempfile("bigwig", fileext = ".bw")
bg1 <- tempfile("bigwig_bg1", fileext = ".bw")
bg2 <- tempfile("bigwig_bg2", fileext = ".bw")
bw_special <- tempfile("bigwig", fileext = ".bw")
bed_with_names <- tempfile("bed", fileext = ".bed")

tiles <- tileGenome(c(chr1 = 200, chr2 = 200),
  tilewidth = 20,
  cut.last.tile.in.chrom = TRUE
)


setup(toy_example(bw1, bw2, bw_special, bed_with_names, bg1, bg2))

teardown({
  unlink(c(bw1, bw2, bg1, bg2, bw_special, bed_with_names))
})

test_that("Setup files exist", {
  expect_true(file_test("-f", bw1))
  expect_true(file_test("-f", bw2))
  expect_true(file_test("-f", bg1))
  expect_true(file_test("-f", bg2))
  expect_true(file_test("-f", bed_with_names))
  expect_true(file_test("-f", bw_special))
})

test_that("bw_ranges returns a GRanges object", {
  bins <- bw_ranges(bw1, tiles, per_locus_stat = "mean")
  expect_is(bins, "GRanges")

})

test_that("bw_ranges returns correct values", {
  bins <- bw_ranges(bw1, tiles, per_locus_stat = "mean")

  expect_equal(bins[1]$score, 1)
  expect_equal(bins[2]$score, 2)
  expect_equal(bins[10]$score, 10)
})

test_that("bw_ranges returns correct values on subset", {
  subset <- GRanges(seqnames = c("chr1"),
                    ranges = IRanges(c(10, 40)))

  bins <- bw_ranges(bw1, tiles, per_locus_stat = "mean", selection = subset)

  expect_equal(bins[1]$score, 1)
  expect_equal(bins[2]$score, 2)
  expect_equal(length(bins), 2)
})

test_that("multi_bw_ranges returns correct values", {
  values <- multi_bw_ranges(c(bw1, bw2), c("bw1", "bw2"), tiles)

  expect_is(values, "GRanges")
  expect_equal(values[1]$bw1, 1)
  expect_equal(values[1]$bw2, 20)
  expect_equal(values[2]$bw1, 2)
  expect_equal(values[2]$bw2, 19)
})

test_that(
  "multi_bw_ranges_norm with bwfiles == background returns all 1 values", {
    values <- multi_bw_ranges_norm(c(bw1, bw2),
                bg_bwfilelist = c(bw1, bw2),
                c("bw1", "bw2"),
                tiles,
                norm_func = identity
              )

  expect_is(values, "GRanges")
  expect_equal(values[1]$bw1, 1)
  expect_equal(values[1]$bw2, 1)
  expect_equal(values[2]$bw1, 1)
  expect_equal(values[2]$bw2, 1)
})

test_that("multi_bw_ranges_norm returns correct values", {
  values <- multi_bw_ranges_norm(c(bw1, bw2),
              bg_bwfilelist = c(bg1, bg2),
              c("bw1", "bw2"),
              tiles,
              norm_func = identity
            )

  expect_is(values, "GRanges")
  expect_equal(values[1]$bw1, 1)
  expect_equal(values[2]$bw1, 2)
  expect_equal(values[1]$bw2, 10)
  expect_equal(values[2]$bw2, 9.5)
})

test_that(
  "multi_bw_ranges_norm fails on background length not matching bwlist", {
    expect_error({
      values <- multi_bw_ranges_norm(c(bw1, bw2),
                  bg_bwfilelist = c(bg1),
                  c("bw1", "bw2"),
                  tiles,
                  norm_func = identity
                )
    },
    "Background and signal bwfile lists must have the same length.")
})

test_that("multi_bw_ranges returns correct values for single bigWig", {
  values <- multi_bw_ranges(bw1, "bw1", tiles)

  expect_is(values, "GRanges")
  expect_equal(values[1]$bw1, 1)
  expect_equal(values[2]$bw1, 2)
})

test_that("multi_bw_ranges returns correct values on subset", {
  subset <- GRanges(seqnames = c("chr1"), ranges = IRanges(c(30, 50)))

  values <- multi_bw_ranges(c(bw1, bw2),
              c("bw1", "bw2"),
              tiles,
              selection = subset
            )

  expect_is(values, "GRanges")
  expect_equal(values[1]$bw1, 2)
  expect_equal(values[1]$bw2, 19)
  expect_equal(values[2]$bw1, 3)
  expect_equal(values[2]$bw2, 18)
})

test_that("bw_bed returns correct per locus values", {
  values <- bw_bed(bw1, bed_with_names, labels = "bw1", per_locus_stat = "mean")
  expect_is(values, "GRanges")
  expect_equal(values[1]$bw1, 2)
  expect_equal(values[2]$bw1, 4.5)
})

test_that("bw_bins returns correct per locus values", {
  values <- bw_bins(bw1,
              selection = import(bed_with_names),
              labels = "bw1",
              per_locus_stat = "mean"
            )

  expect_is(values, "GRanges")
  expect_equal(values[1]$bw1, 5.5)
  expect_equal(values[2]$bw1, 15.5)
})

test_that("bw_bins returns correct per locus values, with bg", {
  values <- bw_bins(bw1,
              bg_bwfiles = bw1,
              selection = import(bed_with_names),
              labels = "bw1",
              per_locus_stat = "mean"
            )

  # FIXME: underlying normalized values in GRanges object are matrix
  # expect_is(values, 'GRanges')
  # expect_equal(values[1]$bw1, 1)
  # expect_equal(values[2]$bw1, 1)
})

test_that("bw_bed returns correct per locus values on multiple files", {
  values <- bw_bed(c(bw1, bw2),
              bed_with_names,
              labels = c("bw1", "bw2"),
              per_locus_stat = "mean"
            )

  expect_is(values, "GRanges")
  expect_equal(values[1]$bw1, 2)
  expect_equal(values[2]$bw1, 4.5)
  expect_equal(values[1]$bw2, 19)
  expect_equal(values[2]$bw2, 16.5)
})

test_that(
  "bw_bed returns correct per locus values on multiple files, with bg", {
  values <- bw_bed(c(bw1, bw2),
              bg_bwfiles = c(bg1, bg2),
              bed_with_names,
              labels = c("bw1", "bw2"),
              per_locus_stat = "mean"
            )

  expect_is(values, "GRanges")
  expect_equal(values[1]$bw2, 9.5)
  expect_equal(values[2]$bw2, 8.25)
})


test_that("bw_bed handles default names with special characters", {
  values <- bw_bed(bw_special,
              bed_with_names,
              aggregate_by = "true_mean"
            )

  expect_is(values, "data.frame")
})

test_that("bw_bed crashes on wrong number of labels for multiple files", {
  expect_error({
    values <- bw_bed(c(bw1, bw2), bed_with_names,
                labels = "bw1",
                per_locus_stat = "mean"
              )
  },
  "BigWig file list and column names must have the same length.")

})

test_that("bw_bed returns correct mean-of-means aggregated values", {
  values <- bw_bed(bw1,
              bed_with_names,
              labels = "bw1",
              per_locus_stat = "mean",
              aggregate_by = "mean"
            )

  expect_is(values, "data.frame")
  expect_equal(values["typeA", "bw1"], 7)
  expect_equal(values["typeB", "bw1"], 13.3333333333)
})

test_that("bw_bed returns correct true_mean aggregated values", {
  values <- bw_bed(bw1, bed_with_names,
              labels = "bw1",
              per_locus_stat = "mean",
              aggregate_by = "true_mean"
            )

  expect_is(values, "data.frame")
  expect_equal(values["typeA", "bw1"], 7)
  expect_equal(values["typeB", "bw1"], 11.125)
})

test_that("bw_bed on an empty list throws an error", {
  expect_error({
    values <- bw_bed(c(), bed_with_names,
                per_locus_stat = "mean",
                aggregate_by = "true_mean"
              )
  },
  "File list provided is empty.")

})

test_that("bw_profile on an empty list throws an error", {
  expect_error({
    values <- bw_profile(c(),
                bed_with_names,
                labels = NULL
              )
  },
  "File list provided is empty.")

})

test_that("bw_bed on non-existing bed file throws an error", {
  expect_error({
    values <- bw_bed(bw1,
                "invalidname.bed",
                per_locus_stat = "mean",
                aggregate_by = "true_mean"
              )
  },
  "Files not found: invalidname.bed")
})

test_that("bw_profile on non-existing bed file throws an error", {
  expect_error({
    values <- bw_profile(bw1,
                bedfile = "invalidname.bed",
                labels = NULL
              )
  },
  "Files not found: invalidname.bed")
})


test_that("bw_bed errors on non existing files on bwlist", {
  expect_error({
    values <- bw_bed(c(bw1, "invalidname.bw"),
                bed_with_names,
                per_locus_stat = "mean",
                aggregate_by = "true_mean"
              )
  },
  "Files not found: invalidname.bw")
})


test_that("bw_profile errors on non existing files on bwlist", {
  expect_error({
    values <- bw_profile(c(bw1, "invalidname.bw"),
                bedfile = bed_with_names,
                labels = NULL
              )
  },
  "Files not found: invalidname.bw")
})


test_that("bw_profile runs quiet on valid parameters", {
  expect_silent({
    values <- bw_profile(c(bw1, bw2),
                bedfile = bed_with_names,
                upstream = 1,
                downstream = 1,
                bin_size = 1
              )
  })

})

test_that("bw_profile runs quiet on valid parameters, mode start", {
  expect_silent({
    values <- bw_profile(c(bw1, bw2),
                bedfile = bed_with_names,
                upstream = 1,
                downstream = 1,
                bin_size = 1,
                mode = "start"
              )
  })

})

test_that(
  "bw_profile runs quiet on valid parameters, mode start, with background", {
  expect_silent({
    values <- bw_profile(c(bw1, bw2),
                bg_bwfiles = c(bg1, bg2),
                bedfile = bed_with_names,
                upstream = 1,
                downstream = 1,
                bin_size = 1,
                mode = "start"
              )
  })

})

test_that("bw_profile normalized returns correct values", {
  values <- bw_profile(c(bw1, bw2),
              bg_bwfiles = c(bw1, bw2),
              bedfile = bed_with_names,
              upstream = 1,
              downstream = 1,
              bin_size = 1,
              mode = "start"
            )

  expect_is(values, "data.frame")
  expect_equal(values[1, "mean"], 1)
  expect_equal(values[2, "mean"], 1)

})

test_that("bw_profile fails if labels and bwfiles have different length", {
  expect_error({
    values <- bw_profile(c(bw1, bw2),
                bedfile = bed_with_names,
                labels = c("bw1"),
                upstream = 1,
                downstream = 1,
                bin_size = 1,
                mode = "start"
              )
  },
  "labels and bwfiles must have the same length")
})


test_that("bw_profile runs quiet on valid parameters, mode end", {
  expect_silent({
    values <- bw_profile(c(bw1, bw2),
                bedfile = bed_with_names,
                upstream = 1,
                downstream = 1,
                bin_size = 1,
                mode = "end"
              )
  })

})


test_that("bw_profile runs quiet on valid parameters with background", {
  expect_silent({
    values <- bw_profile(c(bw1, bw2),
                bg_bwfiles = c(bg1, bg2),
                bedfile = bed_with_names,
                upstream = 1,
                downstream = 1,
                bin_size = 1
              )
  })

})

test_that("bw_profile throws error on flanking region smaller than bin size", {
  expect_error({
    values <- bw_profile(c(bw1, bw2),
                bedfile = bed_with_names,
                upstream = 1,
                downstream = 1,
                bin_size = 10
              )
  },
  "bin size must be smaller than flanking regions")

})

test_that("bw_profile throws error on negative bin size", {
  expect_error({
    values <- bw_profile(c(bw1, bw2),
                bedfile = bed_with_names,
                upstream = 1,
                downstream = 1,
                bin_size = -10
              )
  },
  "bin size must be a positive value: -10")

})

test_that("bw_profile throws error on negative upstream value", {
  expect_error({
    values <- bw_profile(c(bw1, bw2),
                bedfile = bed_with_names,
                upstream = -10,
                downstream = 10,
                bin_size = 10
              )
  },
  "upstream size must be a positive value: -10")

})

test_that("bw_profile throws error on negative downstream value", {
  expect_error({
    values <- bw_profile(c(bw1, bw2),
                bedfile = bed_with_names,
                upstream = 10,
                downstream = -10,
                bin_size = 10
              )
  },
  "downstream size must be a positive value: -10")

})

test_that("bw_bed returns correct median-of-means aggregated values", {
  values <- bw_bed(bw1, bed_with_names,
              labels = "bw1",
              per_locus_stat = "mean",
              aggregate_by = "median"
            )

  expect_is(values, "data.frame")
  expect_equal(values["typeA", "bw1"], 7)
  expect_equal(values["typeB", "bw1"], 16.5)
})

test_that("build_bins crashes on unknown or not included genome", {
  expect_error({
    build_bins(bin_size = 10000, genome = "mm10")
  },
  "Supported genomes: mm9, hg38")
})

test_that("build_bins runs for mm9", {
  values <- build_bins(bin_size = 50000, genome = "mm9")
  expect_is(values, "GRanges")
})

test_that("build_bins creates bins of correct size", {
  values <- build_bins(bin_size = 50000, genome = "mm9")
  expect_is(values, "GRanges")
  expect_equal(width(ranges(values[1])), 50000)
})

test_that("build_bins runs for hg38", {
  values <- build_bins(bin_size = 50000, genome = "hg38")
  expect_is(values, "GRanges")
})

test_that("bw_bed throws error on not implemented aggregate_by", {
  expect_error(
    values <- bw_bed(bw1, bed_with_names,
      labels = "bw1",
      per_locus_stat = "mean",
      aggregate_by = "max"
    )
  )
})
