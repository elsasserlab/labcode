context("Tests functions for bw handling")
library(GenomicRanges)
library(rtracklayer)

toy_example <- function(bw1, bw2, bw_special, bed_with_names) {
  gr <- GRanges(
    seqnames = Rle(c("chr1", "chr2"), c(10, 10)),
    ranges = IRanges(c(seq(1,181, by=20),seq(1,181, by=20)),
                     c(seq(20,200, by=20),seq(20,200, by=20))),
    score = 1:20
  )

  gr2 <- GRanges(
    seqnames = Rle(c("chr1", "chr2"), c(10, 10)),
    ranges = IRanges(c(seq(1,181, by=20),seq(1,181, by=20)),
                     c(seq(20,200, by=20),seq(20,200, by=20))),
    score = 20:1
  )

  labeled_gr <- GRanges(
    seqnames = c('chr1', 'chr1', 'chr2', 'chr2', 'chr2'),
    ranges = IRanges(c(21,  61, 21, 111, 161),
                     c(40, 100, 40, 130, 180))
  )

  labeled_gr$name <- c('typeA', 'typeB', 'typeA', 'typeB', 'typeB')

  chromsizes <- c(200,200)

  seqlengths(gr) <- chromsizes
  seqlengths(gr2) <- chromsizes
  seqlengths(labeled_gr) <- chromsizes

  export(gr, bw1)
  export(gr2, bw2)
  export(gr2, bw_special)
  export(labeled_gr, bed_with_names)

}


bw1 <- tempfile('bigwig', fileext='.bw')
bw2 <- tempfile('bigwig', fileext='.bw')
bw_special <- tempfile('bigwig', fileext='.bw')
bed_with_names <-tempfile('bed', fileext='.bed')

setup(toy_example(bw1, bw2, bw_special, bed_with_names))
teardown({
  unlink(bw1)
  unlink(bw2)
  unlink(bw_special)
  unlink(bed_with_names)
})

test_that("Setup files exist", {
  expect_true(file_test('-f', bw1))
  expect_true(file_test('-f', bw2))
  expect_true(file_test('-f', bed_with_names))
})

test_that("bw_ranges returns a GRanges object", {
  tiles <- tileGenome(c(chr1=200,chr2=200),
                      tilewidth=20,
                      cut.last.tile.in.chrom=T)

  bins <- bw_ranges(bw1, tiles, per.locus.stat='mean')
  expect_is(bins, 'GRanges')

})


test_that("bw_ranges returns correct values", {
  tiles <- tileGenome(c(chr1=200,chr2=200),
                      tilewidth=20,
                      cut.last.tile.in.chrom=T)

  bins <- bw_ranges(bw1, tiles, per.locus.stat='mean')
  expect_equal(bins[1]$score, 1)
  expect_equal(bins[2]$score, 2)
  expect_equal(bins[10]$score, 10)
})

test_that("bw_ranges returns correct values on subset", {
  tiles <- tileGenome(c(chr1=200,chr2=200),
                      tilewidth=20,
                      cut.last.tile.in.chrom=T)

  subset <- GRanges(seqnames=c('chr1'),
                    ranges=IRanges(c(10, 40)))

  bins <- bw_ranges(bw1, tiles, per.locus.stat='mean', selection=subset)

  expect_equal(bins[1]$score, 1)
  expect_equal(bins[2]$score, 2)
  expect_equal(length(bins), 2)
})

test_that("multi_bw_ranges returns correct values", {
  tiles <- tileGenome(c(chr1=200,chr2=200),
                      tilewidth=20,
                      cut.last.tile.in.chrom=T)

  values <- multi_bw_ranges(c(bw1, bw2), c('bw1','bw2'), tiles)

  expect_is(values, 'GRanges')
  expect_equal(values[1]$bw1, 1)
  expect_equal(values[1]$bw2, 20)
  expect_equal(values[2]$bw1, 2)
  expect_equal(values[2]$bw2, 19)

})

test_that("multi_bw_ranges returns correct values for single bigWig", {
  tiles <- tileGenome(c(chr1=200,chr2=200),
                      tilewidth=20,
                      cut.last.tile.in.chrom=T)

  values <- multi_bw_ranges(bw1, 'bw1', tiles)

  expect_is(values, 'GRanges')
  expect_equal(values[1]$bw1, 1)
  expect_equal(values[2]$bw1, 2)

})

test_that("multi_bw_ranges returns correct values on subset", {
  tiles <- tileGenome(c(chr1=200,chr2=200),
                      tilewidth=20,
                      cut.last.tile.in.chrom=T)

  subset <- GRanges(seqnames=c('chr1'),
                    ranges=IRanges(c(30, 50)))

  values <- multi_bw_ranges(c(bw1, bw2),
                            c('bw1','bw2'),
                            tiles,
                            selection=subset)

  expect_is(values, 'GRanges')
  expect_equal(values[1]$bw1, 2)
  expect_equal(values[1]$bw2, 19)
  expect_equal(values[2]$bw1, 3)
  expect_equal(values[2]$bw2, 18)

})

test_that("bw_bed returns correct per locus values", {
  values <- bw_bed(bw1,
                   bed_with_names,
                   colnames='bw1',
                   per.locus.stat='mean')

  expect_is(values, 'GRanges')
  expect_equal(values[1]$bw1, 2)
  expect_equal(values[2]$bw1, 4.5)
})

test_that("bw_bed returns correct per locus values on multiple files", {
  values <- bw_bed(c(bw1, bw2),
                   bed_with_names,
                   colnames=c('bw1','bw2'),
                   per.locus.stat='mean')

  expect_is(values, 'GRanges')
  expect_equal(values[1]$bw1, 2)
  expect_equal(values[2]$bw1, 4.5)
  expect_equal(values[1]$bw2, 19)
  expect_equal(values[2]$bw2, 16.5)
})

test_that("bw_bed handles default names with special characters", {
  values <- bw_bed(bw_special,
                   bed_with_names,
                   aggregate.by='true_mean')

  expect_is(values, 'data.frame')
})

test_that("bw_bed crashes on wrong number of colnames for multiple files", {
  expect_error({ values <- bw_bed(c(bw1, bw2),
                   bed_with_names,
                   colnames='bw1',
                   per.locus.stat='mean') },
               "BigWig file list and column names must have the same length.")

})

test_that("bw_bed returns correct mean-of-means aggregated values", {
  values <- bw_bed(bw1,
                   bed_with_names,
                   colnames='bw1',
                   per.locus.stat='mean',
                   aggregate.by='mean')

  expect_is(values, 'data.frame')
  expect_equal(values['typeA', 'bw1'], 7)
  expect_equal(values['typeB', 'bw1'], 13.3333333333)
})

test_that("bw_bed returns correct true_mean aggregated values", {
  values <- bw_bed(bw1,
                   bed_with_names,
                   colnames='bw1',
                   per.locus.stat='mean',
                   aggregate.by='true_mean')

  expect_is(values, 'data.frame')
  expect_equal(values['typeA', 'bw1'], 7)
  expect_equal(values['typeB', 'bw1'], 11.125)
})

test_that("bw_bed on an empty list throws an error", {
  expect_error({ values <- bw_bed(c(),
                   bed_with_names,
                   per.locus.stat='mean',
                   aggregate.by='true_mean')},
               "File list provided is empty."
  )

})

test_that("bw_profile on an empty list throws an error", {
  expect_error({ values <- bw_profile(c(),
                                      bed_with_names,
                                      colnames=NULL)},
               "File list provided is empty."
  )

})

test_that("bw_bed on non-existing bed file throws an error", {
  expect_error({ values <- bw_bed(bw1,
                                  'invalidname.bed',
                                  per.locus.stat='mean',
                                  aggregate.by='true_mean')},
               "Files not found: invalidname.bed"
  )
})

test_that("bw_profile on non-existing bed file throws an error", {
  expect_error({ values <- bw_profile(bw1,
                                  'invalidname.bed',
                                  colnames=NULL)},
               "Files not found: invalidname.bed"
  )
})


test_that("bw_bed errors on non existing files on bwlist", {
  expect_error({ values <- bw_bed(c(bw1, 'invalidname.bw'),
                                  bed_with_names,
                                  per.locus.stat='mean',
                                  aggregate.by='true_mean')},
               "Files not found: invalidname.bw"
  )
})


test_that("bw_profile errors on non existing files on bwlist", {
  expect_error({ values <- bw_profile(c(bw1, 'invalidname.bw'),
                                      bed_with_names,
                                      colnames=NULL)},
               "Files not found: invalidname.bw"
  )
})


test_that("bw_profile runs quiet on valid parameters", {
  expect_silent({values <- bw_profile(c(bw1), bed_with_names, colnames=NULL)})

})


test_that("bw_bed returns correct median-of-means aggregated values", {
  values <- bw_bed(bw1,
                   bed_with_names,
                   colnames='bw1',
                   per.locus.stat='mean',
                   aggregate.by='median')

  expect_is(values, 'data.frame')
  expect_equal(values['typeA', 'bw1'], 7)
  expect_equal(values['typeB', 'bw1'], 16.5)
})


test_that("build_bins crashes on unknown or not included genome", {
  expect_error(build_bins(bsize=10000, genome='mm10'))
})

test_that("bw_bed throws error on not implemented aggregate.by", {
  expect_error(values <- bw_bed(bw1,
                                bed_with_names,
                                colnames='bw1',
                                per.locus.stat='mean',
                                aggregate.by='max'))
})

