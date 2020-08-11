context("Tests functions for bigWig plots")

bw1 <- system.file("extdata", "ES_H33_00h_rep1_hoxc.bw", package="elsasserlib")
bw2 <- system.file("extdata", "ES_H33_03h_rep1_hoxc.bw", package="elsasserlib")
bg_bw <- system.file("extdata", "ES_H33_inp_rep1_hoxc.bw", package="elsasserlib")
bed <- system.file("extdata", "chromhmm_hoxc.bed", package="elsasserlib")

test_that("Setup files exist", {
  expect_true(file_test('-f', bw1))
  expect_true(file_test('-f', bw2))
  expect_true(file_test('-f', bg_bw))
  expect_true(file_test('-f', bed))
})

