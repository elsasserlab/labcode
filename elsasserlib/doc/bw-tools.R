## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, message=FALSE-----------------------------------------------------
library(elsasserlib)
library(ggplot2)
library(rtracklayer)

## -----------------------------------------------------------------------------
mm9.genome.bins.10000 <- build_bins(bsize=10000, genome="mm9")

head(mm9.genome.bins.10000)

## -----------------------------------------------------------------------------
# Example bw file included in the package.
# Note: this should be a path for an actual file.
bw.file <- system.file("extdata", "example.1000bs.bw", package = "elsasserlib")

# Now we bin it
binned.score <- bw_bins(bw.file, bsize=10000)

head(binned.score)

## ---- warning=FALSE-----------------------------------------------------------
# Convert to dataframe
df <- as.data.frame(binned.score)

ggplot(df, aes(y=score)) + 
  geom_boxplot() + 
  ylim(0,20) +
  theme_elsasserlab_screen()

## -----------------------------------------------------------------------------
# Example bw file included in the package.
# Note: this should be a path for an actual file.
bw.file <- system.file("extdata", "example.1000bs.bw", package = "elsasserlib")
bed.file <- system.file("extdata", "ChromHMM7_subsample.bed", package = "elsasserlib")

# We intersect these two files
summary.score <- bw_bed(bw.file, bed.file, per.locus.stat = "mean")

# We get a GRanges object which one score per locus
head(summary.score)

## -----------------------------------------------------------------------------
# This import is `rtracklayer` import
loci <- import(bed.file)

table(as.factor(loci$name))

## -----------------------------------------------------------------------------
# We intersect these two files
summary.score <- bw_bed(bw.file, bed.file, per.locus.stat = "mean", aggregate.by = mean)

# We get a GRanges object which one score per locus
summary.score

## -----------------------------------------------------------------------------
# We intersect these two files, but this time we summarize by median
median.score <- bw_bed(bw.file,
                       bed.file,
                       per.locus.stat="mean",
                       aggregate.by=median)

# We get a GRanges object which one score group of loci
median.score

## -----------------------------------------------------------------------------
ggplot(summary.score, aes(x=name, y=score, fill=name)) + 
  geom_bar(stat='identity') + 
  ylab('Mean coverage') +
  xlab('ChromHMM7 (subsampled)') +
  theme_elsasserlab_screen()

## ----results = "asis"---------------------------------------------------------
print(?bw_bins)

