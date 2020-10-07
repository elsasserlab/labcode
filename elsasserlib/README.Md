![](https://github.com/elsasserlab/labcode/workflows/R-CMD-check/badge.svg)

## Installation

`elsasserlib` is a package under development. You will need `remotes` to
install it (and `devtools` if you plan to work on it):

    install.packages(c('devtools', 'remotes'))
    
Additionally, there is an issue with installing dependencies that come from 
`BioConductor` repository. Therefore, I recommend installing manually these
dependencies before running the installation:

    install.packages(('BiocManager'))
    BiocManager::install(c('GenomeInfoDbData',
        'GenomeInfoDb',
        'GenomicRanges',
        'rtracklayer',
        'BSgenome.Mmusculus.UCSC.mm9',
        'BSgenome.Hsapiens.UCSC.hg38',
        'DESeq2'))

Then you can install directly from this GitHub repository:

    library(remotes)
    install_github('elsasserlab/labcode',
        subdir='elsasserlib',
        build_vignettes=TRUE)

## Getting started

I recommend looking at the vignettes as they can give a comprehensive
overview of what is available in the package. These are the groups of
functions that are included in `elsasserlib` as of now:

- `bwtools`. Functionality to handle `bigWig` files. Importing, binning
    and intersecting, aggregating `BED` files and so on.
- `normGeneCounts`. Normalization functions on gene counts matrices.
- `themes`. Helper `ggplot2` based functions that can make plots all
    share the same looks for publication, including a print-specific
    one that has larger-than-usual font size.
    
## Troubleshooting

**Q**: When running `install_github` I get the following error:

    Error: package or namespace load failed for ‘GenomeInfoDb’ in loadNamespace(i, c(lib.loc, .libPaths()), versionCheck = vI[[i]]):
    there is no package called ‘GenomeInfoDbData’
    Error: package ‘GenomeInfoDb’ could not be loaded
    Execution halted
    
**A**: This seems to be a problem that comes from installing `Bioconductor`
dependencies. A workaround is installing the `BioConductor` packages manually: 

    if (!requireNamespace('BiocManager', quietly = TRUE))
        install.packages('BiocManager')

    BiocManager::install(c('GenomeInfoDbData',
        'GenomeInfoDb',
        'GenomicRanges',
        'rtracklayer',
        'BSgenome.Mmusculus.UCSC.mm9',
        'BSgenome.Hsapiens.UCSC.hg38'))
 

