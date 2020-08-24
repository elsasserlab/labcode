#' Renders a report with specified bigwig directory and bedfiles
#'
#' @param bwdir Directory where bigwig files are
#' @param bed Bed file that contains some regions of interest to summarize
#' @param highlight Optional bed file to highlight bins overlapping (bin analysis)
#' @param bsize Bin size used.
#' @param minoverlap Minimum overlap between a bin and a highlighted region to be marked.
#' @param genome Genome used (supported: mm9, hg38)
#' @param force Force creation of tables or use them if they exist.
#' @param outdir Output directory
#' @param outfile Output filename (without the full path)
#' @importFrom rmarkdown render
#' @export
render_bw_report <- function(bwdir,
                             bed,
                             outfile='bwreport.html',
                             highlight=NULL,
                             bsize=10000,
                             minoverlap=2000,
                             genome='mm9',
                             outdir='.',
                             force=FALSE) {

    template <- system.file("templates", "bw_report.Rmd", package="elsasserlib")

    render(template,
           params=list(datadir=bwdir,
                       bed=bed,
                       binsize=bsize,
                       highlight=highlight,
                       genome=genome,
                       outdir=outdir,
                       force=force),
           output_file=c(basename(outfile), basename(outfile)),
           output_format=c('html_document', 'md_document')
           )
}

