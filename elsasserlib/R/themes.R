# Data visualization specifics for our plotting in R
# Plotting themes and palettes.

#' ggplot2-resembling default palette for an arbitrary number of values.
#'
#' @param n Number of values to be generated.
#' @return A list of n colors.
#' @importFrom grDevices hcl
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

#' Palette for qualitative datasets.
#'
#' Returns a palette useful for qualitative datasets.
#' Maximum value of 12 colors.
#'
#' @param n Number of values (must be <=12.
#' @return A list of n colors.
#' @export
palette_categorical <- function(n) {
  max <- 12

  # Names are standard browser-supported names for these.
  full.palette <- c(
    '#F08080', #lightcoral
    '#ADD8E6', #lightblue
    '#32CD32', #limegreen
    '#FFD700', #gold
    '#B22222', #firebrick
    '#4682B4', #steelblue
    '#006400', #darkgreen
    '#FF8C00', #darkorange
    '#800000', #maroon
    '#4B0082', #indigo
    '#008B8B', #darkcyan
    '#D3D3D3' #gray
  )
  if (n > max) {
    warning(paste("Qualitative palette is limited to", max, "colors + gray", sep=" "))
    warning("Returning basic ggplot hue colors.")
    full.palette <- gg_color_hue(n)
  }
  full.palette
}

#' Theme for print versions of plots.
#'
#' A theme where font sizes will generally be pretty large compared to screen
#' values, but good to show in combined panels.
#'
#' @param base_size Base font size. Size of elements in the plot are calculated
#'    from this.
#' @param base_family Font family.
#' @param base_line_size Thickness of lines. Default is defined based on base_size.
#' @importFrom ggplot2 %+replace% theme_minimal theme element_text element_line element_blank rel margin unit
#' @importFrom grDevices rgb
#' @export
theme_elsasserlab_print <- function(base_size = 22,
                                    base_family = "",
                                    base_line_size = base_size / 22) {
  half_line <- base_size / 2

  theme_minimal(base_size = base_size,
                base_family = base_family,
                base_line_size = base_line_size) %+replace%
    theme(
      plot.title = element_text(
        color = rgb(25, 43, 65, maxColorValue = 255),
        face = "bold",
        hjust = 0, size=rel(1.5), margin = margin(b = half_line), vjust=4.5),
      axis.title = element_text(
        color = "black",
        size = rel(1.2), face = "bold"),
      axis.text = element_text(
        color = "black",
        size = rel(0.9)),
      strip.text.x = element_text( color = "black", size= rel(1.2), face="bold"),
      axis.line = element_line(
        colour = "black", size = base_line_size, linetype=1, lineend = "butt"
      ),
      axis.ticks = element_line(
        colour = "black"
      ),
      legend.title = element_text( color = "black", face="bold", size = rel(0.9)),
      legend.text = element_text(color="black",size=rel(0.9)),
      plot.margin=unit(c(2.0,2.0,2.0,2.0),"cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      complete = TRUE
    )
}


#' Theme for screen plotting.
#'
#' A theme for using on screen or single plots that plays nice with the print
#' version.
#'
#' @param base_size Base font size. Size of elements in the plot are calculated
#'    from this.
#' @param base_family Font family.
#' @param base_line_size Thickness of lines. Default is defined based on base_size.
#' @importFrom ggplot2 %+replace% theme_minimal theme element_text element_line element_blank rel margin unit
#' @importFrom grDevices rgb
#' @export
theme_elsasserlab_screen <- function(base_size = 12,
                                     base_family = "",
                                     base_line_size = base_size / 22) {

  half_line <- base_size / 2

  theme_minimal(base_size = base_size,
                base_family = base_family,
                base_line_size = base_line_size) %+replace%
    theme(
      plot.title = element_text(
        color = rgb(25, 43, 65, maxColorValue = 255),
        face = "bold",
        hjust = 0, size=rel(1.2), margin = margin(b = half_line)),
      axis.title = element_text(
        color = rgb(105, 105, 105, maxColorValue = 255),
        size = rel(1)),
      axis.text = element_text(
        color = rgb(105, 105, 105, maxColorValue = 255),
        size = rel(0.75)),
      strip.text.x = element_text( color = rgb(105, 105, 105, maxColorValue = 255), size= rel(1), face="bold"),
      axis.line = element_line(
        colour = rgb(105, 105, 105, maxColorValue = 255), size = base_line_size, linetype=1, lineend = "butt"
      ),
      axis.ticks = element_line(
        colour = "black"
      ),
      legend.title = element_text( color = rgb(105, 105, 105, maxColorValue = 255), size= rel(1.1) ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      complete = TRUE
    )
}
