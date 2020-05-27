## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(elsasserlib)
library(ggplot2)
library(scales)

## ----fig.width=7, fig.height=7------------------------------------------------
normal.plot <- ggplot(iris, aes(x=Sepal.Width, y=Sepal.Length)) + 
  geom_point() + 
  ggtitle("Your usual default ggplot scatterplot")

normal.plot

## ----fig.width=7, fig.height=7------------------------------------------------
normal.plot + 
  ggtitle("Different looks for same plot") +
  theme_elsasserlab_screen()


## ----fig.width=7, fig.height=7------------------------------------------------
normal.plot + 
  ggtitle("Large fonts") +
  theme_elsasserlab_print()


## ----fig.width=7, fig.height=7------------------------------------------------

# Same plot with slightly less huge fonts
normal.plot + 
  ggtitle("Not-so-large fonts") +
  theme_elsasserlab_print(base_size=18)


## ----fig.width=7, fig.height=7------------------------------------------------
normal.plot + 
  ggtitle("Back to screen theme + rotated axis labels") +
  theme_elsasserlab_screen() +
  theme(axis.text.x = element_text(angle=90, hjust=1))


## ----fig.width=7, fig.height=7------------------------------------------------
categorical.plot <- ggplot(iris, aes(x=Sepal.Length, y=Sepal.Width, color=Species)) +
  geom_point() +
  theme_elsasserlab_screen() +
  scale_color_manual(values=palette_categorical(3))

categorical.plot

## ----fig.width=7, fig.height=7------------------------------------------------
show_col(palette_categorical(12))

