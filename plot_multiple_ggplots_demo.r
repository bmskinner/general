# Example script creating a collection of Violin plots

library(ggplot2)
library(gridExtra)
library(grid)
library(cowplot)

dfile = "plot_multiple_ggplots_demo.csv"
data  = read.csv(dfile, header=T);

#' @title Make a plot
#' @param yvar the column to plot
#' @param yaxt the y axis title
#' @return a ggplot
makePlot = function(yvar, yaxt){
  ggplot(data, aes(x=Dataset, y=data[[yvar]])) +
    geom_violin(trim=TRUE, lwd=0.9) +
    geom_boxplot(width=0.1, lwd=0.9, outlier.size=0) +
    ylab( yaxt )+
    theme_classic() +
    theme( legend.position="none", axis.title.x = element_blank())
}

# The columns to plot
vars      = c("Area_square_microns", "Variability")
# vars      = c("Area_square_microns", "Variability", "Perimeter_microns")

# The y-axis labels, using expression() to produce Greek letters and superscript
yaxts     = list( expression( paste("Area (", mu, m^2, ")") ), "Variability")

# mapply allows vectors of arguments to be passed to
# the second parameter onwards of functions.
# Setting SIMPLIFY to false prevents mapply coercing the resulting
# gg to a matrix
plot.list = mapply(makePlot, vars, yaxts, SIMPLIFY = F)

# Construct the plot command
plotCmd   = list(plots=plot.list, numRows= length(vars), numCols=1)

# Draw the plots
# Using cowplot::plot_grid to auto align chart axes:
# https://cran.r-project.org/web/packages/cowplot/vignettes/plot_grid.html
do.call(plot_grid, c(plotCmd$plots, align='v', nrow = plotCmd$numRows, ncol = plotCmd$numCols))
