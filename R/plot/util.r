

colourPalette <- c("darkred", "darkblue", "orange", "green", "cyan", "magenta", "yellow", "red",
                   "#F8766D", "#EF7F49", "#E58700", "#D89000", "#C99800", "#B79F00", "#A3A500",
                   "#8AAB00", "#6BB100", "#39B600", "#00BA38", "#00BD5F", "#00BF7D", "#00C097",
                   "#00C0AF", "#00BFC4", "#00BCD8", "#00B7E9", "#00B0F6", "#00A7FF", "#619CFF",
                   "#9590FF", "#B983FF", "#D376FF", "#E76BF3", "#F564E3", "#FD61D1", "#FF62BC",
                   "#FF67A4", "#FE6E8A")


#' Function used for saving the plots.
#' @param name string.
savePlot <- function(filename, plot, plotSaveFormat = "pdf", plotDpi = 300, cols = 1, rows = 1, plotWidth = 10, plotHeight = 3) {
    fp <- str_replace_all(paste(filename, ".", plotSaveFormat, sep = ""), " ", "")
    ggsave(plot = plot, filename = fp, device = plotSaveFormat, dpi = plotDpi, width = cols * plotWidth, height = rows * plotHeight, limitsize = FALSE)
}


