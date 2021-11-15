
#' Retrieve GGPlot object for the Relative Risk (RR) dependend on the Exposure of a Crosspred object.
#' @param cp          Crosspred object to plot.
#' @param filepath    Filepath to save the plot.
#' @param title       Title of plot.
#' @return object     GGPlot object. Execute, i.e. do not save it in a variable, to plot.
plotExposure <- function(cp, filepath = NULL, title = NULL) {

    ## Relative Risk total
    RROverall <- data.frame(row.names = seq(length(cp$at)), tmean = cp$at, RR = cp$allRRfit, LowRR = cp$allRRlow, HighRR = cp$allRRhigh)

    xlab <- pretty(RROverall$tmean)
    ylab <- pretty(c(RROverall$LowRR, RROverall$HighRR))
    if (is.null(title)) title <- paste("Relative Risk for Exposure of", cp$crossbasisName)
    plot <-
        ggplot(RROverall, aes(tmean, RR)) +
        geom_hline(yintercept = 1, size = 0.5, colour = "darkgray") +
        geom_ribbon(aes(ymin = LowRR, ymax = HighRR), fill = "grey80", alpha = 0.7) +
        geom_line(colour = "darkred", size = 1) +
        # geom_point(aes(mmt, 1), shape = 21, fill = "white", size = 2, colour = "#cb181d", show.legend = FALSE) +
        scale_x_continuous(breaks = xlab) +
        scale_y_continuous(breaks = ylab) +
        theme_bw() +
        theme(panel.grid.minor = element_blank()) +
        labs(x = "Exposure", y = "RR", title = title)
    return(plot)
}
