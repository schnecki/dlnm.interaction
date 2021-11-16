
#' Retrieve GGPlot object for the Relative Risk (RR) dependend on the Exposure of a Crosspred object.
#' @param cp          Crosspred object to plot.
#' @param filepath    Filepath to save the plot.
#' @param title       Title of plot.
#' @return object     GGPlot object. Execute, i.e. do not save it in a variable, to plot.
plotExposure <- function(cp, filepath = NULL, title = NULL, colours = colourPalette, legendPos = c(0.5, 0.10)) {

    ## Relative Risk total
    RROverallInteraction <- data.frame(row.names = seq(length(cp$at)), tmean = cp$at, RR = cp$allRRfit, LowRR = cp$allRRlow, HighRR = cp$allRRhigh, idE = colours[1])
    if (!is.null(cp$allRRfitExposure)) {
        RROverallExposure <- data.frame(row.names = seq(length(cp$at)), tmean = cp$at, RR = cp$allRRfitExposure, LowRR = cp$allRRlowExposure, HighRR = cp$allRRhighExposure, idE = colours[2])
        RROverall <- rbind(RROverallInteraction, RROverallExposure)
    } else {
        RROverall <- RROverallInteraction
    }

    xlab <- pretty(RROverall$tmean)
    ylab <- pretty(c(RROverall$LowRR, RROverall$HighRR))
    if (is.null(title)) title <- paste("Relative Risk for Exposure of", cp$crossbasisName)
    legendVals <- base::list(paste("Interaction", cp$crossbasisInteractionName), paste("Exposure", cp$crossbasisExposureName))
    plot <-
        RROverall %>%
        ## ggplot(aes(tmean, RR), group = idE) +
        ggplot(aes(tmean, RR, ymin = LowRR, ymax = HighRR), group = idE) +
        geom_hline(yintercept = 1, size = 0.5, colour = "darkgray") +
        geom_ribbon(aes(ymin = LowRR, ymax = HighRR, colour = idE), fill = "grey80", alpha = 0.4) +
        geom_line(aes(colour = idE), size = 1) +
        # geom_point(aes(mmt, 1), shape = 21, fill = "white", size = 2, colour = "#cb181d", show.legend = FALSE) +
        scale_x_continuous(breaks = xlab) +
        scale_y_continuous(breaks = ylab) +
        ## theme_bw() +
        ## theme(panel.grid.minor = element_blank()) +
        ## labs(x = "Exposure", y = "RR", title = title)
        labs(title = title, x = "Exposure", y = "RR", color = "Legend") +
        theme(legend.position = legendPos, legend.direction = "horizontal") +
        scale_color_manual(values = colours, labels = legendVals)

    return(plot)
}
