

#' Compute GGPlot object for the Relative Risk (RR) dependend on the Lag of a Crosspred object.
#' @param cp          Crosspred object to plot.
#' @param filepath    Filepath to save the plot.
#' @param title       Title of plot.
#' @return object     GGPlot object. Execute, i.e. do not save it in a variable, to plot.
plotLagAt <- function(cp, filepath = NULL, title = NULL, xDist = 0.10, colours = colourPalette, legendPos = c(0.5, 0.10)) {

    lags <- cp$crossbasis$basislag$input
    dists <- purrr::map(seq(length(cp$at)), function(x) (x - 0.5 - length(cp$at) / 2) * xDist)
    minDist <- min(unlist(dists))
    maxDist <- max(unlist(dists))
    xVals <- unlist(purrr::map(dists, function(d) lags + d))
    ids <- unlist(purrr::map(seq(length(cp$at)), function(x) rep(colours[x], length(lags))))
    data <- data.frame(row.names = seq(length(cp$matRRfit)), lag = xVals, RR = c(cp$matRRfit), LowRR = c(cp$matRRlow), HighRR = c(cp$matRRhigh), idE = ids)
    ylab <- pretty(c(data$LowRR, data$HighRR))
    if (is.null(title) && length(cp$at) == 1) title <- paste("Relative Risk for Lag=", round(cp$at, 4), " of ", cp$crossbasisName, sep = "")
    if (is.null(title) && length(cp$at) > 1) title <- paste("Relative Risk for Lags of ", cp$crossbasisName, sep = "")
    legendVals <- unlist(purrr::map(cp$at, function(x) paste("Lag=", round(x, 4))))

    plot <-
        data %>%
        ggplot(aes(lag, RR, ymin = LowRR, ymax = HighRR), group = idE) +
        geom_hline(yintercept = 1, size = 0.5, colour = "darkgray") +
        geom_point(shape = 21, fill = "white", size = 2, aes(colour = idE), show.legend = FALSE) +
        geom_linerange(aes(x = lag, y = RR, ymin = LowRR, ymax = HighRR, colour = idE), size = .8, show.legend = TRUE) +
        scale_x_continuous(breaks = seq(min(data$lag - minDist), max(data$lag) - maxDist, 1)) +
        scale_y_continuous(breaks = ylab) +
        labs(title = title, x = "Lag (days)", y = "RR", color = "Legend") +
        theme(legend.position = legendPos, legend.direction = "horizontal") +
        scale_color_manual(values = colours, labels = legendVals)
    return(plot)
}
