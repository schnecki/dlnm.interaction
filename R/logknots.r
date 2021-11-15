### Reeimplemented function from original dlnm code

#' Logknots
#'
#' @param lag  numeric/vector  maximum lag, lag range or vector of lags
#' @param nk   numeric         number of knots (must be >0)
#' @return     vector
#' @export
logknots <- function(lag, nk) {
    lag <- as.vector(lag)
    if (length(lag) < 1) stop("Expected a lag or vector for 'lag' in lagknots")

    ## If length of x 1 or 2, interpreted as a lag range, otherwise take the range
    range <- if (length(lag) < 2 && lag < 0) c(lag, 0)
             else if (length(lag) < 2) c(0, lag)
             else range(lag, na.rm = TRUE)
    if (diff(range) <= 0) stop("range must be >0")

    ## Define knots at equally-spaced log-values along lag
    if (nk < 1) stop("number of knots 'nk' must be >= 1")
    knots <- range[1] + exp(((1 + log(diff(range))) / (nk + 1)) * seq(nk) - 1)
    return(knots)
}
