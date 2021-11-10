#' Crossbasis class
#'
#' This class defines a Crossbasis
#' @export
Crossbasis <- R6::R6Class(
    classname = "Crossbasis",

    ## Properties
    private = list(
        .x = NULL,        # matrix
        .input = NULL,    # matrix (input vector/matrix)
        .basisvar = NULL, # Basis
        .basislag = NULL, # Basis
        .lag = NULL       # numeric (2)
    ),

    ## Methods
    public = list(
        initialize = function(x, basisvar, basislag) {
            self$input <- as.matrix(x)
            self$basisvar <- basisvar$clone() # Clone, as the object might change
            self$basislag <- basislag$clone() # Clone
            self$lag <- c(min(self$basislag$input), max(self$basislag$input))
            self$x <- computeCrossbasis(self) # Compute crossbasis

        }
    ),

    ## Accessable properties. Active bindings look like fields, but each time they are accessed,
    ## they call a function. They are always publicly visible.
    active = list(
        x = function(value) {
            if (missing(value)) return(private$.x)
            if (!(is.matrix(value)))
                stop("ERROR: Unallowed property value 'x' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.x <- value
            return(self)
        },
        input = function(value) {
            if (missing(value)) return(private$.input)
            if (!(is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'input' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.input <- value
            return(self)
        },
        basisvar = function(value) {
            if (missing(value)) return(private$.basisvar)
            if (!("Basis" %in% class(value)))
                stop("ERROR: Unallowed property ", value, " for 'basisvar' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.basisvar <- value
            return(self)
        },
        basislag = function(value) {
            if (missing(value)) return(private$.basislag)
            if (!("Basis" %in% class(value)))
                stop("ERROR: Unallowed property ", value, " for 'basislag' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.basislag <- value
            return(self)
        },
        lag = function(value) {
            if (missing(value)) return(private$.lag)
            if (!(is.numeric(value) && length(value) == 2))
                stop("ERROR: Unallowed property ", value, " for 'lag' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.lag <- value
            return(self)
        },
        lags = function() {
            return(as.vector(private$.basislag$input))
        }
    )
)

#' Compute the crossbasis by first lagging the variable (basisvar) and then multiplying by the lag.
computeCrossbasis <- function(cb) {

    crossbasis <- matrix(0, nrow = cb$basisvar$dimension[1], ncol = cb$basisvar$dimension[2] * cb$basislag$dimension[2])

    ## After the first transformation, i.e. on the basisvector, create a lagged variable, so we can
    ## transform it again.
    for (v in seq(length = cb$basisvar$dimension[2])) {
        mat <- as.matrix(mkLaggedMatrix(cb$basisvar$x[, v], cb$lags))

        ## make @s(x_t, \eta) = q_{t\cdot}^T * C * \eta@ spline from lagged basis values.
        for (i in seq(length = cb$basislag$dimension[2])) {
            crossbasis[, cb$basislag$dimension[2] * (v - 1) + i] <- mat %*% cb$basislag$x[, i]
        }
    }
    return(crossbasis)
}

#' Lag a matrix
mkLaggedMatrix <- function(x, lags) {
    mkLag <- function(x, lag) {
        if (lag < 0) stop("all lags in basislag have to be >=0. Saw ", lag)
        nx <- c(rep(NA, lag), x[1:(length(x) - lag)])
        return(nx)
    }
    Q <- matrix(0, nrow = length(x), ncol = length(lags))
    i <- 1
    for (l in lags) {
        Q[, i] <- mkLag(x, l)
        i <- i + 1
    }
    return(Q)
}
