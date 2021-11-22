
#' Linear Basis
#'
#' @export
Integer <- R6::R6Class(
    classname = "Integer",

    ## Inheritage
    inherit = Basis,

    ## Private fields
    private = list(
        .values = NULL # vector
    ),

    ## Additional Properties
    public = list(
        initialize = function(x, values = NULL, intercept = FALSE) {
            ns <- names(x)
            x <- as.vector(x)
            self$values <- values
            levels <- if (!is.null(values)) sort(values) else sort(unique(x))
            xfac <- factor(x, levels = levels)
            basis <- as.matrix(outer(xfac, levels, function(a, b) as.integer(a == b)))
            ## If intercept is not required, drop first column
            if (ncol(basis) > 1 && !intercept) basis <- basis[, -1, drop = FALSE]
            else if (ncol(basis) <= 1) intercept <- TRUE

            dimnames(basis) <- list(ns, seq(ncol(basis)))
            className <- get(class(self)[[1]], -1)$classname
            super$initialize(basis, x, className, intercept)
        },
        mkNewWith = function(x) {
            return(Integer$new(x, self$values, self$intercept))
        }
    ),

    active = list(
        values = function(value) {
            if (missing(value)) return(private$.values)
            if (!(base::is.vector(value) || is.null(value)))
                stop("ERROR: Unallowed property ", value, " for 'values' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.values <- value
            return(self)
        }

    )
)
