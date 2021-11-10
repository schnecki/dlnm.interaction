
#' Linear Basis
#'
#' @export
Ns <- R6::R6Class(
    classname = "Ns",

    ## Inheritage
    inherit = Basis,

    private = list(
        .degree = NULL,        # numeric (1)
        .knots = NULL,         # vector
        .boundary.knots = NULL # numeric (2)
    ),

    ## Additional Properties
    public = list(
        initialize = function(x, knots, intercept = FALSE, Boundary.knots = range(x)) {
            ns <- names(x)
            x <- as.vector(x)
            basis <- splines::ns(x, knots = knots, intercept = intercept, Boundary.knots = Boundary.knots)
            dimnames(basis) <- list(ns, seq(ncol(basis)))
            className <- get(class(self)[[1]], -1)$classname
            super$initialize(basis, x, className, intercept)
            self$degree <- attr(basis, "degree")
            self$knots <- attr(basis, "knots")
            self$boundary.knots <- attr(basis, "Boundary.knots")
        }
    ),

    active = list(
        degree = function(value) {
            if (missing(value)) return(private$.degree)
            if (!(is.numeric(value) && length(value) == 1))
                stop("ERROR: Unallowed property ", value, " for 'degree' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.degree <- value
            return(self)
        },

        knots = function(value) {
            if (missing(value)) return(private$.knots)
            if (!(is.numeric(value)))
                stop("ERROR: Unallowed property ", value, " for 'knots' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.knots <- value
            return(self)
        },

        boundary.knots = function(value) {
            if (missing(value)) return(private$.boundary.knots)
            if (!(is.numeric(value) && length(value) == 2))
                stop("ERROR: Unallowed property ", value, " for 'boundary.knots' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.boundary.knots <- value
            return(self)
        }
    )
)
