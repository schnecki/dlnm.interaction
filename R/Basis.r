#' Basis class
#'
#' This class defines a basis
#' @export Basis
#' @exportClass Basis
Basis <- R6::R6Class(
    classname = "Basis",

    ## Properties
    private = list(
        .x = NULL,         # matrix
        .input = NULL,     # matrix (input vector/matrix to Basis)
        .fun = NULL,       # character
        .intercept = NULL, # logical
        .cen = NULL,       # numeric
        .range = NULL,     # numeric
        .dimension = NULL  # numeric (2)
    ),

    ## Methods
    public = list(
        initialize = function(x, input, fun, intercept = FALSE, cen = NULL) {
            self$x <- x
            self$input <- as.matrix(input)
            self$fun <- fun
            self$intercept <- intercept
            self$cen <- cen
        },
        #' Make a new basis with input vector with same settings. Must return matrix object.
        mkNewWith = function(input) {
            stop("mkNewWith must be overwritten by each class. Not done for class '", class(self)[1], "'")
        }

    ),

    ## Accessable properties. Active bindings look like fields, but each time they are accessed,
    ## they call a function. They are always publicly visible.
    active = list(

        x = function(value) {
            if (missing(value)) return(private$.x)
            if (!(is.matrix(value)))
                stop("ERROR: Unallowed property value ", value, " for 'x' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.x <- value
            private$.range <- range(value, na.rm = TRUE)
            private$.dimension <- c(nrow(value), ncol(value))
            return(self)
        },

        input = function(value) {
            if (missing(value)) return(private$.input)
            if (!(is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'input' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.input <- value
            return(self)
        },

        fun = function(value) {
            if (missing(value)) return(private$.fun)
            if (!(is.character(value)))
                stop("ERROR: Unallowed property ", value, " for 'fun' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.fun <- value
            return(self)
        },

        intercept = function(value) {
            if (missing(value)) return(private$.intercept)
            if (!(is.logical(value)))
                stop("ERROR: Unallowed property ", value, " for 'intercept' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.intercept <- value
            return(self)
        },

        cen = function(value) {
            if (missing(value)) return(private$.cen)
            if (!(is.null(value) || is.numeric(value) && length(value) == 1))
                stop("ERROR: Unallowed property ", value, " for 'cen' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.cen <- value
            return(self)
        },
        range = function() return(private$.range),        # Range is automatically set
        dimension = function() return(private$.dimension) # Dimension is automatically set
    )

)


## We could implement the + and other operators
`+.Basis` <- function(b1, b2) {
    stop("The function '+' is on purpose not implemented for objects of type Basis. You can access the matrices with basis$x.")
    ret <- b1$clone()
    ret$x <- b1$x + b2$x
    ret$fun <- paste(b1$fun, "+", b2$fun)
    ret$cen <- NULL
    ret$intercept <- b1$intercept || b2$intercept
    return(ret)
}
