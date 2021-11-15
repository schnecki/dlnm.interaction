
#' Linear Basis
#'
#' @export
Lin <- R6::R6Class(
    classname = "Lin",

    ## Inheritage
    inherit = Basis,

    ## Additional Properties
    public = list(
        initialize = function(x, intercept = FALSE) {
            ns <- names(x)
            x <- as.vector(x)
            basis <- as.matrix(x)
            if (intercept) basis <- cbind(1, basis)
            dimnames(basis) <- list(ns, seq(ncol(basis)))
            className <- get(class(self)[[1]], -1)$classname
            super$initialize(basis, x, className, intercept)
        },
        mkNewWith = function(x) {
            return(Lin$new(x, self$intercept))
        }
    ),

    private = list()
)
