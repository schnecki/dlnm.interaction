#' Crosspred class
#'
#' This class defines and computes a prediction for Crosspred
#' @export
CrosspredInteraction <- R6::R6Class(
    classname = "CrosspredInteraction",
    inherit = Crosspred,

    ## Properties
    private = base::append(base::list(
        .crossbasisInteraction = NULL,     # Crossbasis
        .crossbasisExposure = NULL,        # Crossbasis
        .crossbasisInteractionName = NULL, # character
        .crossbasisExposureName = NULL,    # character
        .coefInteraction = NULL,           # vector
        .coefExposure = NULL,              # vector
        .vcovInteraction = NULL,           # matrix
        .vcovExposure = NULL,              # matrix
        .vcovIntersection = NULL           # matrix
    ), CrosspredInteraction.privateFunctions),

    ## Methods
    public = list(
        initialize = function(crossbasisExposure, crossbasisInteraction, model, at, coefInteraction = NULL, coefExposure = NULL, vcovInteraction = NULL, vcovExposure = NULL, vcovIntersection = NULL
                            , cen = NULL, ci.level=0.95, cumul = FALSE, model.link=NULL) {
            self$crossbasisInteraction <- crossbasisInteraction$clone()
            self$crossbasisExposure <- crossbasisExposure$clone()
            self$crossbasisInteractionName <- deparse(substitute(crossbasisInteraction))
            self$crossbasisExposureName <- deparse(substitute(crossbasisExposure))
            if (!base::is.null(coefInteraction) && !base::is.null(vcovInteraction) && !base::is.null(coefExposure) && !base::is.null(vcovExposure) && !base::is.null(vcovIntersection)) {
                self$coefInteraction <- coefInteraction
                self$vcovInteraction <- vcovInteraction
                self$coefExposure <- coefExposure
                self$vcovExposure <- vcovExposure
                self$vcovIntersection <- vcovIntersection
            } else if (base::is.null(model)) {
                stop("Crosspred expects either a class of type model (", private$.knownModels, ") or `coef{Interaction,Exposure}`, `vcov{Interaction,Exposure}` and `vcovIntersection` have to be specified")
            } else {
                self$model <- model
                self$model.class <- class(model)
                self$link <- private$getLink(model.link)
                self$coefInteraction <- private$calcCoef(self$crossbasisInteraction, self$crossbasisInteractionName)
                self$vcovInteraction <- private$calcVcov(self$crossbasisInteraction, self$crossbasisInteractionName)
                self$coefExposure <- private$calcCoef(self$crossbasisExposure, self$crossbasisExposureName)
                self$vcovInteraction <- private$calcVcov(self$crossbasisExposure, self$crossbasisExposureName)
                self$vcovIntersection <- private$calcVcov(self$crossbasisExposure, self$crossbasisExposureName)
            }
            self$at <- at
            self$cen <- private$mkcen(cen, crossbasis)
            self$cumul <- cumul
            self$ci.level <- ci.level
            mkPredictions(self) # create all predictionsx
        },
        #' Exposure-RelativeRisk plot. See `plotExposure`
        plotExposure = function(filepath = NULL, title = NULL) plotExposure(self, filepath, title),
        #' Plot lag. See `plotLagAt`.
        plotLagAt = function(at, filepath = NULL, title = NULL, xDist = 0.10, colours = colourPalette, legendPos = c(0.5, 0.10)) {
            if (!base::is.numeric(at)) stop("Expecting numeric for at in function plotLagAt")
            at <- as.vector(at)
            cpAt <- Crosspred$new(self$crossbasis, self$model, at, self$coef, self$vcov, self$cen, self$ci.level, self$cumul, self$model.link)
            cpAt$crossbasisName <- self$crossbasisName
            plotLagAt(cpAt, filepath, title, xDist, colours, legendPos)
        }
    ),

    ## Accessable properties. Active bindings look like fields, but each time they are accessed,
    ## they call a function. They are always publicly visible.
    active = list(
        crossbasisInteraction = function(value) {
            if (missing(value)) return(private$.crossbasisInteraction)
            if (!("Crossbasis" %in% class(value)))
                stop("ERROR: Unallowed object (", class(value), ") for 'crossbasisInteraction' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.crossbasisInteraction <- value
            return(self)
        },
        crossbasisExposure = function(value) {
            if (missing(value)) return(private$.crossbasisExposure)
            if (!("Crossbasis" %in% class(value)))
                stop("ERROR: Unallowed object (", class(value), ") for 'crossbasisExposure' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.crossbasisExposure <- value
            return(self)
        },
        crossbasis = function(value) {
            stop("Crossbasis cannot be used for CrosspredInteraction!")
        },
        crossbasisName = function(value) {
            stop("CrossbasisName cannot be used for CrosspredInteraction!")
        },
        crossbasisInteractionName = function(value) {
            if (missing(value)) return(private$.crossbasisInteractionName)
            if (!(base::is.character(value)))
                stop("ERROR: Unallowed property ", value, " for 'crossbasisInteractionName' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.crossbasisInteractionName <- value
            return(self)
        },
        crossbasisExposureName = function(value) {
            if (missing(value)) return(private$.crossbasisExposureName)
            if (!(base::is.character(value)))
                stop("ERROR: Unallowed property ", value, " for 'crossbasisExposureName' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.crossbasisExposureName <- value
            return(self)
        },
        coef = function(value) {
            stop("coef cannot be used for CrosspredInteraction!")
        },
        vcov = function(value) {
            stop("vcov cannot be used for CrosspredInteraction!")
        },
        coefInteraction = function(value) {
            if (missing(value)) return(private$.coefInteraction)
            if (!(base::is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'coefInteraction' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.coefInteraction <- value
            return(self)
        },
        vcovInteraction = function(value) {
            if (missing(value)) return(private$.vcovInteraction)
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'vcovInteraction' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.vcovInteraction <- value
            return(self)
        },
        coefExposure = function(value) {
            if (missing(value)) return(private$.coefExposure)
            if (!(base::is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'coefExposure' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.coefExposure <- value
            return(self)
        },
        vcovExposure = function(value) {
            if (missing(value)) return(private$.vcovExposure)
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'vcovExposure' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.vcovExposure <- value
            return(self)
        },
        vcovIntersection = function(value) {
            if (missing(value)) return(private$.vcovIntersection)
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'vcovIntersection' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.vcovIntersection <- value
            return(self)
        }
    )
)

