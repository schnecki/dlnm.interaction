#' Crosspred class
#'
#' This class defines and computes a prediction for Crosspred
#' @export Crosspred
#' @exportClass Crosspred
Crosspred <- R6::R6Class(
    classname = "Crosspred",

    ## Properties
    private = base::append(base::list(
        .crossbasis = NULL,     # Crossbasis
        .crossbasisName = NULL, # character
        .model = NULL,          # Model
        .model.class = NULL,    # Character
        .link = NULL,           # Character
        .at = NULL,             # Numeric Vector
        .cen = NULL,            # Bool or Numeric
        .coef = NULL,           # vector
        .vcov = NULL,           # matrix
        .cumul = NULL,          # Bool
        .matfit = NULL,         # Matrix
        .matse = NULL,          # Matrix
        .allfit = NULL,         # Vector
        .allse = NULL,          # Vector
        .cumfit = NULL,         # Matrix
        .cumse = NULL,          # Matrix
        .ci.level = NULL,       # numeric
        .matRRfit = NULL,       # Matrix
        .matRRlow = NULL,       # Matrix
        .matRRhigh = NULL,      # Matrix
        .allRRfit = NULL,       # Vector
        .allRRlow = NULL,       # Vector
        .allRRhigh = NULL,      # Vector
        .cumRRfit = NULL,       # Matrix
        .cumRRlow = NULL,       # Matrix
        .cumRRhigh = NULL,      # Matrix
        .knownModels = c("lm", "glm", "gam", "coxph", "lme", "lmerMod", "glmerMod", "lmerModLmerTest", "gee", "geeglm")
    ), Crosspred.privateFunctions),

    ## Methods
    public = list(
        initialize = function(crossbasis, model, at, coef = NULL, vcov = NULL, cen = NULL, ci.level=0.95, cumul = FALSE, model.link=NULL) {
            self$crossbasis <- crossbasis$clone()
            self$crossbasisName <- deparse(substitute(crossbasis))
            if (!base::is.null(coef) && !base::is.null(vcov)) {
                self$coef <- coef
                self$vcov <- vcov
            } else if (base::is.null(model)) {
                stop("Crosspred expects either a class of type model (", private$.knownModels, ") or `coef` and `vcov` have to be specified")
            } else {
                self$model <- model
                self$model.class <- class(model)
                self$link <- private$getLink(model.link)
                self$coef <- private$calcCoef(self$crossbasisName, self$crossbasis)
                self$vcov <- private$calcVcov(self$crossbasisName, self$crossbasis)
            }
            self$at <- at
            self$cen <- private$mkcen(cen, crossbasis)
            self$cumul <- cumul
            self$ci.level <- ci.level
            private$mkPredictions() # create all predictions
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
        crossbasis = function(value) {
            if (missing(value)) return(private$.crossbasis)
            if (!("Crossbasis" %in% class(value)))
                stop("ERROR: Unallowed object (", class(value), ") for 'crossbasis' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.crossbasis <- value
            return(self)
        },
        crossbasisName = function(value) {
            if (missing(value)) return(private$.crossbasisName)
            if (!(base::is.character(value)))
                stop("ERROR: Unallowed property ", value, " for 'crossbasisName' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.crossbasisName <- value
            return(self)
        },
        model = function(value) {
            if (base::missing(value)) return(private$.model)
            if (!(base::any(class(value) %in% private$.knownModels)))
                stop("ERROR: Unallowed property ", value, " for 'model' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.model <- value
            return(self)
        },
        model.class = function(value) {
            if (missing(value)) return(private$.model.class)
            if (!(base::is.character(value)))
                stop("ERROR: Unallowed property ", value, " for 'model.class' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.model.class <- value
            return(self)
        },
        link = function(value) {
            if (missing(value)) return(private$.link)
            if (!(base::is.character(value)))
                stop("ERROR: Unallowed property ", value, " for 'link' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.link <- value
            return(self)
        },
        at = function(value) {
            if (missing(value)) return(private$.at)
            if (!(base::is.numeric(value) && base::is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'at' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.at <- value
            return(self)
        },
        coef = function(value) {
            if (missing(value)) return(private$.coef)
            if (!(base::is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'coef' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.coef <- value
            return(self)
        },
        vcov = function(value) {
            if (missing(value)) return(private$.vcov)
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'vcov' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.vcov <- value
            return(self)
        },
        cen = function(value) {
            if (missing(value)) return(private$.cen)
            if (!(base::is.numeric(value) || base::is.null(value)))
                stop("ERROR: Unallowed property ", value, " for 'cen' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.cen <- value
            return(self)
        },
        matfit = function(value) {
            if (missing(value)) return(private$.matfit)
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'matfit' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.matfit <- value
            return(self)
        },
        matse = function(value) {
            if (missing(value)) return(private$.matse)
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'matse' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.matse <- value
            return(self)
        },
        cumul = function(value) {
            if (missing(value)) return(private$.cumul)
            if (!(base::is.logical(value)))
                stop("ERROR: Unallowed property ", value, " for 'cumul' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.cumul <- value
            return(self)
        },
        allfit = function(value) {
            if (missing(value)) return(private$.allfit)
            if (!(base::is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'allfit' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.allfit <- value
            return(self)
        },
        allse = function(value) {
            if (missing(value)) return(private$.allse)
            if (!(base::is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'allse' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.allse <- value
            return(self)
        },
        cumfit = function(value) {
            if (missing(value)) {
                private$ensureCumulValues()
                return(private$.cumfit)
            }
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'cumfit' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.cumfit <- value
            return(self)
        },
        cumse = function(value) {
            if (missing(value)) {
                private$ensureCumulValues()
                return(private$.cumse)
            }
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'cumse' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.cumse <- value
            return(self)
        },
        ci.level = function(value) {
            if (missing(value)) return(private$.ci.level)
            if (!(base::is.numeric(value)))
                stop("ERROR: Unallowed property ", value, " for 'ci.level' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            if (value >= 1 || value <= 0) stop("'ci.level' must be between 0 and 1")
            private$.ci.level <- value
            return(self)
        },
        matRRfit = function(value) {
            if (missing(value)) return(private$.matRRfit)
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'matRRfit' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.matRRfit <- value
            return(self)
        },
        matRRlow = function(value) {
            if (missing(value)) return(private$.matRRlow)
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'matRRlow' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.matRRlow <- value
            return(self)
        },
        matRRhigh = function(value) {
            if (missing(value)) return(private$.matRRhigh)
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'matRRhigh' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.matRRhigh <- value
            return(self)
        },
        allRRfit = function(value) {
            if (missing(value)) return(private$.allRRfit)
            if (!(base::is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'allRRfit' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.allRRfit <- value
            return(self)
        },
        allRRlow = function(value) {
            if (missing(value)) return(private$.allRRlow)
            if (!(base::is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'allRRlow' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.allRRlow <- value
            return(self)
        },
        allRRhigh = function(value) {
            if (missing(value)) return(private$.allRRhigh)
            if (!(base::is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'allRRhigh' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.allRRhigh <- value
            return(self)
        },
        cumRRfit = function(value) {
            if (missing(value)) {
                private$ensureCumulValues
                return(private$.cumRRfit)
            }
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'cumRRfit' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.cumRRfit <- value
            return(self)
        },
        cumRRlow = function(value) {
            if (missing(value)) {
                private$ensureCumulValues
                return(private$.cumRRlow)
            }
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'cumRRlow' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.cumRRlow <- value
            return(self)
        },
        cumRRhigh = function(value) {
            if (missing(value)) {
                private$ensureCumulValues
                return(private$.cumRRhigh)
            }
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'cumRRhigh' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.cumRRhigh <- value
            return(self)
        }
    )
)

