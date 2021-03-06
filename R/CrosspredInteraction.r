#' Crosspred class
#'
#' This class defines and computes a prediction for Crosspred
#' @export CrosspredInteraction
#' @exportClass CrosspredInteraction
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
        .vcovIntersection = NULL,          # matrix
        .matfitInteraction = NULL,         # Matrix
        .matfitExposure = NULL,            # Matrix
        .matseInteraction = NULL,          # Matrix
        .matseExposure = NULL,             # Matrix
        .allfitInteraction = NULL,         # Vector
        .allseInteraction = NULL,          # Vector
        .cumfitInteraction = NULL,         # Matrix
        .cumseInteraction = NULL,          # Matrix
        .allfitExposure = NULL,            # Vector
        .allseExposure = NULL,             # Vector
        .cumfitExposure = NULL,            # Matrix
        .cumseExposure = NULL,             # Matrix
        .matRRfitInteraction = NULL,       # Matrix
        .matRRlowInteraction = NULL,       # Matrix
        .matRRhighInteraction = NULL,      # Matrix
        .allRRfitInteraction = NULL,       # Vector
        .allRRlowInteraction = NULL,       # Vector
        .allRRhighInteraction = NULL,      # Vector
        .cumRRfitInteraction = NULL,       # Matrix
        .cumRRlowInteraction = NULL,       # Matrix
        .cumRRhighInteraction = NULL,      # Matrix
        .matRRfitExposure = NULL,          # Matrix
        .matRRlowExposure = NULL,          # Matrix
        .matRRhighExposure = NULL,         # Matrix
        .allRRfitExposure = NULL,          # Vector
        .allRRlowExposure = NULL,          # Vector
        .allRRhighExposure = NULL,         # Vector
        .cumRRfitExposure = NULL,          # Matrix
        .cumRRlowExposure = NULL,          # Matrix
        .cumRRhighExposure = NULL          # Matrix


    ), CrosspredInteraction.privateFunctions),

    ## Methods
    public = list(
        initialize = function(crossbasisExposure, crossbasisInteraction, model, at = NULL, quantilesExposure1 = NULL, quantilesExposure2 = NULL, coefInteraction = NULL, coefExposure = NULL, vcovInteraction = NULL, vcovExposure = NULL, vcovIntersection = NULL
                            , cen = NULL, ci.level=0.95, cumul = FALSE, model.link=NULL) {
            if (identical(crossbasisExposure$basislag, crossbasisInteraction$basislag)) stop("The crossbases have to have the same basislag!")
            if (base::is.null(at) && (base::is.null(quantilesExposure1) || base::is.null(quantilesExposure2)))
                stop("You must either specify a value for `at` or a value for both, `quantilesExposure1` and `quantilesExposure2`. Note that the value for `at` should be computed in regard to both dimensions w/o unique!")
            else if (base::is.null(at)) self$at <- unique(quantilesExposure1 * quantilesExposure2)

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
                self$coefInteraction <- private$calcCoef(self$crossbasisInteractionName, self$crossbasisInteraction)
                self$vcovInteraction <- private$calcVcov(self$crossbasisInteractionName, self$crossbasisInteraction)
                self$coefExposure <- private$calcCoef(self$crossbasisExposureName, self$crossbasisExposure)
                self$vcovExposure <- private$calcVcov(self$crossbasisExposureName, self$crossbasisExposure)
                self$vcovIntersection <- private$calcVcov(self$crossbasisExposureName, self$crossbasisExposure, self$crossbasisInteractionName, self$crossbasisInteraction)
            }
            ## Build coef and vcov matrices
            coefNu <- c(self$coefExposure, self$coefInteraction)
            vcovNu <- kronecker(matrix(c(1, 0, 0, 0), nrow = 2), self$vcovExposure) +
                      kronecker(matrix(c(0, 1, 0, 0), nrow = 2), self$vcovIntersection) +
                      kronecker(matrix(c(0, 0, 1, 0), nrow = 2), t(self$vcovIntersection)) +
                      kronecker(matrix(c(0, 0, 0, 1), nrow = 2), self$vcovInteraction)

            ## Calculate beta value
            R <- self$crossbasisInteraction$basislag$x
            Gamma_l <- self$crossbasisInteraction$basislag$dimension[[1]]
            Gamma_x <- self$crossbasisInteraction$basisvar$dimension[[2]]
            nrLags <- length(self$crossbasisInteraction$lags)
            I <- diag(2 * Gamma_x)
            IR <- kronecker(I, R)
            beta <- IR %*% coefNu
            self$vcov <- IR %*% vcovNu %*% t(IR)
            ## beta1Sel <- c(rep(1, Gamma_l * Gamma_x), rep(0, Gamma_l * Gamma_x))
            ## beta12Sel <- c(rep(1, Gamma_l * Gamma_x), rep(1, Gamma_l * Gamma_x))
            beta1Sel <- c(kronecker(diag(Gamma_x), rep(1, nrLags)), c(kronecker(diag(Gamma_x), rep(0, nrLags))))
            beta12Sel <- c(kronecker(diag(Gamma_x), rep(1, nrLags)), c(kronecker(diag(Gamma_x), rep(1, nrLags))))

            beta1Star <- beta1Sel %*% beta
            beta2Star <- beta12Sel %*% beta
            self$coef <- c(beta1Star, beta2Star)


            ## Lo mismo:
            ## beta1Star_ <- c(kronecker(diag(Gamma_x), rep(1, Gamma_l)), kronecker(diag(Gamma_x), rep(0, Gamma_l))) %*% beta
            ## beta2Star_ <- c(kronecker(diag(Gamma_x), rep(1, Gamma_l)), kronecker(diag(Gamma_x), rep(1, Gamma_l))) %*% beta


            ## selLowerBottom <- t(matrix(c(rep(0, 2 * Gamma_l), rep(1, Gamma_l), rep(0, Gamma_l)), nrow = 2 * Gamma_l))

            ## betas <- rep(0, Gamma_l)
            ## gammas <- rep(0, Gamma_l)
            ## offset <- Gamma_x * Gamma_l
            ## for (i in seq_len(Gamma_x)) {
            ##     betas <- betas + beta[seq(1 + (i - 1) * Gamma_l, (i - 1) * Gamma_l + Gamma_l)]
            ##     gammas <- gammas + beta[seq(offset + (i - 1) * Gamma_l + 1, offset + (i - 1) * Gamma_l + Gamma_l)]
            ## }

            ## self$coef <- c(betas, gammas)
            ## attr(self$coef, "names") <- unlist(list(attr(self$coefExposure, "names"), attr(self$coefInteraction, "names")))

            ## sel <- c(rep(c(1, rep(0, Gamma_l - 1)), Gamma_x), rep(0, Gamma_x * Gamma_l))
            # sumOverLagSel <- kronecker(diag(2 * Gamma_x)  + selLowerBottom, t(rep(1, Gamma_l)))
            ## sumOverLagSel <- kronecker(diag(2 * Gamma_l)  + selLowerBottom, t(rep(1, Gamma_l))

            ## self$coef <- as.vector(sumOverLagSel %*% beta)
            # attr(self$coef, "names") <- unlist(list(attr(self$coefExposure, "names"), attr(self$coefInteraction, "names")))

            ## Sum betas
            ## beta1Sel <- c(rep(1, Gamma_l * Gamma_x), rep(0, Gamma_l * Gamma_x))  # Only Exposure
            ## beta12Sel <- c(rep(1, Gamma_l * Gamma_x), rep(1, Gamma_l * Gamma_x)) # Exposure and Interaction
            ## beta1Star <- beta1Sel %*% beta
            ## beta2Star <- beta12Sel %*% beta
            ##
            ## beta1Star_ <- c(kronecker(diag(Gamma_x), rep(1, Gamma_l))), kronecker(diag(Gamma_x), rep(0, Gamma_l))) %*% beta
            ## beta2Star_ <- c(kronecker(diag(Gamma_x), rep(1, Gamma_l)), kronecker(diag(Gamma_x), rep(1, Gamma_l))) %*% beta
            ##
            ## self$coef <- c(beta1Star, beta2Star)
            ## attr(self$coef, "names") <- unlist(list(attr(self$coefExposure, "names"), attr(self$coefInteraction, "names")))

            ## self$vcov <- IR %*% vcovNu %*% t(IR)
            ## dimnames(self$vcov)[[1]] <- unlist(list(dimnames(self$vcovExposure)[[1]], dimnames(self$vcovInteraction)[[1]]))
            ## dimnames(self$vcov)[[2]] <- unlist(list(dimnames(self$vcovExposure)[[2]], dimnames(self$vcovInteraction)[[2]]))


            self$at <- at
            self$cen <- private$mkcen(cen, crossbasisInteraction)
            self$cumul <- cumul
            self$ci.level <- ci.level
            private$mkPredictions() # create all predictionsx
        },
        #' Exposure-RelativeRisk plot. See `plotExposure`
        plotExposure = function(filepath = NULL, title = NULL) plotExposure(self, filepath, title),
        #' Plot lag. See `plotLagAt`.
        plotLagAt = function(at, filepath = NULL, title = NULL, xDist = 0.10, colours = colourPalette, legendPos = c(0.5, 0.10)) {
            if (!base::is.numeric(at)) stop("Expecting numeric for at in function plotLagAt")
            at <- as.vector(at)
            cpAtInteraction <- Crosspred$new(self$crossbasisInteraction, self$model, at, NULL, NULL, self$vcovIntersection, self$cen, self$ci.level, self$cumul, self$model.link)
            cpAtExposure <- Crosspred$new(self$crossbasisExposure, self$model, at, self$coefExposure, self$vcovExposure, self$cen, self$ci.level, self$cumul, self$model.link)
            cpAtInteraction$crossbasisName <- self$crossbasisInteractionName
            cpAtExposure$crossbasisName <- self$crossbasisExposureName
            plotLagAt(cpAtInteraction, filepath, title, xDist, colours, legendPos)
        }
    ),

    ## Accessable properties. Active bindings look like fields, but each time they are accessed,
    ## they call a function. They are always publicly visible.
    active = list(
        crossbasis =     function(value) if (missing(value)) self$crossbasisInteraction else private$.crossbasis <- self$crossbasisInteraction <- value,
        crossbasisName = function(value) if (missing(value)) self$crossbasisInteractionName else private$.crossbasisName <- self$crossbasisInteractionName <- value,
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
        },
        ## matfit = function(value) if (missing(value)) self$matfitInteraction else private$.matfit <- self$matfitInteraction <- value,
        ## matse =  function(value) if (missing(value)) self$matseInteraction else private$.matse <- self$matseInteraction <- value,
        matfitInteraction = function(value) {
            if (missing(value)) return(private$.matfitInteraction)
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'matfitInteraction' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.matfitInteraction <- value
            return(self)
        },
        matseInteraction = function(value) {
            if (missing(value)) return(private$.matseInteraction)
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'matseInteraction' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.matseInteraction <- value
            return(self)
        },
        matfitExposure = function(value) {
            if (missing(value)) return(private$.matfitExposure)
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'matfitExposure' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.matfitExposure <- value
            return(self)
        },
        matseExposure = function(value) {
            if (missing(value)) return(private$.matseExposure)
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'matseExposure' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.matseExposure <- value
            return(self)
        },
        ## allfit = function(value) if (missing(value)) self$allfitInteraction else private$.allfit <- self$allfitInteraction <- value,
        ## cumfit = function(value) if (missing(value)) self$cumfitInteraction else private$.cumfit <- self$cumfitInteraction <- value,
        ## allse =  function(value) if (missing(value)) self$allseInteraction else private$.allse <- self$allseInteraction <- value,
        ## cumse =  function(value) if (missing(value)) self$cumseInteraction else private$.cumse <- self$cumseInteraction <- value,
        allfitInteraction = function(value) {
            if (missing(value)) return(private$.allfitInteraction)
            if (!(base::is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'allfitInteraction' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.allfitInteraction <- value
            return(self)
        },
        allseInteraction = function(value) {
            if (missing(value)) return(private$.allseInteraction)
            if (!(base::is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'allseInteraction' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.allseInteraction <- value
            return(self)
        },
        cumfitInteraction = function(value) {
            if (missing(value)) {
                private$ensureCumulValues()
                return(private$.cumfitInteraction)
            }
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'cumfitInteraction' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.cumfitInteraction <- value
            return(self)
        },
        cumseInteraction = function(value) {
            if (missing(value)) {
                private$ensureCumulValues()
                return(private$.cumseInteraction)
            }
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'cumseInteraction' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.cumseInteraction <- value
            return(self)
        },
        allfitExposure = function(value) {
            if (missing(value)) return(private$.allfitExposure)
            if (!(base::is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'allfitExposure' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.allfitExposure <- value
            return(self)
        },
        allseExposure = function(value) {
            if (missing(value)) return(private$.allseExposure)
            if (!(base::is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'allseExposure' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.allseExposure <- value
            return(self)
        },
        cumfitExposure = function(value) {
            if (missing(value)) {
                private$ensureCumulValues()
                return(private$.cumfitExposure)
            }
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'cumfitExposure' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.cumfitExposure <- value
            return(self)
        },
        cumseExposure = function(value) {
            if (missing(value)) {
                private$ensureCumulValues()
                return(private$.cumseExposure)
            }
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'cumseExposure' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.cumseExposure <- value
            return(self)
        },
        ## matRRfit =  function(value) if (missing(value)) self$matRRfitInteraction  else private$.matRRfit <- self$matRRfitInteraction <- value,
        ## matRRlow =  function(value) if (missing(value)) self$matRRlowInteraction  else private$.matRRlow <- self$matRRlowInteraction <- value,
        ## matRRhigh = function(value) if (missing(value)) self$matRRhighInteraction else private$.matRRhigh <- self$matRRhighInteraction <- value,
        ## allRRfit =  function(value) if (missing(value)) self$allRRfitInteraction  else private$.allRRfit <- self$allRRfitInteraction <- value,
        ## allRRlow =  function(value) if (missing(value)) self$allRRlowInteraction  else private$.allRRlow <- self$allRRlowInteraction <- value,
        ## allRRhigh = function(value) if (missing(value)) self$allRRhighInteraction else private$.allRRhigh <- self$allRRhighInteraction <- value,
        ## cumRRfit =  function(value) if (missing(value)) self$cumRRfitInteraction  else private$.cumRRfit <- self$cumRRfitInteraction <- value,
        ## cumRRlow =  function(value) if (missing(value)) self$cumRRlowInteraction  else private$.cumRRlow <- self$cumRRlowInteraction <- value,
        ## cumRRhigh = function(value) if (missing(value)) self$cumRRhighInteraction else private$.cumRRhigh <- self$cumRRhighInteraction <- value,
        matRRfitInteraction = function(value) {
            if (missing(value)) return(private$.matRRfitInteraction)
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'matRRfitInteraction' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.matRRfitInteraction <- value
            return(self)
        },
        matRRlowInteraction = function(value) {
            if (missing(value)) return(private$.matRRlowInteraction)
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'matRRlowInteraction' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.matRRlowInteraction <- value
            return(self)
        },
        matRRhighInteraction = function(value) {
            if (missing(value)) return(private$.matRRhighInteraction)
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'matRRhighInteraction' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.matRRhighInteraction <- value
            return(self)
        },
        allRRfitInteraction = function(value) {
            if (missing(value)) return(private$.allRRfitInteraction)
            if (!(base::is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'allRRfitInteraction' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.allRRfitInteraction <- value
            return(self)
        },
        allRRlowInteraction = function(value) {
            if (missing(value)) return(private$.allRRlowInteraction)
            if (!(base::is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'allRRlowInteraction' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.allRRlowInteraction <- value
            return(self)
        },
        allRRhighInteraction = function(value) {
            if (missing(value)) return(private$.allRRhighInteraction)
            if (!(base::is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'allRRhighInteraction' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.allRRhighInteraction <- value
            return(self)
        },
        cumRRfitInteraction = function(value) {
            if (missing(value)) {
                private$ensureCumulValues
                return(private$.cumRRfitInteraction)
            }
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'cumRRfitInteraction' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.cumRRfitInteraction <- value
            return(self)
        },
        cumRRlowInteraction = function(value) {
            if (missing(value)) {
                private$ensureCumulValues
                return(private$.cumRRlowInteraction)
            }
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'cumRRlowInteraction' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.cumRRlowInteraction <- value
            return(self)
        },
        cumRRhighInteraction = function(value) {
            if (missing(value)) {
                private$ensureCumulValues
                return(private$.cumRRhighInteraction)
            }
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'cumRRhighInteraction' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.cumRRhighInteraction <- value
            return(self)
        },
        matRRfitExposure = function(value) {
            if (missing(value)) return(private$.matRRfitExposure)
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'matRRfitExposure' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.matRRfitExposure <- value
            return(self)
        },
        matRRlowExposure = function(value) {
            if (missing(value)) return(private$.matRRlowExposure)
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'matRRlowExposure' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.matRRlowExposure <- value
            return(self)
        },
        matRRhighExposure = function(value) {
            if (missing(value)) return(private$.matRRhighExposure)
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'matRRhighExposure' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.matRRhighExposure <- value
            return(self)
        },
        allRRfitExposure = function(value) {
            if (missing(value)) return(private$.allRRfitExposure)
            if (!(base::is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'allRRfitExposure' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.allRRfitExposure <- value
            return(self)
        },
        allRRlowExposure = function(value) {
            if (missing(value)) return(private$.allRRlowExposure)
            if (!(base::is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'allRRlowExposure' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.allRRlowExposure <- value
            return(self)
        },
        allRRhighExposure = function(value) {
            if (missing(value)) return(private$.allRRhighExposure)
            if (!(base::is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'allRRhighExposure' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.allRRhighExposure <- value
            return(self)
        },
        cumRRfitExposure = function(value) {
            if (missing(value)) {
                private$ensureCumulValues
                return(private$.cumRRfitExposure)
            }
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'cumRRfitExposure' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.cumRRfitExposure <- value
            return(self)
        },
        cumRRlowExposure = function(value) {
            if (missing(value)) {
                private$ensureCumulValues
                return(private$.cumRRlowExposure)
            }
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'cumRRlowExposure' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.cumRRlowExposure <- value
            return(self)
        },
        cumRRhighExposure = function(value) {
            if (missing(value)) {
                private$ensureCumulValues
                return(private$.cumRRhighExposure)
            }
            if (!(base::is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'cumRRhighExposure' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.cumRRhighExposure <- value
            return(self)
        }
    )
)
