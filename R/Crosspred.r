#' Crosspred class
#'
#' This class defines and computes a prediction for Crosspred
#' @export
Crosspred <- R6::R6Class(
    classname = "Crosspred",

    ## Properties
    private = list(
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
    ),

    ## Methods
    public = list(
        initialize = function(crossbasis, model, at, coef = NULL, vcov = NULL, cen = NULL, ci.level=0.95, cumul = FALSE, model.link=NULL) {
            self$crossbasis <- crossbasis$clone()
            self$crossbasisName <- deparse(substitute(crossbasis))
            if (!is.null(coef) && !is.null(vcov)) {
                self$coef <- coef
                self$vcov <- vcov
            } else if (is.null(model)) {
                stop("Crosspred expects either a class of type model (", private$.knownModels, ") or `coef` and `vcov` have to be specified")
            } else {
                self$model <- model
                self$model.class <- class(model)
                self$link <- getLink(self, model.link)
                self$coef <- calcCoef(self)
                self$vcov <- calcVcov(self)
            }
            self$at <- at
            self$cen <- mkcen(cen, crossbasis)
            self$cumul <- cumul
            self$ci.level <- ci.level
            mkPredictions(self) # create all predictionsx
        },
        #' Exposure-RelativeRisk plot. See `plotExposure`
        plotExposure = function(filepath = NULL, title = NULL) plotExposure(self, filepath, title),
        #' Plot lag. See `plotLagAt`.
        plotLagAt = function(at, filepath = NULL, title = NULL, xDist = 0.10, colours = colourPalette, legendPos = c(0.5, 0.10)) {
            if (!is.numeric(at)) stop("Expecting numeric for at in function plotLagAt")
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
            if (!(is.character(value)))
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
            if (!(is.character(value)))
                stop("ERROR: Unallowed property ", value, " for 'model.class' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.model.class <- value
            return(self)
        },
        link = function(value) {
            if (missing(value)) return(private$.link)
            if (!(is.character(value)))
                stop("ERROR: Unallowed property ", value, " for 'link' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.link <- value
            return(self)
        },
        at = function(value) {
            if (missing(value)) return(private$.at)
            if (!(is.numeric(value) && is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'at' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.at <- value
            return(self)
        },
        coef = function(value) {
            if (missing(value)) return(private$.coef)
            if (!(is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'coef' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.coef <- value
            return(self)
        },
        vcov = function(value) {
            if (missing(value)) return(private$.vcov)
            if (!(is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'vcov' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.vcov <- value
            return(self)
        },
        cen = function(value) {
            if (missing(value)) return(private$.cen)
            if (!(is.numeric(value) || is.null(value)))
                stop("ERROR: Unallowed property ", value, " for 'cen' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.cen <- value
            return(self)
        },
        matfit = function(value) {
            if (missing(value)) return(private$.matfit)
            if (!(is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'matfit' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.matfit <- value
            return(self)
        },
        matse = function(value) {
            if (missing(value)) return(private$.matse)
            if (!(is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'matse' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.matse <- value
            return(self)
        },
        cumul = function(value) {
            if (missing(value)) return(private$.cumul)
            if (!(is.logical(value)))
                stop("ERROR: Unallowed property ", value, " for 'cumul' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.cumul <- value
            return(self)
        },
        allfit = function(value) {
            if (missing(value)) return(private$.allfit)
            if (!(is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'allfit' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.allfit <- value
            return(self)
        },
        allse = function(value) {
            if (missing(value)) return(private$.allse)
            if (!(is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'allse' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.allse <- value
            return(self)
        },
        cumfit = function(value) {
            if (missing(value)) {
                ensureCumulValues(self)
                return(private$.cumfit)
            }
            if (!(is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'cumfit' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.cumfit <- value
            return(self)
        },
        cumse = function(value) {
            if (missing(value)) {
                ensureCumulValues(self)
                return(private$.cumse)
            }
            if (!(is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'cumse' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.cumse <- value
            return(self)
        },
        ci.level = function(value) {
            if (missing(value)) return(private$.ci.level)
            if (!(is.numeric(value)))
                stop("ERROR: Unallowed property ", value, " for 'ci.level' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            if (value >= 1 || value <= 0) stop("'ci.level' must be between 0 and 1")
            private$.ci.level <- value
            return(self)
        },
        matRRfit = function(value) {
            if (missing(value)) return(private$.matRRfit)
            if (!(is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'matRRfit' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.matRRfit <- value
            return(self)
        },
        matRRlow = function(value) {
            if (missing(value)) return(private$.matRRlow)
            if (!(is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'matRRlow' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.matRRlow <- value
            return(self)
        },
        matRRhigh = function(value) {
            if (missing(value)) return(private$.matRRhigh)
            if (!(is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'matRRhigh' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.matRRhigh <- value
            return(self)
        },
        allRRfit = function(value) {
            if (missing(value)) return(private$.allRRfit)
            if (!(is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'allRRfit' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.allRRfit <- value
            return(self)
        },
        allRRlow = function(value) {
            if (missing(value)) return(private$.allRRlow)
            if (!(is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'allRRlow' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.allRRlow <- value
            return(self)
        },
        allRRhigh = function(value) {
            if (missing(value)) return(private$.allRRhigh)
            if (!(is.vector(value)))
                stop("ERROR: Unallowed property ", value, " for 'allRRhigh' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.allRRhigh <- value
            return(self)
        },
        cumRRfit = function(value) {
            if (missing(value)) {
                ensureCumulValues(self)
                return(private$.cumRRfit)
            }
            if (!(is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'cumRRfit' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.cumRRfit <- value
            return(self)
        },
        cumRRlow = function(value) {
            if (missing(value)) {
                ensureCumulValues(self)
                return(private$.cumRRlow)
            }
            if (!(is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'cumRRlow' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.cumRRlow <- value
            return(self)
        },
        cumRRhigh = function(value) {
            if (missing(value)) {
                ensureCumulValues(self)
                return(private$.cumRRhigh)
            }
            if (!(is.matrix(value)))
                stop("ERROR: Unallowed property ", value, " for 'cumRRhigh' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.cumRRhigh <- value
            return(self)
        }
    )
)

#' Calculate coef
calcCoef <- function(cp) {
    if (!is.null(cp$coef)) return(cp$coef)
    coef <- getcoef(cp$model, class(cp$model))
    possNames <- possibleNames(cp$crossbasisName, cp$crossbasis)
    indices <- mkIndexVector(names(coef), possNames)
    coef <- coef[indices]
    return(coef)
}

#' Calculate vcov
calcVcov <- function(cp) {
    if (!is.null(cp$vcov)) return(cp$vcov)
    vcov <- getvcov(cp$model, class(cp$model))
    possNames <- possibleNames(cp$crossbasisName, cp$crossbasis)
    indices <- mkIndexVector(rownames(vcov), possNames)
    vcov <- vcov[indices, indices, drop = FALSE]
    return(vcov)
}

#' Get Link
getLink <- function(cp, model.link=NULL) {
    if (!is.null(model.link)) return(model.link)
    class <- class(cp$model)
    model <- cp$model

    ## OTHERWISE, EXTRACT FROM MODEL (IF AVAILABLE)
    link <- if (base::all(class %in% c("lm")) || base::all(class %in% c("lme")) || base::any(class %in% "nlme") || base::any(class %in% "lmerMod")) "identity"
            else if (base::any(class %in% c("clogit"))) "logit"
            else if (base::any(class %in% c("coxph"))) "log"
            else if (base::any(class %in% c("glm")) || base::any(class %in% c("glmmPQL"))) model$family$link
            else if (base::any(class %in% c("glmerMod"))) model@resp$family$link
            else NA
    return(link)
}


#' Generate list of possible names
possibleNames <- function(name, crossbasis) {
    possNames <- list()
    if (crossbasis$dimension[2] == 1) {
        possNames <- append(possNames, name)
    } else {
        for (i in seq(1, crossbasis$dimension[2]))
            possNames <- append(possNames, paste(name, "$x", i, sep = ""))
    }
    return(possNames)
}

#' Vector for finding the index of the column and row names of coef/vcov
mkIndexVector <- function(names, possNames) {
    ls <- list()
    for (idx in seq(1, length(names)))
        if (base::any(names[idx] %in% possNames)) ls <- append(ls, idx)
    return(unlist(ls))
}

# Extract coef from model
getcoef <- function(model, class) {
    ## NB: gam, gee AND geeglm HAVE CLASS glm AS WELL
    coef <-
        if (base::any(class %in% c("glm", "gam", "coxph"))) coef(model)
        else if (base::any(class %in% c("lme", "lmerMod", "glmerMod", "lmerModLmerTest"))) fixef(model)
        else tryCatch(coef(model), error = function(w) "error")
    if (identical(coef, "error"))
        stop("methods for coef() and vcov() must exist for the class of object 'model'. If not, provide the arguments 'coef' and 'vcov'")
    return(coef)
}

## Extract vcov from model
getvcov <- function(model, class) {
    ## NB: gam, gee AND geeglm HAVE CLASS glm AS WELL
    vcov <-
        if (base::any(class %in% c("lm", "glm", "lme", "coxph")) && !base::any(class %in% c("gee"))) vcov(model)
        else if (identical(class, c("gee", "glm"))) model$robust.variance
        else if (base::any(class %in% c("geeglm"))) summary(model)$cov.scaled
        else if (base::any(class %in% c("lmerMod", "glmerMod", "lmerModLmerTest"))) as.matrix(vcov(model))
        else tryCatch(vcov(model), error = function(w) "error")
    if (identical(vcov, "error"))
        stop("methods for coef() and vcov() must exist for the class of object 'model'. If not, provide the arguments 'coef' and 'vcov'")
    return(vcov)
}

#' Reimplementation of original function. Trying to keep the behaviour of dlnm.
mkcen <- function(cen, crossbasis) {

    ## If null, try to extract it from crossbasis
    if (nocen <- is.null(cen)) cen <- crossbasis$basisvar$cen

    ## set centering depending on function and cen type
    if (crossbasis$basisvar$fun %in% c("Thr", "Strata", "Integer", "Lin")) {
        if (is.logical(cen)) cen <- NULL
    } else {
        ## if null or true, set to (approximately) mid-range and message
        if (is.null(cen) || (is.logical(cen) && cen)) cen <- median(pretty(crossbasis$range))
        ## if false, null
        if (is.logical(cen) && !cen) cen <- NULL
    }

    ## however, if intercept is present, set to null
    int <- crossbasis$basisvar$intercept
    if (is.logical(int) && int) cen <- NULL

    ## message
    if (nocen && !is.null(cen))
        message("centering value unspecified. Automatically set to ", cen)
    return(cen)
}

#' Tensor multiplication
mkXpred <- function(crossbasis, predvar, predlag, cen) {

    ## Create vectorized lagged values
    varvec <- rep(predvar, length(predlag))
    lagvec <- rep(predlag, each = length(predvar))

    ## create marginal basis and call tensor
    ## NB: order of basis matrices in tensor changed since version 2.2.4 centering applied
    ##     only marginally to var dimension
    basisvar <- crossbasis$basisvar$mkNewWith(varvec)$x
    basislag <- crossbasis$basislag$mkNewWith(lagvec)$x

    if (!is.null(cen)) {
        basiscen <- crossbasis$basisvar$mkNewWith(cen)$x
        basisvar <- scale(basisvar, center = basiscen, scale = FALSE)
    }
    Xpred <- tensor.prod.model.matrix(list(basisvar, basislag))
    return(Xpred)
}

#' Prediction of lag-specific effects
predOfLag <- function(self) {
    ## create the matrix of transformed centred variables
    predlag <- as.vector(self$crossbasis$basislag$input)
    Xpred <- mkXpred(self$crossbasis, self$at, predlag, self$cen)

    ## Create lag-specific effects and SE
    self$matfit <- matrix(Xpred %*% self$coef, length(self$at), length(predlag))
    self$matse <- matrix(sqrt(pmax(0, rowSums((Xpred %*% self$vcov) * Xpred))), length(self$at), length(predlag))

    ## Names
    rownames(self$matfit) <- rownames(self$matse) <- self$at
    colnames(self$matfit) <- colnames(self$matse) <- outer("lag", predlag, paste, sep = "")


}

#' Prediction of overall+cumulative effects
predOfOverall <- function(self) {

    ## Create the matrix of transformed variables (dependent on type)
    predlag <- as.vector(self$crossbasis$basislag$input)
    Xpred <- mkXpred(self$crossbasis, self$at, predlag, self$cen)

    ## Create overall and (optional) cumulative effects and se
    Xpredall <- 0
    if (self$cumul) {
        self$cumfit <- self$cumse <- matrix(0, length(self$at), length(predlag))
    }
    for (i in seq(length(predlag))) {
        ind <- seq(length(self$at)) + length(self$at) * (i - 1)
        Xpredall <- Xpredall + Xpred[ind, , drop = FALSE]
        if (self$cumul) {
            self$cumfit[, i] <- Xpredall %*% self$coef
            self$cumse[, i] <- sqrt(pmax(0, rowSums((Xpredall %*% self$vcov) * Xpredall)))
      }
    }
    self$allfit <- as.vector(Xpredall %*% self$coef)
    self$allse <- sqrt(pmax(0, rowSums((Xpredall %*% self$vcov) * Xpredall)))

    ## Names
    names(self$allfit) <- names(self$allse) <- self$at
    if (self$cumul) {
      rownames(self$cumfit) <- rownames(self$cumse) <- self$at
      colnames(self$cumfit) <- colnames(self$cumse) <- outer("lag", predlag, paste, sep = "")
    }
}

computeNormedValues <- function(self) {
    z <- qnorm(1 - (1 - self$ci.level) / 2)
    if (!is.null(self$model.link) && self$model.link %in% c("log", "logit")) {
        self$matRRfit <- exp(self$matfit)
        self$matRRlow <- exp(self$matfit - z * self$matse)
        self$matRRhigh <- exp(self$matfit + z * self$matse)
        self$allRRfit <- exp(self$allfit)
        self$allRRlow <- exp(self$allfit - z * self$allse)
        names(self$allRRlow) <- names(self$allfit)
        self$allRRhigh <- exp(self$allfit + z * self$allse)
        names(self$allRRhigh) <- names(self$allfit)
        if (self$cumul) {
            self$cumRRfit <- exp(self$cumfit)
            self$cumRRlow <- exp(self$cumfit - z * self$cumse)
            self$cumRRhigh <- exp(self$cumfit + z * self$cumse)
        }
    } else {
        ## RR or No RR?
        self$matRRfit <- self$matfit
        self$matRRlow <- self$matfit - z * self$matse
        self$matRRhigh <- self$matfit + z * self$matse
        self$allRRfit <- self$allfit
        self$allRRlow <- self$allfit - z * self$allse
        names(self$allRRlow) <- names(self$allfit)
        self$allRRhigh <- self$allfit + z * self$allse
        names(self$allRRhigh) <- names(self$allfit)
        if (self$cumul) {
            self$cumRRfit <- self$cumfit
            self$cumRRlow <- self$cumfit - z * self$cumse
            self$cumRRhigh <- self$cumfit + z * self$cumse
        }
    }
}

ensureCumulValues <- function(self) {
    if (!self$cumul) { # Compute if it was set to FALSE
        self$cumul <- TRUE
        mkPredictions(self)
    }
}

mkPredictions <- function(self) {
    ## Prediction of lag-specific effects
    predOfLag(self)
    predOfOverall(self)
    computeNormedValues(self)
}
