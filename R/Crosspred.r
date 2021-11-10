#' Crosspred class
#'
#' This class defines a Crosspred
#' @export
Crosspred <- R6::R6Class(
    classname = "Crosspred",

    ## Properties
    private = list(
        .basis = NULL,         # Basis
        .basisName = NULL,     # character
        .model = NULL,         # Model
        .at = NULL,            # Numeric Vector
        .coef = NULL,          # vector
        .vcov = NULL,          # matrix
        .knownModels = c("lm", "glm", "gam", "coxph", "lme", "lmerMod", "glmerMod", "lmerModLmerTest", "gee", "geeglm")
    ),

    ## Methods
    public = list(
        initialize = function(basis, model, at, coef = NULL, vcov = NULL) {
            self$basis <- basis
            self$basisName <- deparse(substitute(basis))
            if (!is.null(coef) && !is.null(vcov)) {
                self$coef <- coef
                self$vcov <- vcov
            } else if (is.null(model)) {
                stop("Crosspred expects either a class of type model (", private$.knownModels, ") or `coef` and `vcov` have to be specified")
            } else {
                self$model <- model
                self$coef <- calcCoef(self)
                self$vcov <- calcVcov(self)
            }
            self$at <- at
            stop("TODO")
        }
    ),

    ## Accessable properties. Active bindings look like fields, but each time they are accessed,
    ## they call a function. They are always publicly visible.
    active = list(
        basis = function(value) {
            if (missing(value)) return(private$.basis)
            if (!("Crossbasis" %in% class(value)))
                stop("ERROR: Unallowed object (", class(value), ") for 'basis' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.basis <- value
            return(self)
        },
        basisName = function(value) {
            if (missing(value)) return(private$.basisName)
            if (!(is.character(value)))
                stop("ERROR: Unallowed property ", value, " for 'basisName' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.basisName <- value
            return(self)
        },
        model = function(value) {
            if (base::missing(value)) return(private$.model)
            if (!(base::any(class(value) %in% private$.knownModels)))
                stop("ERROR: Unallowed property ", value, " for 'model' at ", getSrcFilename(function(){}), ":", getSrcLocation(function(){}))
            private$.model <- value
            return(self)
        },
        at = function(value) {
            if (missing(value)) return(private$.at)
            if (!(is.numeric(value)))
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
        }
    )
)

#' Calculate coef
calcCoef <- function(cp) {
    if (!is.null(cp$coef)) return(cp$coef)
    coef <- getcoef(cp$model, class(cp$model))
    possNames <- possibleNames(cp$basisName, cp$basis)
    indices <- mkIndexVector(names(coef), possNames)
    coef <- coef[indices]
    return(coef)
}

#' Calculate vcov
calcVcov <- function(cp) {
    if (!is.null(cp$vcov)) return(cp$vcov)
    vcov <- getvcov(cp$model, class(cp$model))
    possNames <- possibleNames(cp$basisName, cp$basis)
    indices <- mkIndexVector(rownames(vcov), possNames)
    vcov <- vcov[indices, indices, drop = FALSE]
    return(vcov)

}

#' Generate list of possible names
possibleNames <- function(name, basis) {
    possNames <- list()
    if (basis$dimension[2] == 1) {
        possNames <- append(possNames, name)
    } else {
        for (i in seq(1, basis$dimension[2]))
            possNames <- append(possNames, paste(name, "$x", i, sep = ""))
    }
    return(possNames)
}


mkIndexVector <- function(names, possNames) {
    ls <- list()
    for (idx in seq(1, length(names)))
        if (base::any(names[idx] %in% possNames)) ls <- append(ls, idx)
    return(unlist(ls))
}


getcoef <- function(model, class) {
    ## Extract coef
    ## NB: gam, gee AND geeglm HAVE CLASS glm AS WELL
    coef <-
        if (base::any(class %in% c("glm", "gam", "coxph"))) coef(model)
        else if (base::any(class %in% c("lme", "lmerMod", "glmerMod", "lmerModLmerTest"))) fixef(model)
        else tryCatch(coef(model), error = function(w) "error")
    if (identical(coef, "error"))
        stop("methods for coef() and vcov() must exist for the class of object 'model'. If not, provide the arguments 'coef' and 'vcov'")
    return(coef)
}

getvcov <- function(model, class) {
    ## Extract vcov
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
