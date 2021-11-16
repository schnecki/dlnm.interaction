
#' @rdname Crosspred
Crosspred.privateFunctions <- base::list(

    #' Generate list of possible names
    possibleNames = function(name, crossbasis) {
        possNames <- base::list()
        if (crossbasis$dimension[2] == 1) {
            possNames <- base::append(possNames, name)
        } else {
            for (i in seq(1, crossbasis$dimension[2]))
                possNames <- base::append(possNames, paste(name, "$x", i, sep = ""))
        }
        return(possNames)
    },

    #' Calculate coef
    calcCoef = function(cbName, cb) {
        coef <- private$getcoef(self$model, class(self$model))
        possNames <- private$possibleNames(cbName, cb)
        indices <- private$mkIndexVector(names(coef), possNames)
        coef <- coef[indices]
        return(coef)
    },

    ## #' Calculate vcov
    ## calcVcov = function(cbName, cb) {
    ##     vcov <- private$getvcov(self$model, class(self$model))
    ##     possNames <- private$possibleNames(cbName, cb)
    ##     indices <- private$mkIndexVector(rownames(vcov), possNames)
    ##     vcov <- vcov[indices, indices, drop = FALSE]
    ##     return(vcov)
    ## },
    #'  Calculate vcov, with possible to select different rows/cols.
    calcVcov = function(cbNameCols, cbCols, cbNameRows = NULL, cbRows = NULL) {
        if (is.null(cbNameRows)) cbNameRows <- cbNameCols
        if (is.null(cbRows)) cbRows <- cbCols
        vcov <- private$getvcov(self$model, class(self$model))
        possNamesCols <- private$possibleNames(cbNameCols, cbCols)
        possNamesRows <- private$possibleNames(cbNameRows, cbRows)
        indicesRows <- private$mkIndexVector(rownames(vcov), possNamesRows)
        indicesCols <- private$mkIndexVector(colnames(vcov), possNamesCols)
        vcov <- vcov[indicesRows, indicesCols, drop = FALSE]
        return(vcov)
    },


    #' Get Link
    getLink = function(model.link = NULL) {
        if (!base::is.null(model.link)) return(model.link)
        class <- class(self$model)
        model <- self$model

        ## OTHERWISE, EXTRACT FROM MODEL (IF AVAILABLE)
        link <- if (base::all(class %in% c("lm")) || base::all(class %in% c("lme")) || base::any(class %in% "nlme") || base::any(class %in% "lmerMod")) "identity"
                else if (base::any(class %in% c("clogit"))) "logit"
                else if (base::any(class %in% c("coxph"))) "log"
                else if (base::any(class %in% c("glm")) || base::any(class %in% c("glmmPQL"))) model$family$link
                else if (base::any(class %in% c("glmerMod"))) model@resp$family$link
                else NA
        return(link)
    },


    #' Vector for finding the index of the column and row names of coef/vcov
    mkIndexVector = function(names, possNames) {
        ls <- base::list()
        for (idx in seq(1, base::length(names)))
            if (base::any(names[idx] %in% possNames)) ls <- base::append(ls, idx)
        return(unlist(ls))
    },

    # Extract coef from model
    getcoef = function(model, class) {
        ## NB: gam, gee AND geeglm HAVE CLASS glm AS WELL
        coef <-
            if (base::any(class %in% c("glm", "gam", "coxph"))) coef(model)
            else if (base::any(class %in% c("lme", "lmerMod", "glmerMod", "lmerModLmerTest"))) fixef(model)
            else tryCatch(coef(model), error = function(w) "error")
        if (base::identical(coef, "error"))
            stop("methods for coef() and vcov() must exist for the class of object 'model'. If not, provide the arguments 'coef' and 'vcov'")
        return(coef)
    },

    ## Extract vcov from model
    getvcov = function(model, class) {
        ## NB: gam, gee AND geeglm HAVE CLASS glm AS WELL
        vcov <-
            if (base::any(class %in% c("lm", "glm", "lme", "coxph")) && !base::any(class %in% c("gee"))) vcov(model)
            else if (base::identical(class, c("gee", "glm"))) model$robust.variance
            else if (base::any(class %in% c("geeglm"))) summary(model)$cov.scaled
            else if (base::any(class %in% c("lmerMod", "glmerMod", "lmerModLmerTest"))) as.matrix(vcov(model))
            else tryCatch(vcov(model), error = function(w) "error")
        if (base::identical(vcov, "error"))
            stop("methods for coef() and vcov() must exist for the class of object 'model'. If not, provide the arguments 'coef' and 'vcov'")
        return(vcov)
    },

    #' Reimplementation of original function. Trying to keep the behaviour of dlnm.
    mkcen = function(cen, crossbasis) {

        ## If null, try to extract it from crossbasis
        if (nocen <- base::is.null(cen)) cen <- crossbasis$basisvar$cen

        ## set centering depending on function and cen type
        if (crossbasis$basisvar$fun %in% c("Thr", "Strata", "Integer", "Lin")) {
            if (base::is.logical(cen)) cen <- NULL
        } else {
            ## if null or true, set to (approximately) mid-range and message
            if (base::is.null(cen) || (base::is.logical(cen) && cen)) cen <- median(pretty(crossbasis$range))
            ## if false, null
            if (base::is.logical(cen) && !cen) cen <- NULL
        }

        ## however, if intercept is present, set to null
        int <- crossbasis$basisvar$intercept
        if (base::is.logical(int) && int) cen <- NULL

        ## message
        if (nocen && !base::is.null(cen))
            message("centering value unspecified. Automatically set to ", cen)
        return(cen)
    },

    #' Tensor multiplication
    mkXpred = function(crossbasis, predvar, predlag, cen) {

        ## Create vectorized lagged values
        varvec <- rep(predvar, base::length(predlag))
        lagvec <- rep(predlag, each = base::length(predvar))

        ## create marginal basis and call tensor
        ## NB: order of basis matrices in tensor changed since version 2.2.4 centering applied
        ##     only marginally to var dimension
        basisvar <- crossbasis$basisvar$mkNewWith(varvec)$x
        basislag <- crossbasis$basislag$mkNewWith(lagvec)$x

        if (!base::is.null(cen)) {
            basiscen <- crossbasis$basisvar$mkNewWith(cen)$x
            basisvar <- scale(basisvar, center = basiscen, scale = FALSE)
        }
        Xpred <- tensor.prod.model.matrix(list(basisvar, basislag))
        return(Xpred)
    },

    #' Prediction of lag-specific effects
    predOfLag = function() {
        ## create the matrix of transformed centred variables
        predlag <- as.vector(self$crossbasis$basislag$input)
        Xpred <- private$mkXpred(self$crossbasis, self$at, predlag, self$cen)

        ## Create lag-specific effects and SE
        self$matfit <- matrix(Xpred %*% self$coef, base::length(self$at), base::length(predlag))
        self$matse <- matrix(sqrt(pmax(0, rowSums((Xpred %*% self$vcov) * Xpred))), base::length(self$at), base::length(predlag))

        ## Names
        rownames(self$matfit) <- rownames(self$matse) <- self$at
        colnames(self$matfit) <- colnames(self$matse) <- outer("lag", predlag, paste, sep = "")


    },

    #' Prediction of overall+cumulative effects
    predOfOverall = function() {

        ## Create the matrix of transformed variables (dependent on type)
        predlag <- as.vector(self$crossbasis$basislag$input)
        Xpred <- private$mkXpred(self$crossbasis, self$at, predlag, self$cen)

        ## Create overall and (optional) cumulative effects and se
        Xpredall <- 0
        if (self$cumul) {
            self$cumfit <- self$cumse <- matrix(0, base::length(self$at), base::length(predlag))
        }
        for (i in seq(base::length(predlag))) {
            ind <- seq(base::length(self$at)) + base::length(self$at) * (i - 1)
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
    },

    computeNormedValues = function() {
        z <- qnorm(1 - (1 - self$ci.level) / 2)
        if (!base::is.null(self$model.link) && self$model.link %in% c("log", "logit")) {
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
    },

    ensureCumulValues = function() {
        if (!self$cumul) { # Compute if it was set to FALSE
            self$cumul <- TRUE
            private$mkPredictions()
        }
    },

    mkPredictions = function() {
        ## Prediction of lag-specific effects
        private$predOfLag()
        private$predOfOverall()
        private$computeNormedValues()
    }
)
