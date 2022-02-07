#' @rdname CrosspredInteraction
CrosspredInteraction.privateFunctions <- list(

    calcVcov = function(cbNameCols, cbCols, cbNameRows = NULL, cbRows = NULL) {
        if (is.null(cbNameRows)) cbNameRows <- cbNameCols
        if (is.null(cbRows)) cbRows <- cbCols
        vcov <- private$getvcov(self$model, class(self$model))
        possNamesCols <- private$possibleNames(cbNameCols, cbCols)
        possNamesRows <- private$possibleNames(cbNameRows, cbRows)
        indicesRows <- private$mkIndexVector(rownames(vcov), possNamesRows)
        indicesCols <- private$mkIndexVector(colnames(vcov), possNamesCols)
        vcovCut <- vcov[indicesRows, indicesCols, drop = FALSE]
        return(vcovCut)
    },

    #' Prediction of lag-specific effects
    predOfLag = function() {

        ## create the matrix of transformed centred variables. basislags are equal for both crossbases
        predlag <- as.vector(self$crossbasisInteraction$basislag$input)
        XpredInteraction <- private$mkXpred(self$crossbasisInteraction, self$at, predlag, self$cen)
        XpredExposure <- private$mkXpred(self$crossbasisExposure, self$at, predlag, self$cen)
        ## Xpred <- cbind(XpredExposure, XpredInteraction)
        Xpred <- XpredInteraction


        ## coefNu <- c(self$coefExposure, self$coefInteraction)
        ## R <- self$crossbasisInteraction$basislag$x
        ## Gamma_l <- self$crossbasisInteraction$basislag$dimension[[1]]
        ## Gamma_x <- self$crossbasisInteraction$basisvar$dimension[[2]]
        ## nrLags <- length(self$crossbasisInteraction$lags)
        ## I <- diag(2 * Gamma_x)
        ## IR <- kronecker(I, R)
        ## beta <- IR %*% coefNu


        ## Create lag-specific effects and SE
        self$matfitInteraction <- matrix(XpredInteraction %*% self$coefInteraction, base::length(self$at), base::length(predlag))
        self$matfitExposure <- matrix(XpredExposure %*% self$coefExposure, base::length(self$at), base::length(predlag))
        ##:ess-bp-start::conditional@:##
browser(expr={TRUE})##:ess-bp-end:##
        self$matfit <- matrix((Xpred %*% self$coef), base::length(self$at), base::length(predlag))

        self$matseInteraction <- matrix(sqrt(pmax(0, rowSums((XpredInteraction %*% self$vcovInteraction) * XpredInteraction))), base::length(self$at), base::length(predlag))
        self$matseExposure <- matrix(sqrt(pmax(0, rowSums((XpredExposure %*% self$vcovExposure) * XpredExposure))), base::length(self$at), base::length(predlag))
        self$matse <- matrix(sqrt(pmax(0, rowSums((Xpred %*% self$vcov) * Xpred))), base::length(self$at), base::length(predlag))

        ## Names
        rownames(self$matfitInteraction) <- rownames(self$matseInteraction) <- rownames(self$matfitExposure) <- rownames(self$matseExposure) <- self$at
        rownames(self$matfit) <- rownames(self$matse) <- self$at
        colnames(self$matfitInteraction) <- colnames(self$matseInteraction) <- colnames(self$matfitExposure) <- colnames(self$matseExposure) <- outer("lag", predlag, paste, sep = "")
        colnames(self$matfit) <- colnames(self$matse) <- outer("lag", predlag, paste, sep = "")
    },

    #' Prediction of overall+cumulative effects
    predOfOverall = function() {

        ## Create the matrix of transformed variables (dependent on type)
        predlag <- as.vector(self$crossbasisInteraction$basislag$input)
        XpredInteraction <- private$mkXpred(self$crossbasisInteraction, self$at, predlag, self$cen)
        XpredExposure <- private$mkXpred(self$crossbasisExposure, self$at, predlag, self$cen)
        Xpred <- cbind(XpredExposure, XpredInteraction)

        ## Create overall and (optional) cumulative effects and se
        XpredallInteraction <- 0
        XpredallExposure <- 0
        Xpredall <- 0
        if (self$cumul) {
            self$cumfitInteraction <- self$cumseInteraction <- self$cumfitExposure <- self$cumseExposure <- self$cumfit <- self$cumse <- matrix(0, base::length(self$at), base::length(predlag))
        }
        for (i in seq(base::length(predlag))) {
            ind <- seq(base::length(self$at)) + base::length(self$at) * (i - 1)
            XpredallInteraction <- XpredallInteraction + XpredInteraction[ind, , drop = FALSE]
            XpredallExposure <- XpredallExposure + XpredExposure[ind, , drop = FALSE]
            Xpredall <- Xpredall + Xpred[ind, , drop = FALSE]
            if (self$cumul) {
                self$cumfitInteraction[, i] <- XpredallInteraction %*% self$coefInteraction
                self$cumseInteraction[, i] <- sqrt(pmax(0, rowSums((XpredallInteraction %*% self$vcovInteraction) * XpredallInteraction)))
                self$cumfitExposure[, i] <- XpredallExposure %*% self$coefExposure
                self$cumseExposure[, i] <- sqrt(pmax(0, rowSums((XpredallExposure %*% self$vcovExposure) * XpredallExposure)))
                self$cumfit[, i] <- Xpredall %*% self$coef
                self$cumse[, i] <- sqrt(pmax(0, rowSums((Xpredall %*% self$vcov) * Xpredall)))
          }
        }
        self$allfitInteraction <- as.vector(XpredallInteraction %*% self$coefInteraction)
        self$allseInteraction <- sqrt(pmax(0, rowSums((XpredallInteraction %*% self$vcovInteraction) * XpredallInteraction)))
        self$allfitExposure <- as.vector(XpredallExposure %*% self$coefExposure)
        self$allseExposure <- sqrt(pmax(0, rowSums((XpredallExposure %*% self$vcovExposure) * XpredallExposure)))
        self$allfit <- as.vector(Xpredall %*% self$coef)
        self$allse <- sqrt(pmax(0, rowSums((Xpredall %*% self$vcov) * Xpredall)))

        ## Names
        names(self$allfitInteraction) <- names(self$allseInteraction) <- names(self$allfitExposure) <- names(self$allseExposure) <- self$at
        if (self$cumul) {
          rownames(self$cumfitInteraction) <- rownames(self$cumseInteraction) <- rownames(self$cumfitExposure) <- rownames(self$cumseExposure) <- rownames(self$cumfit) <- rownames(self$cumse) <- self$at
          colnames(self$cumfitInteraction) <- colnames(self$cumseInteraction) <- colnames(self$cumfitExposure) <- colnames(self$cumseExposure) <- colnames(self$cumfit) <- colnames(self$cumse) <- outer("lag", predlag, paste, sep = "")
        }
    },

    computeNormedValues = function() {
        z <- qnorm(1 - (1 - self$ci.level) / 2)
        if (!base::is.null(self$model.link) && self$model.link %in% c("log", "logit")) {
            self$matRRfitInteraction <- exp(self$matfitInteraction)
            self$matRRlowInteraction <- exp(self$matfitInteraction - z * self$matseInteraction)
            self$matRRhighInteraction <- exp(self$matfitInteraction + z * self$matseInteraction)
            self$allRRfitInteraction <- exp(self$allfitInteraction)
            self$allRRlowInteraction <- exp(self$allfitInteraction - z * self$allseInteraction)
            names(self$allRRlowInteraction) <- names(self$allfitInteraction)
            self$allRRhighInteraction <- exp(self$allfitInteraction + z * self$allseInteraction)
            names(self$allRRhighInteraction) <- names(self$allfitInteraction)
            if (self$cumul) {
                self$cumRRfitInteraction <- exp(self$cumfitInteraction)
                self$cumRRlowInteraction <- exp(self$cumfitInteraction - z * self$cumseInteraction)
                self$cumRRhighInteraction <- exp(self$cumfitInteraction + z * self$cumseInteraction)
            }
            self$matRRfitExposure <- exp(self$matfitExposure)
            self$matRRlowExposure <- exp(self$matfitExposure - z * self$matseExposure)
            self$matRRhighExposure <- exp(self$matfitExposure + z * self$matseExposure)
            self$allRRfitExposure <- exp(self$allfitExposure)
            self$allRRlowExposure <- exp(self$allfitExposure - z * self$allseExposure)
            names(self$allRRlowExposure) <- names(self$allfitExposure)
            self$allRRhighExposure <- exp(self$allfitExposure + z * self$allseExposure)
            names(self$allRRhighExposure) <- names(self$allfitExposure)
            if (self$cumul) {
                self$cumRRfitExposure <- exp(self$cumfitExposure)
                self$cumRRlowExposure <- exp(self$cumfitExposure - z * self$cumseExposure)
                self$cumRRhighExposure <- exp(self$cumfitExposure + z * self$cumseExposure)
            }
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
            self$matRRfitInteraction <- self$matfitInteraction
            self$matRRlowInteraction <- self$matfitInteraction - z * self$matseInteraction
            self$matRRhighInteraction <- self$matfitInteraction + z * self$matseInteraction
            self$allRRfitInteraction <- self$allfitInteraction
            self$allRRlowInteraction <- self$allfitInteraction - z * self$allseInteraction
            names(self$allRRlowInteraction) <- names(self$allfitInteraction)
            self$allRRhighInteraction <- self$allfitInteraction + z * self$allseInteraction
            names(self$allRRhighInteraction) <- names(self$allfitInteraction)
            if (self$cumul) {
                self$cumRRfitInteraction <- self$cumfitInteraction
                self$cumRRlowInteraction <- self$cumfitInteraction - z * self$cumseInteraction
                self$cumRRhighInteraction <- self$cumfitInteraction + z * self$cumseInteraction
            }
            self$matRRfitExposure <- self$matfitExposure
            self$matRRlowExposure <- self$matfitExposure - z * self$matseExposure
            self$matRRhighExposure <- self$matfitExposure + z * self$matseExposure
            self$allRRfitExposure <- self$allfitExposure
            self$allRRlowExposure <- self$allfitExposure - z * self$allseExposure
            names(self$allRRlowExposure) <- names(self$allfitExposure)
            self$allRRhighExposure <- self$allfitExposure + z * self$allseExposure
            names(self$allRRhighExposure) <- names(self$allfitExposure)
            if (self$cumul) {
                self$cumRRfitExposure <- self$cumfitExposure
                self$cumRRlowExposure <- self$cumfitExposure - z * self$cumseExposure
                self$cumRRhighExposure <- self$cumfitExposure + z * self$cumseExposure
            }
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
)
