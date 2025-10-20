#Prediction functions

################## Brier score prediction function for casebase #######################
predictRisk.CompRisk<- function(object, newdata, times, cause, ...) {
    #get all covariates excluding intercept and time
    
    coVars=colnames(object@x[, c(grepl("X", colnames(object@x)))])
    #coVars is used in lines 44 and 50
    newdata=data.matrix(drop(subset(newdata, select=coVars)))
    
    if (missing(cause)) stop("Argument cause should be the event type for which we predict the absolute risk.")
    # the output of absoluteRisk is an array with dimension depending on the length of the requested times:
    # case 1: the number of time points is 1
    #         dim(array) =  (length(time), NROW(newdata), number of causes in the data)
    if (length(times) == 1) {
        a <- absoluteRisk.CompRisk(object, newdata = newdata, time = times, addZero = FALSE)
        p <- matrix(a, ncol = 1)
    } else {
        # case 2 a) zero is included in the number of time points
        if (0 %in% times) {
            # dim(array) =  (length(time)+1, NROW(newdata)+1, number of causes in the data)
            a <- casebase::absoluteRisk.CompRisk(object, newdata = newdata, time = times)
            p <- t(a)
        } else {
            # case 2 b) zero is not included in the number of time points (but the absoluteRisk function adds it)
            a <- casebase::absoluteRisk.CompRisk(object, newdata = newdata, time = times)
            ### we need to invert the plot because, by default, we get cumulative incidence
            #a[, -c(1)] <- 1 - a[, -c(1)]
            ### we remove time 0 for everyone, and remove the time column
            a <- a[-c(1), -c(1)] ### a[-c(1), ] to keep times column, but remove time 0 probabilities
            # now we transpose the matrix because in riskRegression we work with number of
            # observations in rows and time points in columns
            p <- t(a)
        }
    }
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) {
        stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ", 
                   NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ", 
                   NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
    }
    p
}

##################### Prediction function for casebase penalized model ################
predict_CompRisk <- function(object, newdata = NULL) {
    ttob <- terms(object)
    X <- model.matrix(delete.response(ttob), newdata,
                      contrasts = if (length(object@contrasts)) {
                          object@contrasts
                      } else NULL,
                      xlev = object@xlevels)
    coeffs <- matrix(coef(object), nrow = ncol(X),
                     byrow = TRUE)
    preds <- X %*% coeffs
    colnames(preds) <- paste0("log(mu[,",
                              seq(2, length(object@typeEvents)),
                              "]/mu[,1])")
    return(preds)
}

################## Prediction function for coxBoost ###########################
#this version predicts with replicates
# predictRisk.iCoxBoost <- function(object, newdata, times, cause,...){
#     p <- predict(object, newdata= newdata,type="CIF",times= newdata$ftime)
#     # if (is.null(dim(p))) {
#     #    if (length(p)!=length(times))
#     #     stop("Prediction failed (wrong number of times)")
#     # }
#     #  else{
#     #    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
#     #      stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
#     #  }
#     p
# }

predictRisk.iCoxBoost <- function(object, newdata, times, cause, ...) {
    p <- predict(object, newdata = newdata, type = "CIF", times = times)
    
    # handle possible return shapes
    if (is.list(p)) {
        key <- if (!is.null(names(p)) && as.character(cause) %in% names(p)) as.character(cause) else cause
        p <- p[[key]]
    }
    if (length(dim(p)) == 3L) {              
        p <- p[,, cause, drop = TRUE]
    }
    if (is.vector(p)) {                       
        p <- matrix(p, nrow = NROW(newdata), ncol = length(times), byrow = FALSE)
    }
    if (nrow(p) == length(times) && ncol(p) == NROW(newdata)) p <- t(p)
    
    stopifnot(nrow(p) == NROW(newdata), ncol(p) == length(times))
    colnames(p) <- format(times)
    p
}



# predictRisk.penalizedCompRisk <- function(object, newdata, times, cause, ...) {
#     #get all covariates excluding intercept and time
#     #newdata = newdata[, .SD, .SDcols = grep("X", colnames(newdata), value = TRUE)]
#     if (missing(cause)) stop("Argument cause should be the event type for which we predict the absolute risk.")
#     # the output of absoluteRisk is an array with dimension depending on the length of the requested times:
#     # case 1: the number of time points is 1
#     #         dim(array) =  (length(time), NROW(newdata), number of causes in the data)
#     if (length(times) == 1) {
#         a <- absoluteRisk.penalized(object, newdata = newdata, time = times, addZero = FALSE)
#         p <- matrix(a, ncol = 1)
#     } else {
#         # case 2 a) zero is included in the number of time points
#         if (0 %in% times) {
#             # dim(array) =  (length(time)+1, NROW(newdata)+1, number of causes in the data)
#             a <- absoluteRisk.penalized(object, newdata = newdata, time = times)
#             p <- t(a)
#         } else {
#             # case 2 b) zero is not included in the number of time points (but the absoluteRisk function adds it)
#             a <- absoluteRisk.penalized(object, newdata = newdata, time = times)
#             ### we need to invert the plot because, by default, we get cumulative incidence
#             #a[, -c(1)] <- 1 - a[, -c(1)]
#             ### we remove time 0 for everyone, and remove the time column
#             a <- a[-c(1), -c(1)] ### a[-c(1), ] to keep times column, but remove time 0 probabilities
#             # now we transpose the matrix because in riskRegression we work with number of
#             # observations in rows and time points in columns
#             p <- t(a)
#         }
#     }
#     if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) {
#         stop(paste("\nPrediction matrix has wrong dimensions:\nRequested newdata x times: ", 
#                    NROW(newdata), " x ", length(times), "\nProvided prediction matrix: ", 
#                    NROW(p), " x ", NCOL(p), "\n\n", sep = ""))
#     }
#     p
# }

`%||%` <- function(a,b) if (!is.null(a)) a else b

predictRisk.penalizedCompRisk <- function(object, newdata, times, cause, ...) {
    stopifnot(cause == 1)
    
    cb <- object$cb_data
    cn <- colnames(cb$covariates)
    Xnew <- as.matrix(newdata[, cn, drop = FALSE])
    N <- nrow(Xnew); Tt <- length(times)
    
    # --- grab and align coefficients by name ---
    beta_mat <- object$coefficients
    # drop intercept if present
    if (!is.null(rownames(beta_mat)) && "(Intercept)" %in% rownames(beta_mat)) {
        beta_mat <- beta_mat[setdiff(rownames(beta_mat), "(Intercept)"), , drop = FALSE]
    }
    # choose the column you want (here first)
    if (is.matrix(beta_mat)) beta_mat <- beta_mat[, 1, drop = FALSE]
    
    # build beta vector matching covariate columns; fill absent ones with 0
    if (!is.null(rownames(beta_mat))) {
        beta <- setNames(rep(0, length(cn)), cn)
        rn <- intersect(rownames(beta_mat), cn)
        beta[rn] <- as.numeric(beta_mat[rn, 1])
        beta <- as.numeric(beta)  # now length(beta) == length(cn)
    } else {
        # no names available → last resort: require matching length
        stopifnot(length(beta_mat) == length(cn))
        beta <- as.numeric(beta_mat)
    }
    
    # --- map requested times to the cb grid + offsets ---
    gTimes <- sort(unique(cb$time))
    off_by_time <- tapply(cb$offset, cb$time, mean)
    pos <- prodlim::sindex(jump.times = gTimes, eval.times = pmin(times, max(gTimes)))
    off_t <- as.numeric(off_by_time)[pos]  # length Tt
    
    # --- hazard and CIF (single cause) ---
    linp <- as.vector(Xnew %*% beta)
    hmat <- sapply(off_t, function(o) plogis(linp + o))  # N x Tt
    
    CIF <- matrix(0, N, Tt); S <- matrix(1, N, Tt)
    for (j in seq_len(Tt)) {
        add <- if (j == 1) hmat[, j] else S[, j-1] * hmat[, j]
        CIF[, j] <- if (j == 1) add else CIF[, j-1] + add
        S[, j]   <- if (j == 1) (1 - hmat[, j]) else S[, j-1] * (1 - hmat[, j])
    }
    CIF <- pmin(1, pmax(0, CIF))
    
    return(CIF)
}