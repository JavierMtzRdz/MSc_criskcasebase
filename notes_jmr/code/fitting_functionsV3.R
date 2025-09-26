library(future)
library(furrr)
library(progressr)

same <- function (x, y, tolerance = .Machine$double.eps) {
    abs(x - y) < tolerance
}


##################################################################
##                         Project: fit                         ##
##################################################################
 


#' Prepare Penalty Parameters for the Fitting Function
#'
#' A helper function to calculate the lambda parameters required by the cdSCRIP
#' fitting function based on the chosen regularization method.
#'
#' @param regularization A string, either 'elastic-net' or 'SCAD'.
#' @param lambda The primary shrinkage parameter.
#' @param alpha The mixing parameter for elastic-net or shape parameter for SCAD.
#' @return A list containing the calculated lambda1 and lambda2 values.
prepare_penalty_params <- function(regularization, lambda, alpha) {
    # Ensure regularization is a character string
    if (!is.character(regularization) || length(regularization) != 1) {
    }
    switch(
        regularization,
        `elastic-net` = {
            # Default alpha for elastic-net is 0.5
            if (is.null(alpha)) {
                alpha <- 0.5
            }

            list(lambda1 = lambda * alpha, 
                 lambda2 = 0.5 * lambda * (1 - alpha))
        },
        
        
        SCAD = {
            # Default 'a' parameter for SCAD is 3.7
            if (is.null(alpha))  {
                alpha <- 3.7 
                
            } else if (alpha < 2)  {
                alpha <- 3.7 
                
            }
            # Return lambda and the shape parameter
            list(lambda1 = lambda, lambda2 = alpha)
        },
        
       
    )
}

penalty_params_notif <- function(regularization, alpha) {
    # Ensure regularization is a character string
    if (!is.character(regularization) || length(regularization) != 1) {
        cli::cli_abort("Argument 'regularization' must be a single string, like 'elastic-net' or 'SCAD'.")
    }
    switch(
        regularization,
        `elastic-net` = {
            # Default alpha for elastic-net is 0.5
            if (is.null(alpha)) {
                cli::cli_inform("Using default alpha = 0.5 for elastic-net.")
            } else if (alpha < 0 || alpha > 1) {
                cli::cli_abort("'alpha' for elastic-net must be between 0 and 1.")
            }
            
        },
        
        
        SCAD = {
            if (is.null(alpha))  {
                cli::cli_inform("Using default alpha (shape parameter) = 3.7 for SCAD.")
                
            } else if (alpha < 2)  {
                cli::cli_inform("Using default alpha (shape parameter) = 3.7 for SCAD.")
                
            }
        },
        
 
    )
}


#' Create a Case-Base Dataset for Competing Risks
#'
#' Transforms a survival dataset into the case-base format required for fitting
#' a multinomial logistic model for competing risks.
#'
#' @param formula A survival formula, e.g., `Surv(time, event) ~ cov1 + cov2`.
#' @param data A data frame containing the variables in the formula.
#' @param ratio Numeric. The ratio of base samples to case samples.
#' @param ratio_event A numeric vector specifying which event categories to use
#'   for determining the size of the base series. Defaults to "all", which uses
#'   the total number of all non-censored events.
#'
#' @return A list containing the components of the case-base dataset: `time`,
#'   `event`, `covariates`, and `offset`.
#' @export
create_cb_data <- function(formula, data, ratio = 20, ratio_event = "all") {
    
    data <- data.frame(data)
    
    response_vars <- all.vars(formula[[2]])
    
    if (length(response_vars) != 2) {
        cli::cli_abort("LHS of formula must be Surv(time, event).")
    }
    time_var <- response_vars[1]
    status_var <- response_vars[2]
    
    if (!all(c(time_var, status_var) %in% names(data))) {
        cli::cli_abort("Time and/or status variables from formula not found in data.")
    }
    
    # time <- data[[time_var]]
    # status <- data[[status_var]]
    time <- data[,time_var]
    
    status <- data[,status_var]
    
    cov_matrix <- stats::model.matrix(as.formula(formula[-2]),
                                       data[!(names(data) %in% 
                                                 c(status_var,
                                                   time_var))])[, -1, 
                                                                drop = FALSE]
    
    # Determine which events to use for the base series ratio
    
    all_event_types <- unique(status[status != 0])
    if (identical(ratio_event, "all")) {
        event_types_for_cases <- unique(status[status != 0])
    } else if (is.numeric(ratio_event) && all(ratio_event %in% all_event_types)) {
        event_types_for_cases <- ratio_event
        
    } else {
        cli::cli_abort("'ratio_event' must be 'all' or a numeric vector of valid event types.")
    }
    
    case_count <- length(which(status %in% event_types_for_cases))
    
    if (case_count == 0) cli::cli_abort("No events found for the specified 'ratio_event'.")
    
    # Create Base Series
    n <- nrow(data)
    B <- sum(time)
    b_size <- ratio * case_count
    offset <- log(B / b_size)
    
    cli::cli_alert_info("Created base series with {b_size} samples based on {case_count} event(s).")
    
    prob_select <- time / B
    sampled_indices <- sample(n, size = b_size, 
                              replace = TRUE, prob = prob_select)
    
    
    case_indices <- which(status %in% all_event_types)
    
    # Combine and return
    n_b <- length(sampled_indices)
    n_c <- length(case_indices)
    total_rows <- n_b + n_c
    
    # Pre-allocate the final objects 
    final_time <- numeric(total_rows)
    final_event <- numeric(total_rows)
    final_covs <- matrix(0, nrow = total_rows, ncol = ncol(cov_matrix))
    colnames(final_covs) <- colnames(cov_matrix)
    
    # base series
    final_time[1:n_b] <- runif(n_b) * time[sampled_indices]
    final_event[1:n_b] <- 0 
    final_covs[1:n_b, ] <- cov_matrix[sampled_indices, , drop = FALSE]
    
    # case series
    final_time[(n_b + 1):total_rows] <- time[case_indices]
    final_event[(n_b + 1):total_rows] <- status[case_indices]
    final_covs[(n_b + 1):total_rows, ] <- cov_matrix[case_indices, , drop = FALSE]
    
    # Pre-allocated and filled list
    list(
        time = final_time,
        event = final_event,
        covariates = final_covs,
        offset = rep(offset, total_rows)
    )
}

#' Fit a Penalized Multinomial Model on Case-Base Data
#'
#' A wrapper function to fit a penalized model (Elastic-Net or SCAD) using
#' an optimizer from the `cdSCRIP` package.
#'
#' @param cb_data A case-base dataset from `create_cb_data`.
#' @param regularization A string, either 'elastic-net' or 'SCAD'.
#' @param lambda The primary shrinkage parameter.
#' @param alpha The mixing/shape parameter. See `prepare_penalty_params`.
#' @param unpen_cov Integer. The number of leading covariates to leave unpenalized.
#' @param fit_fun The fitting function from `cdSCRIP` to use.
#' @param param_start Optional starting values for the coefficients.
#' @param standardize Logical. If TRUE, covariates are scaled to have mean 0 and SD 1.
#' @return The fitted model object from the specified `fit_fun`.
#' @export
fit_cb_model <- function(cb_data,
                         regularization = c('elastic-net', 'SCAD'), 
                         lambda, alpha = NULL,
                         fit_fun = MNlogisticAcc,
                         param_start = NULL,
                         n_unpenalized = 2,
                         standardize = TRUE, 
                         all_event_levels = NULL,
                         ...) {
    
    regularization <- rlang::arg_match(regularization)
    penalty_params <- prepare_penalty_params(regularization, lambda, alpha)
    
    if (is.null(all_event_levels)) all_event_levels <- sort(unique(cb_data$event))
    Y_factor <- factor(cb_data$event, levels = all_event_levels)
    
    penalized_covs <- cb_data$covariates
    scaler <- NULL
    if (standardize) {
        #  center and scale
        center <- colMeans(penalized_covs, na.rm = TRUE)
        scale <- apply(penalized_covs, 2, sd, na.rm = TRUE)
        
        # Check for zero variance
        if (any(scale == 0, na.rm = TRUE)) {
            cli::cli_warn("Some covariates have zero variance and were not scaled.")
            scale[scale == 0] <- 1
        }
        penalized_covs <- scale(penalized_covs, center = center, scale = scale)
        scaler <- list(center = center, scale = scale)
        
        if(!is.null(param_start)){
            
            original_cov_coefs <- param_start[1:ncol(penalized_covs), , drop = FALSE]
            
            intercept_adjustment <- crossprod(scaler$center, original_cov_coefs)
            param_start[ncol(penalized_covs) + 1, ] <- param_start[ncol(penalized_covs) + 1, ] + intercept_adjustment
            
            param_start_scaled <- sweep(original_cov_coefs, 1, scaler$scale, FUN = "*")
            
            param_start[1:ncol(penalized_covs), ] <- param_start_scaled
        }
        
    }
    
    X <- stats::model.matrix(~., 
                             data = data.frame(cbind(penalized_covs, 
                                                     time = log(cb_data$time))))
    
    X <- cbind(X[,-1], X[,1])
    
    # design matrix
    # X <- as.matrix(cbind(penalized_covs, time = log(cb_data$time), 1))
    
    opt_args <- list(X = X, Y = Y_factor, offset = cb_data$offset,
                     N_covariates = n_unpenalized, 
                     regularization = regularization,
                     transpose = FALSE, lambda1 = penalty_params$lambda1,
                     lambda2 = penalty_params$lambda2, lambda3 = 0,
                     param_start = param_start,
                     ...)
    
    fit <- do.call(fit_fun, opt_args)
    
    # Re-scaling
    if (standardize) {
        # Re-scale both matrices
        for (coef_type in c("coefficients", "coefficients_sparse")) {
            if (!is.null(fit[[coef_type]])) {
                coefs_scaled <- fit[[coef_type]]
                p_penalized <- ncol(penalized_covs)
                
                # Separate penalized from unpenalized coefficients
                beta_penalized_scaled <- coefs_scaled[1:p_penalized, , drop = FALSE]
                
                # Re-scale coefficients
                beta_penalized_orig <- sweep(beta_penalized_scaled, 1, scaler$scale, FUN = "/")
                
                # Adjust the intercept
                intercept_adjustment <- colSums(sweep(beta_penalized_scaled, 1, scaler$center / scaler$scale, FUN = "*"), na.rm = TRUE)
                intercept_orig <- coefs_scaled[nrow(coefs_scaled), ] - intercept_adjustment
                
                # full coefficient matrix on the original scale
                coefs_orig <- rbind(
                    beta_penalized_orig,
                    coefs_scaled[p_penalized + 1, , drop = FALSE], 
                    intercept_orig
                )
                
                # Update row names to match the original structure
                rownames(coefs_orig) <- c(colnames(cb_data$covariates), 
                                          "log(time)", "(Intercept)")
                
                # Replace the coefficients in the fit object
                fit[[coef_type]] <- round(coefs_orig, 8)
            }
        }
        return(c(fit, scaler = list(scaler)))
    }
    return(fit)
}


#################################################################
##                         Project: CV                         ##
#################################################################


#' Calculate Multinomial Deviance for a Fitted Case-Base Model
#' @param cb_data A case-base dataset from `create_cb_data`.
#' @param fit_object The output from `fit_cb_model`.
#' @param all_event_levels A vector of all possible event levels to ensure
#'   dimensional consistency.
#' @return The multinomial deviance value.
calc_multinom_deviance <- function(cb_data, fit_object, all_event_levels) {
    # Reconstruct the design matrix exactly as in fit_cb_model
    X <- as.matrix(cbind(cb_data$covariates, time = log(cb_data$time), 1))
    
    fitted_vals <- X %*% fit_object$coefficients
    if (is.vector(fitted_vals)) fitted_vals <- as.matrix(fitted_vals)
    
    pred_mat <- VGAM::multilogitlink(fitted_vals, inverse = TRUE)
    
    Y_fct <- factor(cb_data$event, levels = all_event_levels)
    
    Y_mat <- matrix(0, ncol = length(all_event_levels), nrow = nrow(X))
    
    valid_indices <- !is.na(Y_fct)
    Y_mat[cbind(which(valid_indices), as.integer(Y_fct)[valid_indices])] <- 1
    
    VGAM::multinomial()@deviance(mu = pred_mat, y = Y_mat, w = rep(1, nrow(X)))
}



#' Find the Maximum Lambda (lambda_max)
#' @param ... Other arguments passed to find_lambda_max.
#' @return The estimated lambda_max value.
find_lambda_max <- function(cb_data,
                            n_unpenalized,
                            alpha,
                            regularization,
                            ...) {
    
    n_event_types <- length(unique(cb_data$event))
    null_model_coefs <- n_unpenalized * (n_event_types - 1)
    search_grid <- round(exp(seq(log(7), log(0.005), length.out = 5)), 6)
    
    progressr::handlers("cli")
    
    with_progress({
        fine_grid_size <- 5
        
        p <- progressr::progressor(steps = length(search_grid) + fine_grid_size)
        
        # Helper function to fit the model and advance the progress bar
        fit_and_get_count <- function(lambda_val) {
            fit_model <- fit_cb_model(
                cb_data,
                lambda = lambda_val,
                regularization = regularization,
                alpha = alpha,
                n_unpenalized = n_unpenalized,
                ...
            )
            result <- sum(!same(fit_model$coefficients_sparse, 0))
            
            p()
            
            return(result)
        }
        
        # Coarse Search
        cli::cli_alert_info("Searching for lambda_max (coarse grid)...")
        coarse_results <- furrr::future_map_dbl(
            .x = search_grid,
            .f = fit_and_get_count,
            .options = furrr::furrr_options(seed = TRUE)
        )
        
        upper_idx <- which(coarse_results <= null_model_coefs)
        lower_idx <- which(coarse_results > null_model_coefs)
        
        # grid is too narrow.
        if (length(upper_idx) == 0 || length(lower_idx) == 0) {
            cli::cli_warn("Could not bracket lambda_max with the initial grid. Using largest value.")
            
            p(steps = fine_grid_size)
            
            return(max(search_grid))
        }
        
        #  bounds for the finer search.
        upper_bound <- min(search_grid[upper_idx])
        lower_bound <- max(search_grid[lower_idx])
        
        # Finer Search
        fine_grid <- round(seq(lower_bound, upper_bound, length.out = fine_grid_size), 6)
        cli::cli_alert_info("Searching for lambda_max (fine grid)...")
        fine_results <- furrr::future_map_dbl(
            .x = fine_grid,
            .f = fit_and_get_count, 
            .options = furrr::furrr_options(seed = TRUE)
        )
        

        # This is our best estimate for lambda_max.
        first_null_model_idx <- which(fine_results <= null_model_coefs)[1]
        lambda_max <- fine_grid[first_null_model_idx]
        
        # Fallback if the fine search fails for some reason.
        if (is.na(lambda_max)) {
            cli::cli_warn("Fine grid search failed to find a suitable lambda_max. Returning upper bound.")
            return(upper_bound)
        }
    }) 
    
    lambda_max<- lambda_max*1.2
    
    cli::cli_alert_success("Found lambda_max: {round(lambda_max, 4)}")
    
    return(lambda_max)
}

#' Generate a Lambda Grid for Regularization
#' @return A numeric vector of lambda values.
create_lambda_grid <- function(cb_data,
                               lambda,
                               nlambda,
                               lambda_max,
                               lambda.min.ratio,
                               regularization,
                               alpha,
                               n_unpenalized,
                               ...){
    
    nobs <- length(cb_data$event)
    
    nvars <- ncol(cb_data$covariates)
    
    if (is.null(lambda.min.ratio)) lambda.min.ratio <- ifelse(nobs < nvars, 0.01, 5e-04)
    
    if(is.null(lambda)){
        
        penalty_params_notif(regularization = regularization,
                             alpha = alpha)
        
        if(is.null(lambda_max)){
            lambda_max <- find_lambda_max(
                cb_data,
                regularization = regularization,
                alpha = alpha,
                n_unpenalized = n_unpenalized,
                ...)
            
            grid <- rev(exp(seq(log(lambda_max * lambda.min.ratio), 
                                log(lambda_max), length.out = nlambda)))
            
        }} else {
            
            grid <- sort(unique(lambda[lambda >= 0]), decreasing = TRUE) 
            
            if(length(grid) == 0) cli::cli_abort("Provided lambda values invalid.")
            
            
        }
    
            cli::cli_alert_info("Using {length(grid)} lambdas. Range: {signif(min(grid), 3)} to {signif(max(grid), 3)}")
    
    return(round(grid, 5))
}

#' Run a Single Fold of Cross-Validation
#' @return A vector of multinomial deviances for the fold.
run_cv_fold <- function(fold_indices, cb_data, 
                        lambdagrid, 
                        all_event_levels,
                        regularization,
                        alpha,
                        n_unpenalized = 2,
                        warm_start = T,
                        update_f = NULL,
                        ...) {
    # Split data into training and validation folds
    train_cv_data <- lapply(cb_data, function(x) if(is.matrix(x)) x[-fold_indices, , drop = FALSE] else x[-fold_indices])
    test_cv_data  <- lapply(cb_data, function(x) if(is.matrix(x)) x[fold_indices, , drop = FALSE] else x[fold_indices])
    
    # Initialize for loop
    deviances <- numeric(length(lambdagrid))
    non_zero <- numeric(length(lambdagrid))
    
    param_start <- NULL 
    
    for (i in seq_along(lambdagrid)) {
    
        model_info <- fit_cb_model(
            train_cv_data,
            lambda = lambdagrid[i],
            all_event_levels = all_event_levels,
            regularization = regularization,
            alpha = alpha,
            param_start = param_start, # Pass warm start
            n_unpenalized = n_unpenalized,
            ...
        )
        
        # Calculate deviance on the original, unscaled test fold
        deviances[i] <- calc_multinom_deviance(
            test_cv_data,
            model_info, # Contains re-scaled coefficients
            all_event_levels = all_event_levels
        )
        
        non_zero[i] <- sum(!same(model_info$coefficients_sparse, 0))
        
        update_f()
        
        # Update the warm start for the next iteration
        if(warm_start) param_start <- model_info$coefficients
    }
    
    return(list(deviances = deviances,
                non_zero = non_zero))
}


#' Cross-Validation for Penalized Multinomial Case-Base Models
#' @return An object of class `cb.cv` containing the results.
#' @export
cv_cbSCRIP <- function(formula, data, regularization = 'elastic-net',
                       cb_data = NULL,
                       alpha = NULL,
                       lambda = NULL,
                       nfold = 5, nlambda = 50, 
                       ncores = parallel::detectCores() / 2,
                       n_unpenalized = 2,
                       ratio = 25, ratio_event = 1,
                       lambda_max = NULL,
                       lambda.min.ratio = NULL,
                       warm_start = T,
                       ...) {
    
    if(is.null(cb_data)){
        cb_data <- create_cb_data(formula, data, ratio = ratio,
                                  ratio_event = ratio_event)
    }
    
    all_event_levels <- sort(unique(cb_data$event))
    
    n_event_types <- length(all_event_levels)
  
    
    if (ncores > 1) {
        future::plan(multisession, workers = ncores)
    } else {
        future::plan(sequential)
    }
    
    on.exit(future::plan(future::sequential), add = TRUE)
        
    lambdagrid <- create_lambda_grid(cb_data = cb_data,
                                     lambda = lambda,
                                     nlambda = nlambda,
                                     lambda_max = lambda_max,
                                     lambda.min.ratio = lambda.min.ratio,
                                     regularization = regularization,
                                     alpha = alpha,
                                     n_unpenalized = n_unpenalized,
                                     ...)
    
    folds <- caret::createFolds(factor(cb_data$event), k = nfold, list = TRUE)
    cli::cli_alert_info("Starting {nfold}-fold cross-validation...")
    
    progressr::handlers("cli")
    
    with_progress({
        p <- progressr::progressor(steps = nfold * nlambda)
        
        fold_list <- furrr::future_map(
            .x = folds, 
            .f = function(fold) {
                res <- run_cv_fold(
                    fold_indices = fold,
                    cb_data = cb_data,
                    lambdagrid = lambdagrid,
                    n_unpenalized = n_unpenalized,
                    regularization = regularization,
                    alpha = alpha,
                    all_event_levels = all_event_levels,
                    update_f = p,
                    ...
                )
                return(res)
            },
            .options = furrr::furrr_options(seed = TRUE)
        )
    })
    
    deviance_matrix <- lapply(fold_list, function(.x){.x$deviances})
    
    coeffs_list <- lapply(fold_list, function(.x){.x$coeffs})
    
    deviance_matrix <- do.call(cbind, deviance_matrix)
    
    # rownames(deviance_matrix) <- lambdagrid
    
    non_zero_matrix <- lapply(fold_list, function(.x){.x$non_zero})
    
    non_zero_matrix <- do.call(cbind, non_zero_matrix)
    
    mean_dev <- rowMeans(deviance_matrix)
    
    se_dev <- apply(deviance_matrix, 1, sd) / sqrt(nfold)
    
    lambda.min <- lambdagrid[which.min(mean_dev)]
    
    min_dev_upper_bound <- min(mean_dev) + se_dev[which.min(mean_dev)]
    
    lambda.1se <- max(lambdagrid[mean_dev <= min_dev_upper_bound])
    
    fit.min <- fit_cb_model(
        cb_data,
        regularization = regularization,
        lambda = lambda.min,
        alpha = alpha,
        n_unpenalized = n_unpenalized,
        ...
    )
    
    result <- list(
        fit.min = fit.min,
        lambdagrid = lambdagrid,
        deviance_matrix = deviance_matrix,
        non_zero_matrix = non_zero_matrix,
        deviance_mean = mean_dev,
        deviance_se = se_dev,
        lambda.min = lambda.min,
        lambda.1se = lambda.1se,
        cb_data = cb_data,
        call = match.call()
    )
    
    class(result) <- "cbSCRIP.cv"
    cli::cli_alert_success("Cross-validation complete.")
    return(result)
}


#' Fit a Penalized Multinomial Model over a Regularization Path
#'
#' Fits the case-base model for a sequence of lambda values using warm starts,
#' returning the path of coefficients. Includes a progress bar.
#'
#' @inheritParams cv_cb_model
#' @return An object of class `cb.path` containing coefficient paths.
#' @export
cbSCRIP <- function(formula, data, regularization = 'elastic-net', 
                    cb_data = NULL,
                    alpha = NULL,
                    lambda = NULL,
                    nlambda = 100, n_unpenalized = 2,
                    ratio = 25, ratio_event = 1,
                    lambda_max = NULL,
                    lambda.min.ratio = NULL,
                    warm_start = T,
                    ...) {
    
    # Create the full case-base dataset 
    cli::cli_alert_info("Creating case-base dataset...")
    
    if(is.null(cb_data)){
        
    cb_data <- create_cb_data(formula, data, ratio = ratio,
                              ratio_event = ratio_event)
    }
    
    all_event_levels <- sort(unique(cb_data$event))
    
    n_event_types <- length(all_event_levels)
    
    lambdagrid <- create_lambda_grid(cb_data = cb_data,
                                     lambda = lambda,
                                     nlambda = nlambda,
                                     lambda_max = lambda_max,
                                     lambda.min.ratio = lambda.min.ratio,
                                     regularization = regularization,
                                     alpha = alpha,
                                     n_unpenalized = n_unpenalized,
                                     ...)
    nlambda <- length(lambdagrid)
    
    
    cli::cli_alert_info("Fitting model path for {nlambda} lambda values...")
    
    path_fits <- vector("list", nlambda)
    
    param_start <- NULL #
    
    cli::cli_progress_bar("Fitting Path", total = nlambda)
    for (i in seq_along(lambdagrid)) {
        
        model_info <- fit_cb_model(
            cb_data = cb_data,
            lambda = lambdagrid[i],
            regularization = regularization,
            alpha = alpha,
            n_unpenalized = n_unpenalized,
            param_start = param_start, # Pass warm start
            ...
        )
        path_fits[[i]] <- model_info
        
        # Update the warm start for the next iteration
        if(warm_start) {
            param_start <- model_info$coefficients
        }
        cli::cli_progress_update()
    }
    # Extract the re-scaled coefficients 
    coefficients <- purrr::map(path_fits, ~.x$coefficients)
    
    non_zero <- purrr::map(path_fits, ~sum(!same(.x$coefficients, 0)))
    
    result <- list(
        coefficients = coefficients,
        non_zero = non_zero,
        lambdagrid = lambdagrid,
        models_info = path_fits,
        cb_data = cb_data,
        call = match.call()
    )
    
    class(result) <- "cbSCRIP"
    cli::cli_alert_success("Path fitting complete.")
    return(result)
}

#################################################################
##                        Project: Plot                        ##
#################################################################
 


#' Plot Cross-Validation Results
#' @return A ggplot object.
#' @export
plot.cbSCRIP.cv <- function(x, ...) {
    
    mean_n_vars <- rowMeans(x$non_zero_matrix, na.rm = TRUE)
    
    plot_data <- data.frame(
        lambda = x$lambdagrid,
        mean_dev = x$deviance_mean,
        upper = x$deviance_mean + x$deviance_se,
        lower = x$deviance_mean - x$deviance_se,
        n_vars = round(mean_n_vars)
    )
    
    n_total <- nrow(plot_data)
    n_labels <- min(n_total, 10) # Show at most 10 labels
    # Select ~10 evenly spaced indices from the data
    label_indices <- round(seq(1, n_total, length.out = n_labels))
    axis_labels_data <- plot_data[label_indices, ]
    
    ggplot(plot_data, aes(x = lambda, y = mean_dev)) +
        geom_errorbar(aes(ymin = lower, ymax = upper), 
                      width = 0.05, 
                      color = "grey80") +
        geom_point(color = "#f94144",
                   alpha = 1) +
        geom_vline(aes(xintercept = x$lambda.min,
                       color = "Lambda.min",
                       linetype = "Lambda.min")) +
        geom_vline(aes(xintercept = x$lambda.1se,
                       color = "Lambda.1se",
                       linetype = "Lambda.1se")) +
        scale_color_manual(
            name = NULL,
            values = c("Lambda.min" = "#277DA1", "Lambda.1se" = "#264653")
        ) +
        scale_linetype_manual(
            name = NULL,
            values = c("Lambda.min" = "dashed", "Lambda.1se" = "dotted")
        ) +
        labs(
            x = "Lambda",
            y = "Multinomial Deviance",
            title = "Cross-Validation Performance",
            color = "",
            linetype = ""
        ) +
        scale_x_log10(
            sec.axis = sec_axis(
                trans = ~.,
                name = "Mean Number of Selected Variables",
                breaks = axis_labels_data$lambda,  # Position ticks at these lambda values
                labels = axis_labels_data$n_vars   # Use n_vars as the label text
            )
        ) +
        theme_minimal()
}


#' Plot Coefficient Paths from a cb.path Object
#'
#' S3 method to plot the regularization path of coefficients.
#'
#' @param x An object of class `cb.path`.
#' @param plot_intercept Logical. Whether to include the intercept in the plot.
#' @param ... Not used.
#'
#' @return A ggplot object.
#' @export
plot.cbSCRIP <- function(x, plot_intercept = FALSE, ...) {
    
    # Wrangle the list of coefficient matrices into a long-format tibble
    x$lambdagrid
    plot_data <- imap_dfr(x$coefficients, ~{
        .x %>%
            as.data.frame() %>%
            tibble::rownames_to_column("variable") %>%
            mutate(lambda = as.numeric(x$lambdagrid[.y]))
    }) %>%
        pivot_longer(
            cols = -c(variable, lambda),
            names_to = "event_type",
            values_to = "coefficient"
        )
    
    if (!plot_intercept) {
        plot_data <- filter(plot_data, variable != "(Intercept)")
    }
    
    ggplot(plot_data, aes(x = lambda, y = coefficient, group = variable, color = variable)) +
        geom_line(alpha = 0.8) +
        facet_wrap(~event_type, scales = "free_y") +
        theme_minimal() +
        guides(color = "none") + # Hide legend for clarity if many variables
        labs(
            x = "Lambda",
            y = "Coefficient Value",
            title = "Coefficient Regularization Paths",
            subtitle = "Each line represents a variable's coefficient as penalty increases"
        ) +
        scale_x_log10()
}


weibull_hazard <- Vectorize(function(gamma, lambda, t) {
    return(gamma * lambda * t^(gamma - 1))
})

#' Simulate Competing Risks Data from Cause-Specific Hazards
#'
#' This function generates competing risks survival data using the cause-specific
#' hazards (CSH) framework. It implements the inverse transform sampling method
#' described by Binder et al. (2009) assuming Weibull baseline hazards for each cause.
#'
#' @param p Integer, total number of covariates.
#' @param n Integer, number of subjects to simulate.
#' @param beta1 Numeric vector of length `p`, coefficients for cause 1.
#' @param beta2 Numeric vector of length `p`, coefficients for cause 2.
#' @param nblocks Integer, number of blocks for block-diagonal correlation.
#' @param cor_vals Numeric vector of length `nblocks`, correlation for each block.
#' @param num.true Integer, number of non-zero ("true") covariates.
#' @param lambda01 Numeric, the baseline rate parameter for the Weibull hazard of cause 1.
#' @param lambda02 Numeric, the baseline rate parameter for the Weibull hazard of cause 2.
#' @param gamma1 Numeric, the baseline shape parameter for the Weibull hazard of cause 1.
#' @param gamma2 Numeric, the baseline shape parameter for the Weibull hazard of cause 2.
#' @param max_time Numeric, the maximum follow-up time (administrative censoring).
#' @param noise_cor Numeric, the correlation for noise variables.
#' @param rate_cens Numeric, the rate parameter for the exponential censoring distribution.
#' @param min_time Numeric, the minimum possible event time.
#' @param exchangeable Logical, if TRUE, use an exchangeable correlation structure
#'   for true covariates instead of a block-diagonal one.
#'
#' @return A data.frame with `n` rows and `p+2` columns ('fstatus', 'ftime',
#'   and covariates X1...Xp).
#'
cause_hazards_sim <- function(p, n, beta1, beta2,
                              nblocks = 4, cor_vals = c(0.7, 0.4, 0.6, 0.5), num.true = 20,
                              lambda01 = 0.55, lambda02 = 0.10,
                              gamma1 = 1.5, gamma2 = 1.5, max_time = 1.5, noise_cor = 0.1,
                              rate_cens = 0.05, min_time = 1/365, exchangeable = FALSE) {
    
    
    if(length(beta1) != p || length(beta2) != p) stop("Length of beta1 and beta2 must match p.")
    if(!exchangeable && nblocks != length(cor_vals)) stop("Length of cor_vals must match nblocks.")
    
    # Covariate Generation
    if(isTRUE(exchangeable)) {
        # Exchangeable correlation structure
        mat <- matrix(noise_cor, nrow = p, ncol = p)
        cor_exchangeable <- 0.5
        mat[1:num.true, 1:num.true] <- cor_exchangeable
        diag(mat) <- 1
        X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = mat)
    } else {
        # Block-diagonal correlation structure
        vpb <- num.true / nblocks
        correlation_matrix <- matrix(noise_cor, nrow = p, ncol = p)
        for (i in 1:nblocks) {
            start_index <- (i - 1) * vpb + 1
            end_index <- i * vpb
            correlation_matrix[start_index:end_index, start_index:end_index] <- cor_vals[i]
        }
        diag(correlation_matrix) <- 1
        X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = correlation_matrix)
    }
    
    # vent Time Generation 
    # Calculate individual-specific rate parameters 
    lambda1_i <- as.vector(lambda01 * exp(X %*% beta1))
    lambda2_i <- as.vector(lambda02 * exp(X %*% beta2))
    
    # Define the root-finding function: F(t) - u = 0
    cdf_solver <- function(t, g1, l1, g2, l2, u) {
        H1 <- l1 * t^g1 # Cumulative hazard for cause 1
        H2 <- l2 * t^g2 # Cumulative hazard for cause 2
        return((1 - exp(-(H1 + H2))) - u)
    }
    
    # Generate a uniform random variable for each subject
    u <- stats::runif(n)
    
    # For each subject, find the event time 't' by solving cdf_solver for 0
    times <- sapply(1:n, function(i) {
        stats::uniroot(
            cdf_solver,
            interval = c(0, max_time * 2),
            extendInt = "upX",
            g1 = gamma1, l1 = lambda1_i[i],
            g2 = gamma2, l2 = lambda2_i[i],
            u = u[i]
        )$root
    })
    
    # At the generated event time, determine the cause based on relative hazards
    hazard1 <- gamma1 * lambda1_i * times^(gamma1 - 1)
    hazard2 <- gamma2 * lambda2_i * times^(gamma2 - 1)
    prob_cause1 <- hazard1 / (hazard1 + hazard2)
    
    # Handle cases where total hazard is zero
    prob_cause1[is.nan(prob_cause1)] <- 0
    
    event_type <- stats::rbinom(n = n, size = 1, prob = prob_cause1)
    c.ind <- ifelse(event_type == 1, 1, 2)
    
    # Generate censoring times from an exponential distribution
    cens_times <- stats::rexp(n = n, rate = rate_cens)
    
    # Apply censoring: if censoring time is earlier, status is 0
    c.ind[cens_times < times] <- 0
    times <- pmin(times, cens_times)
    
    # Apply administrative censoring and winsorize time
    c.ind[times >= max_time] <- 0
    times <- pmin(times, max_time)
    times[times < min_time] <- min_time
    
    sim.data <- data.frame(fstatus = c.ind, ftime = times)
    X_df <- as.data.frame(X)
    colnames(X_df) <- paste0("X", seq_len(p))
    sim.data <- cbind(sim.data, X_df)
    
    return(sim.data)
}

#' Simulate Competing Risks Data from a Mixture Model
#'
#' @description
#' This function generates competing risks survival data from a mixture model framework.
#' A subject is first assigned a latent cause of failure, and the event time is
#' then drawn from a cause-specific Weibull distribution.
#'
#' **Note:** This method is distinct from and does **not** necessarily produce data
#' that follows a proportional sub-distribution hazards (Fine & Gray) model.
#'
#' @param n Integer, number of subjects to simulate.
#' @param p Integer, total number of covariates.
#' @param beta1 Numeric vector of length `p`, coefficients for cause 1.
#' @param beta2 Numeric vector of length `p`, coefficients for cause 2.
#' @param num.true Integer, number of non-zero ("true") covariates.
#' @param mix_p Numeric (0-1), base probability for the mixture assignment.
#' @param cor_vals Numeric vector, correlation for each block in block-diagonal structure.
#' @param noise_cor Numeric, the correlation for noise variables.
#' @param nblocks Integer, number of blocks for block-diagonal correlation.
#' @param lambda1 Numeric, the baseline rate parameter for the Weibull distribution of cause 1.
#' @param rho1 Numeric, the baseline shape parameter for the Weibull distribution of cause 1.
#' @param lambda2 Numeric, the baseline rate parameter for the Weibull distribution of cause 2.
#' @param rho2 Numeric, the baseline shape parameter for the Weibull distribution of cause 2.
#' @param cens_max Numeric, the maximum time for the uniform censoring distribution.
#' @param max_time Numeric, the maximum follow-up time (administrative censoring).
#' @param min_time Numeric, the minimum possible event time.
#' @param exchangeable Logical, if TRUE, use an exchangeable correlation structure.
#'
#' @return A data.frame with `n` rows and `p+2` columns ('fstatus', 'ftime',
#'   and covariates X1...Xp).
#'
cause_subdist_sim <- function(n, p, beta1, beta2, num.true = 20, mix_p = 0.5,
                              cor_vals = c(0.7, 0.4, 0.6, 0.5), noise_cor = 0.1,
                              nblocks = 4, lambda1 = 1, rho1 = 4,
                              lambda2 = 0.8, rho2 = 10, cens_max = 1.5,
                              max_time = 1.5, min_time = 1/365, exchangeable = FALSE) {
    
    if(length(beta1) != p || length(beta2) != p) stop("Length of beta1 and beta2 must match p.")
    
    if(isTRUE(exchangeable)) {
        # Exchangeable correlation structure
        mat <- matrix(noise_cor, nrow = p, ncol = p)
        cor_exchangeable <- 0.5
        mat[1:num.true, 1:num.true] <- cor_exchangeable
        diag(mat) <- 1
        X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = mat)
    } else {
        # Block-diagonal correlation structure
        vpb <- num.true / nblocks
        correlation_matrix <- matrix(noise_cor, nrow = p, ncol = p)
        for (i in 1:nblocks) {
            start_index <- (i - 1) * vpb + 1
            end_index <- i * vpb
            correlation_matrix[start_index:end_index, start_index:end_index] <- cor_vals[i]
        }
        diag(correlation_matrix) <- 1
        X <- mvtnorm::rmvnorm(n, mean = rep(0, p), sigma = correlation_matrix)
    }
    
    eta1_prob <- X %*% beta1
    prob_not_cause1 <- (1 - mix_p)^exp(eta1_prob)
    prob_cause1 <- 1 - prob_not_cause1
    c.ind <- 1 + rbinom(n, 1, prob = prob_cause1) # 1 = cause 2, 2 = cause 1
    
    # To match description: beta1 affects event 1, beta2 affects event 2
    c.ind <- ifelse(c.ind == 1, 2, 1)
    
    
    ftime <- numeric(n)
    
    # Subjects assigned to cause 1
    is_cause1 <- which(c.ind == 1)
    n1 <- length(is_cause1)
    if (n1 > 0) {
        eta1_time <- X[is_cause1, ] %*% beta1
        u1 <- runif(n1)
        t1 <- (-log(u1) / (lambda1 * exp(eta1_time)))^(1 / rho1)
        ftime[is_cause1] <- t1
    }
    
    # Subjects assigned to cause 2
    is_cause2 <- which(c.ind == 2)
    n2 <- length(is_cause2)
    if (n2 > 0) {
        eta2_time <- X[is_cause2, ] %*% beta2
        u2 <- runif(n2)
        t2 <- (-log(u2) / (lambda2 * exp(eta2_time)))^(1 / rho2)
        ftime[is_cause2] <- t2
    }
    
    cens_times <- stats::runif(n, min = 0, max = cens_max)
    
    # Apply censoring
    fstatus <- c.ind # Start with original cause
    fstatus[cens_times < ftime] <- 0
    ftime <- pmin(ftime, cens_times)
    
    # Apply administrative censoring and winsorize
    fstatus[ftime >= max_time] <- 0
    ftime <- pmin(ftime, max_time)
    ftime[ftime < min_time] <- min_time
    
    sim.data <- data.frame(fstatus = fstatus, ftime = ftime)
    X_df <- as.data.frame(X)
    colnames(X_df) <- paste0("X", seq_len(p))
    sim.data <- cbind(sim.data, X_df)
    
    return(sim.data)
}

#' Generate Competing Risks Survival Data for Simulation Studies
#'
#' @description
#' This function generates complex competing risks data based on five distinct
#' settings described in the simulation study. It handles the creation of
#' coefficient vectors, covariate correlation structures, and calls the appropriate
#' underlying simulation engine (either Cause-Specific Hazards or a Mixture Model).
#'
#' The five settings are:
#' 1.  **CSH: Single effects on endpoint 1.**
#' 2.  **CSH: Single effects on both endpoints (block structure).**
#' 3.  **CSH: Opposing effects.**
#' 4.  **CSH: Mixture of single and opposing effects.**
#' 5.  **Mixture Model: Opposing effects (violates CSH proportionality).**
#'
#' @param n Integer, total number of subjects to simulate.
#' @param p Integer, total number of covariates.
#' @param num_true Integer, number of non-zero ("true") covariates.
#' @param setting Integer (1-5), the simulation setting to use.
#' @param iter Integer, the seed for the simulation run for reproducibility.
#' @param sims Integer, optional, the total number of simulations for display purposes.
#'
#' @return A list containing:
#' \item{train}{A data.frame for the training set (75% of data).}
#' \item{test}{A data.frame for the test set (25% of data).}
#' \item{beta1}{The true coefficient vector for cause 1.}
#' \item{beta2}{The true coefficient vector for cause 2.}
#' \item{call}{The function call.}
#' \item{cen.prop}{The proportion of observations for each status (0=censored).}
#'
gen_data <- function(n = 400, p = 300,
                     num_true = 20, setting = 1,
                     iter = runif(1, 0, 9e5), sims = NULL) {
    
    cli::cli_alert_info("Setting: {setting} | Iteration {iter}/{sims} | p = {p} | k = {num_true}")
    set.seed(iter)
    set.seed(sample.int(5))
    
    beta1 <- rep(0, p)
    beta2 <- rep(0, p)
    nu_ind <- seq_len(num_true)
    k <- num_true
    
    # Define coefficient patterns based on the setting
    if (setting == 1) {
        beta1[nu_ind] <- 1
        beta2[nu_ind] <- 0
    } else if (setting == 2) {
        beta1[nu_ind] <- rep(c(1, 0, 1, 0), each = k / 4)
        beta2[nu_ind] <- rep(c(0, 1, 0, 1), each = k / 4)
    } else if (setting == 3) {
        beta1[nu_ind] <- rep(c(0.5, -0.5), times = k / 2)
        beta2[nu_ind] <- rep(c(-0.5, 0.5), times = k / 2)
    } else if (setting == 4) {
        beta1_true <- c(rep(1, k / 4),
                        rep(c(0.5, -0.5), times = k / 8),
                        rep(1, k / 4),
                        rep(0, k / 4))
        beta2_true <- c(rep(0, k / 4),
                        rep(c(-0.5, 0.5), times = k / 8),
                        rep(0, k / 4),
                        rep(1, k / 4))
        beta1[nu_ind] <- beta1_true
        beta2[nu_ind] <- beta2_true
    } else if (setting == 5) {
        beta1[nu_ind] <- 1
        beta2[nu_ind] <- -1
    } else {
        stop("'setting' must be an integer between 1 and 5.")
    }
    
    # Data Simulation 
    # Correctly choose simulation function and correlation structure based on 
    if (setting %in% c(1, 2, 3, 4)) {
        # CSH framework for settings 1-4
        sim.data <- cause_hazards_sim(
            n = n, p = p,
            beta1 = beta1, beta2 = beta2,
            num.true = k,
            exchangeable = (setting == 1), # Exchangeable for setting 1
            lambda01 = 0.55, lambda02 = 0.35,
            gamma1 = 1.5, gamma2 = 1.5
        )
    } else if (setting == 5) {
        # Mixture Model framework for setting 5
        sim.data <- cause_subdist_sim(
            n = n, p = p,
            beta1 = beta1, beta2 = beta2,
            num.true = k,
            exchangeable = TRUE, # Exchangeable for setting 5
            cens_max = 1.5
        )
    }
    
    # Train-Test Split 
    train.index <- caret::createDataPartition(sim.data$fstatus, p = 0.75, list = FALSE)
    train <- sim.data[train.index, ]
    test <- sim.data[-train.index, ]
    
    return(list(
        train = train,
        test = test,
        beta1 = beta1,
        beta2 = beta2,
        call = match.call(),
        cen.prop = prop.table(table(factor(sim.data$fstatus, levels = 0:2)))
    ))
}

#################################################################
##                        Project: loss                        ##
#################################################################
 

#' Calculate Multinomial Negative Log-Likelihood (Loss)
#'
#' This function computes the total negative log-likelihood for a multinomial
#' logistic regression model with a baseline class.
#'
#' @param coefficients A p x K matrix of model coefficients.
#' @param X A n x p design matrix of predictors.
#' @param Y A numeric vector of length 'n' containing the true class labels.
#'   The baseline class should be coded as 0.
#' @param offset A numeric vector of length 'n' for the offset term.
#' @return The total negative log-likelihood (a single numeric value).
calculate_multinomial_loss <- function(coefficients, X, Y, offset = NULL) {
    # Calculate linear scores (eta)
    eta <- X %*% coefficients
    if (!is.null(offset)) {
        eta <- eta + offset
    }
    
    #  log-sum-exp trick
    eta_with_baseline <- cbind(eta, 0)
    max_scores <- apply(eta_with_baseline, 1, max)
    log_denominators <- max_scores + log(rowSums(exp(sweep(eta_with_baseline, 1, max_scores, "-"))))
    
    # Get the linear score
    n <- nrow(X)
    numerators <- numeric(n)
    rows_with_event <- which(Y > 0)
    
    if (length(rows_with_event) > 0) {
        true_class_indices <- Y[rows_with_event]
        numerators[rows_with_event] <- eta[cbind(rows_with_event, true_class_indices)]
    }
    
    total_loss <- -sum(numerators - log_denominators)
    
    return(total_loss/n)
}


#' Calculate Penalized Multinomial Negative Log-Likelihood
#'
#' Adds the Elastic Net penalty to the negative log-likelihood loss.
#'
#' @param coefficients A p x K matrix of model coefficients. Assumes the first
#'   row corresponds to the non-penalized intercept.
#' @param X A n x p design matrix.
#' @param Y A numeric vector of length n with class labels (0 for baseline).
#' @param offset A numeric vector of length n for the offset term.
#' @param lambda A single numeric value for the overall regularization strength.
#' @param alpha The Elastic Net mixing parameter (1 for Lasso, 0 for Ridge).
#' @param penalty_weights A numeric vector of weights for each penalized covariate,
#'   of length p-1. Defaults to 1 for all variables.
#' @return The total penalized negative log-likelihood.
calculate_penalized_multinomial_loss <- function(
        coefficients, X, Y, offset = NULL,
        lambda, alpha, penalty_weights = NULL
) {
    
    # Calculate the Negative Log-Likelihood 
    nll_loss <- calculate_multinomial_loss(coefficients, X, Y, offset)
    
    # Calculate the Regularization Penalty
    
    # Isolate coefficients to be penalized
    penalized_coefs <- coefficients
    
    # default penalty
    if (is.null(penalty_weights)) {
        penalty_weights <- rep(1, nrow(penalized_coefs))
    }
    
    if (length(penalty_weights) != nrow(penalized_coefs)) {
        stop("Length of 'penalty_weights' must match the number of penalized covariates.")
    }
    
    l1_component <- sum(penalty_weights * rowSums(abs(penalized_coefs)))
    l2_component <- sum(penalty_weights * rowSums(penalized_coefs^2))
    
    en_penalty <- lambda * (alpha * l1_component + (1 - alpha) / 2 * l2_component)
    
    total_loss <- nll_loss + en_penalty
    
    return(total_loss)
}