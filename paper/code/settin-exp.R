# library(casebase)
# library(future.apply)
# library(glmnet)
# library(mtool)
# library(parallel)
# library(tictoc)
# library(tidyverse)
# library(foreach)
# library(survival)
# library(cmprsk)
# library(glue)
# library(pec)
# library(survminer)
# library(here)
pacman::p_load(casebase, 
               glmnet, 
               parallel, 
               tictoc, 
               dplyr, tidyr,
               foreach,
               survival,
               cmprsk,
               here,
               caret,
               # survminer,
               pec,
               glue)

# Fitting functions 
source(here("paper",
            "code", "fitting_functionsV2.R"))

# sims <- 100
# setting <- 1
# p <- 120
# num_true <- round(p*0.08333333)*2
if(!exists("start_sims")) start_sims <- 1

for (i in start_sims:sims) {
    
    # Setup
    n <- 400
    data <- gen_data(n = n, 
                     p = p, 
                     num_true = num_true,
                     setting = setting,
                     iter = i, 
                     sims = sims)
    beta1 <- data$beta1
    beta2 <- data$beta2
    train <- data$train
    test <- data$test
    
    ##############################################################
    # We have two competitor models for variable selection:
    # 1) Independent cox-regression model 
    # 2) penCR cox regression model - where the lambda penalties are trained together 
    ######################## Fit indepedent cox-regression model ###############################
    ######################### Cause-1 #########################################
    # Censor competing event
    y_train <- Surv(time = train$ftime, event = train$fstatus == 1)
    
    x_train <- model.matrix(~ . -ftime -fstatus, data = train)[, -1] 
    
    # Censor competing event
    y_test <- Surv(time = test$ftime, event = test$fstatus == 1)
    
    x_test <- model.matrix(~ . -ftime -fstatus, data = test)[, -1] 
    
    # Fit cause-specific cox model with glmnet on training set 
    cox_mod <- cv.glmnet(x = x_train, y = y_train, family = "cox", alpha = 0.7)
    
    # Fit on validation set 
    cox_val_min <- glmnet(x = x_test, y = y_test, family = "cox", alpha = 0.7, 
                          lambda = cox_mod$lambda.min)
    
    cc_min <- coef(cox_val_min)
    
    res_cox_min1 <- varsel_perc(cc_min, beta1)
    
    # let's calculate the bias for all the competitors as well (a task could be turning this one line into a function as well)
    # Only for the true non-zero variables
    cox_mse1 <- mse_bias(cc_min, beta1)
    
    cli::cli_alert_success("Cox-regression model ran")
    
    ########################## Fit PenCR model ##################################
    penCR = cv.glmnet.CR(data = train, family="cox", alpha= 0.7, standardize= TRUE,
                         nlambda = 20, t.BS = median(train$ftime), seed = 115, causeOfInt = 1,
                         nfold = 5)
    
    cc_min_penCR1 <- penCR$glmnet.fits$models$`Cause 1`$glmnet.res$lambda[penCR$min.index[1]]
    cc_min_penCR2 <- penCR$glmnet.fits$models$`Cause 2`$glmnet.res$lambda[penCR$min.index[2]]
    
    # Fit on validation set 
    penCR_val_min1 <- glmnet(x = x_test, y = y_test, family = "cox", alpha = 0.7, 
                             lambda = cc_min_penCR1)
    
    cc_min_penCR1 <- coef(penCR_val_min1)
    
    penCR_val_min2 <- glmnet(x = x_test, y = y_test, family = "cox", alpha = 0.7, 
                             lambda = cc_min_penCR2)
    
    res_pencr_min1 <- varsel_perc(cc_min_penCR1, beta1)
    
    
    # Calculate MSE here as well (try and fill it out!)
    penCR_mse1 <-mse_bias(cc_min_penCR1, beta1)
    
    cli::cli_alert_success("PenCR model ran")
    ########################## Fit casebase model #################################
    # Train case-base model through cross-validation 
    
    # tic()
    # cv.lambda <- mtool.multinom.cv(train, seed = 1, alpha = 0.7, nfold = 5)
    # toc()
    # 
    # # Test set 
    # surv_obj_val <- with(test, Surv(ftime, as.numeric(fstatus), type = "mstate"))
    # 
    # # Covariance matrix
    # cov_val <- cbind(test[, c(grepl("X", colnames(test)))], time = log(test$ftime))
    # 
    # # Case-base dataset
    # cb_data_val <- create_cbDataset(surv_obj_val, as.matrix(cov_val), ratio = 10)
    # 
    # # Case-base fits 
    # # Lambda.min
    # fit_val_min <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
    #                            lambda = cv.lambda$lambda.min, alpha = 0.7,
    #                            unpen_cov = 2)
    # 
    # 
    # res_cb_min1 <- varsel_perc(fit_val_min$coefficients[1:eval(parse(text="p")), 1], beta1)
    # 
    # res_cb_min2 <- tryCatch(
    #     {varsel_perc(fit_val_min$coefficients[1:eval(parse(text="p")), 2], beta2)}, 
    #     error = function(msg){
    #         return(list(TP = NA, TN = NA, FP = NA, FN = NA,
    #                     Sensitivity = NA, Specificity = NA,
    #                     MCC = NA))
    #     })
    # 
    # # Calculate MSE here as well
    # casebase_mse1 <-mse_bias(fit_val_min$coefficients[1:eval(parse(text="p")), 1], beta1)
    # ############################ Case-base with post-"elastic net" ###########################
    # cb_postenet <- tryCatch(
    #     {multinom.post_enet(fit_val_min, cause = 1)
    #     }, error = function(msg){
    #         a <- list()
    #         a$coefs_all <- NA
    #         return(a)
    #     })
    # 
    # # Calculate MSE here as well
    # casebase_mse_enet1 <-mse_bias(coef = cb_postenet$coefs_all, true_coefs = beta1)
    
    ########################## Fit casebase Acc model #################################
    # Train case-base model through cross-validation 
    
    mtool_cv <- purrr::partial(mtool::mtool.MNlogisticAcc,
                             niter_inner_mtplyr = 1.5,
                             maxit = 100,
                             momentum_gamma = 0,
                             tolerance = 1e-8,
                             learning_rate = 1e-4,
                             verbose = F
                             )
    
    # mtool_cv <- function(..., niter_inner_mtplyr = 1.5,
    #                          maxit = 100,
    #                          momentum_gamma = 0.3,
    #                          tolerance = 5e-3,
    #                          learning_rate = 1e-4,
    #                          verbose = F) {
    #     args <- list(...,
    #                  niter_inner_mtplyr = niter_inner_mtplyr,
    #                  maxit = maxit,
    #                  momentum_gamma = momentum_gamma,
    #                  tolerance = tolerance,
    #                  verbose = verbose,
    #                  learning_rate = learning_rate
    #                  )
    #     do.call(mtool::mtool.MNlogisticAcc, args)
    # }
    
    
    tic()
    cv.lambdaAcc <- mtool.multinom.cv(train, seed = 1, alpha = 0.7, 
                                         nfold = 5, 
                                         fit_fun = mtool_cv,
                                         # ncores = parallelly::availableCores()/2,
                                         # train_ratio = 5,
                                         # lambda_max = 0.5,
                                         # ws = F,
                                         grid_size = 60,
                                      )
    toc()
    
    plot_cv.multinom(cv.lambdaAcc)
    
    
    # Test set 
    surv_obj_val <- with(test, Surv(ftime, as.numeric(fstatus), type = "mstate"))
    
    # Covariance matrix
    cov_val <- cbind(test[, c(grepl("X", colnames(test)))], time = log(test$ftime))
    
    # Case-base dataset
    cb_data_val <- create_cbDataset(surv_obj_val, as.matrix(cov_val))
    
    # Case-base fits 
    # Lambda.min
    
    cli::cli_alert_info("Lambda.min: {cv.lambdaAcc$lambda.min}")
    fit_val_min_acc <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                                   lambda = cv.lambdaAcc$lambda.min, alpha = 0.7,
                                   unpen_cov = 2,
                                   fit_fun = mtool_cv)
    
    res_cb_min1_acc <- varsel_perc(fit_val_min_acc$coefficients[1:eval(parse(text="p")), 1], beta1)
    
    res_cb_min2_acc <- tryCatch(
        {varsel_perc(fit_val_min_acc$coefficients[1:eval(parse(text="p")), 2], beta2)}, 
        error = function(msg){
            return(list(TP = NA, TN = NA, FP = NA, FN = NA,
                        Sensitivity = NA, Specificity = NA,
                        MCC = NA))
        })
    
    # Calculate MSE here as well
    casebase_mse1_acc <-mse_bias(fit_val_min_acc$coefficients[1:eval(parse(text="p")), 1], beta1)
    
    ############################ Case-base with post-"elastic net" ###########################
    cb_postenet_acc <- tryCatch(
        {multinom.post_enet(fit_val_min_acc, cause = 1)
        }, error = function(msg){
            a <- list()
            a$coefs_all <- NA
            return(a)
        })
    
    # Calculate MSE here as well
    casebase_mse_enet1_acc <- mse_bias(coef = cb_postenet_acc$coefs_all, true_coefs = beta1)
    
    
    cli::cli_alert_success("casebase model ran")
    
    
    ########################## Fit casebase SCAD model #################################
    # Train case-base model through cross-validation 
    
    mtool_cv <- purrr::partial(mtool::mtool.MNlogisticAcc,
                               niter_inner_mtplyr = 1.5,
                               maxit = 100,
                               momentum_gamma = 0,
                               tolerance = 1e-3,
                               learning_rate = 1e-4,
                               verbose = F
    )
    
    # mtool_cv <- function(..., niter_inner_mtplyr = 1.5,
    #                          maxit = 100,
    #                          momentum_gamma = 0.3,
    #                          tolerance = 5e-3,
    #                          learning_rate = 1e-4,
    #                          verbose = F) {
    #     args <- list(...,
    #                  niter_inner_mtplyr = niter_inner_mtplyr,
    #                  maxit = maxit,
    #                  momentum_gamma = momentum_gamma,
    #                  tolerance = tolerance,
    #                  verbose = verbose,
    #                  learning_rate = learning_rate
    #                  )
    #     do.call(mtool::mtool.MNlogisticAcc, args)
    # }
    
    
    tic()
    cv.lambdaSCAD <- mtool.multinom.cv(train, seed = 1, 
                                      nfold = 5, 
                                      fit_fun = mtool_cv,
                                      epsilon = 0.05,
                                      constant_covariates = 2,
                                      # ncores = parallelly::availableCores()/2,
                                      # train_ratio = 5,
                                      # lambda_max = 0.5,
                                      # ws = F,
                                      grid_size = 60,
                                      regularization = 'SCAD'
    )
    toc()
    
    plot_cv.multinom(cv.lambdaSCAD)
    
    
    # Test set 
    surv_obj_val <- with(test, Surv(ftime, as.numeric(fstatus), type = "mstate"))
    
    # Covariance matrix
    cov_val <- cbind(test[, c(grepl("X", colnames(test)))], time = log(test$ftime))
    
    # Case-base dataset
    cb_data_val <- create_cbDataset(surv_obj_val, as.matrix(cov_val))
    
    # Case-base fits 
    # Lambda.min
    
    cli::cli_alert_info("Lambda.min: {cv.lambdaSCAD$lambda.min}")
    fit_val_min_SCAD <- fit_cbmodel(cb_data_val, regularization = 'SCAD',
                                   lambda = cv.lambdaSCAD$lambda.min,
                                   unpen_cov = 2,
                                   fit_fun = mtool_cv)
    
    res_cb_min1_SCAD <- varsel_perc(fit_val_min_SCAD$coefficients[1:eval(parse(text="p")), 1], beta1)
    
    res_cb_min2_SCAD <- tryCatch(
        {varsel_perc(fit_val_min_SCAD$coefficients[1:eval(parse(text="p")), 2], beta2)}, 
        error = function(msg){
            return(list(TP = NA, TN = NA, FP = NA, FN = NA,
                        Sensitivity = NA, Specificity = NA,
                        MCC = NA))
        })
    
    # Calculate MSE here as well
    casebase_mse1_SCAD <- mse_bias(fit_val_min_SCAD$coefficients[1:eval(parse(text="p")), 1], beta1)
    
    
    
    cli::cli_alert_success("casebase SCAD model ran")
    
    
    ### Oracle ---
    true_vars <- which(beta1 != 0 | beta2 != 0)
    
    cb_data_val_true <- cb_data_val
    
    cb_data_val_true$covariates <- cb_data_val$covariates[,c(true_vars, p + 1)]
    
    fit_val_min_ora <- fit_cbmodel(cb_data_val_true, regularization = 'elastic-net',
                                   lambda = 1e-20, alpha = 0.5,
                                   unpen_cov = 2,
                                   fit_fun = mtool_cv)
    
    res_cb_min1_ora <- varsel_perc(fit_val_min_ora$coefficients[1:length(true_vars), 1], beta1)
    
    casebase_mse_enet1_ora <- mse_bias(coef = fit_val_min_ora$coefficients, true_coefs = beta1)
    
    
    

    
    res_sel <- rbind(res_cox_min1, res_pencr_min1, res_cb_min1_acc,
                     res_cb_min1_SCAD,
                     res_cb_min1_ora)
    
    rownames(res_sel) <- c("enet-iCR", "enet-penCR", 
                            "enet-casebase-Acc",
                           "SCAD-casebase",
                           "enet-casebase-oracle")
    
    ############################# Format and export results ###############################
    
    coef_all <- rbind(as.vector(cc_min), as.vector(cc_min_penCR1), 
                      # as.vector(fit_val_min$coefficients[1:eval(parse(text="p")), 1]), cb_postenet$coefs_all,
                      as.vector(fit_val_min_acc$coefficients[1:eval(parse(text="p")), 1]), 
                      cb_postenet_acc$coefs_all,
                      as.vector(fit_val_min_SCAD$coefficients[1:eval(parse(text="p")), 1]),
                      as.vector(fit_val_min_ora$coefficients[1:length(true_vars), 1]))
    
    colnames(coef_all) <- paste("X", 1:eval(parse(text="p")), sep = "")
    
    rownames(coef_all) <- c("enet-iCR", "enet-penCR", 
                            # "enet-casebase", "postCB-bias",
                            "enet-casebase-Acc", "postenet-casebase-Acc",
                            "SCAD-casebase",
                            "enet-casebase-oracle")
    
    Res_bias <- rbind(cox_mse1, penCR_mse1, #casebase_mse1, casebase_mse_enet1,
                      casebase_mse1_acc, casebase_mse_enet1_acc,
                      casebase_mse1_SCAD,
                      casebase_mse_enet1_ora)
    
    rownames(Res_bias) <- c("enet-iCR", "enet-penCR",
                            # "enet-casebase", "postCB-bias",
                            "enet-casebase-Acc", "postenet-casebase-Acc",
                            "SCAD-casebase",
                            "enet-casebase-oracle")
    
    # Models
    models <- list("cox_val_min" = cox_val_min,
                   "cox_val_min_cv" = cox_mod,
                   "penCR_val_min_cv" = penCR,
                   "penCR_val_min1" = penCR_val_min1,
                   "penCR_val_min2" = penCR_val_min2,
                   "fit_val_min_acc_cv" = cv.lambdaAcc,
                   "fit_val_min_acc" = fit_val_min_acc)
    
    saveRDS(models,
            here("paper",
                "results",
                glue("models_setting-{setting}_iter-{i}_p-{p}_k-{num_true}.rds")))
    
    # colnames(Res_bias) <- "Coefficient Bias"
    
    # write.csv(Res_bias, file = here("paper", "results",
    #                                 glue("bias_setting-{setting}_iter-{i}_p-{p}.csv")))
    
    print(res_sel)
    
    print(Res_bias)
    
    write.csv(coef_all, file = here("paper",
                                    "results",
                                    glue("coefficients_setting-{setting}_iter-{i}_p-{p}_k-{num_true}.csv")))
}


