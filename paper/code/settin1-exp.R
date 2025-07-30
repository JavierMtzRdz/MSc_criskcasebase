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

for (i in 1:sims) {
    
    # Set seed
    set.seed(i)
    the_seed <- as.numeric(paste0(round(sample(runif(10, 0,9), 5)), collapse = ""))
    cli::cli_alert_info("Setting: {setting} | Iteration {i}/{sims} | seed: {the_seed} | p = {p} | k = {num_true}")
    set.seed(the_seed)
    
    # Setup
    n <- 400
    # Let's start yours with p = 20
    beta1 <- c(rep(0, p))
    beta2 <- c(rep(0, p))
    nu_ind <- seq(num_true)
    # Here out of 20 predictors, 10 should be non-zero 
    beta1[nu_ind] <- c(rep(1, num_true/4), rep(1, num_true/4), rep(1, num_true/4), rep(1, num_true/4))
    beta2[nu_ind] <- c(rep(0, num_true/4), rep(0, num_true/4), rep(0, num_true/4), rep(0, num_true/4))
    
    
    
    # Simulate data
    sim.data <- cause_hazards_sim(n = n, p = p, 
                                  nblocks = 4, num.true = num_true, 
                                  beta1 = beta1, beta2 = beta2, rate_cens = 0.25, 
                                  h1 = 0.55, h2 = 0.35, gamma1 = 1.5, gamma2 = 1.5, exchangeable = TRUE)
    
    
    # Censoring proportion
    cen.prop <- c(prop.table(table(sim.data$fstatus)), 0, 0, 0, 0)
    
    
    # Let us plot and visualize the competing risks curves
    #cif <- cuminc(ftime = sim.data$ftime, fstatus = sim.data$fstatus)
    
    #ggcompetingrisks(cif)
    
    # Training-test split 
    # We only do this (instead of generating datasets for train and test like Anthony mentioned because it is faster computationally 
    # as casebase resamples) + proportion of censoring can be quite random in each run of the simulation so we want to maintain the same in validation and test set
    train.index <- caret::createDataPartition(sim.data$fstatus, p = 0.75, list = FALSE)
    train <- sim.data[train.index,]
    test <- sim.data[-train.index,]
    
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
    
    mtool3 <- purrr::partial(mtool::mtool.MNlogisticAcc,
                             niter_inner_mtplyr = 1,
                             maxit = 100,
                             momentum_gamma = .96,
                             tolerance = 1e-2,
                             learning_rate = 1e-4,
                             verbose = F)
    
    tic()
    cv.lambdaAcc <- mtool.multinom.cv.ws(train, seed = 1, alpha = 0.7, 
                                         nfold = 5, 
                                         fit_fun = mtool3)
    toc()
    
    # Test set 
    surv_obj_val <- with(test, Surv(ftime, as.numeric(fstatus), type = "mstate"))
    
    # Covariance matrix
    cov_val <- cbind(test[, c(grepl("X", colnames(test)))], time = log(test$ftime))
    
    # Case-base dataset
    cb_data_val <- create_cbDataset(surv_obj_val, as.matrix(cov_val), ratio = 10)
    
    # Case-base fits 
    # Lambda.min
    cli::cli_alert_info("Lambda.min: {cv.lambdaAcc$lambda.min}")
    
    fit_val_min_acc <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                    lambda = cv.lambdaAcc$lambda.min, alpha = 0.7,
                    unpen_cov = 2,
                    fit_fun = mtool3)
    
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
    casebase_mse_enet1_acc <-mse_bias(coef = cb_postenet_acc$coefs_all, true_coefs = beta1)
    
    cli::cli_alert_success("casebase model ran")
    res_cb_min1_acc
    
    res_sel <- rbind(res_cox_min1, res_pencr_min1, res_cb_min1_acc)
    
    rownames(res_sel) <- c("enet-iCR", "enet-penCR", 
                            "enet-casebase-Acc")
    
    ############################# Format and export results ###############################
    Res_bias <- rbind(cox_mse1, penCR_mse1, #casebase_mse1, casebase_mse_enet1,
                      casebase_mse1_acc, casebase_mse_enet1_acc)
    
    coef_all <- rbind(as.vector(cc_min), as.vector(cc_min_penCR1), 
                      # as.vector(fit_val_min$coefficients[1:eval(parse(text="p")), 1]), cb_postenet$coefs_all,
                      as.vector(fit_val_min_acc$coefficients[1:eval(parse(text="p")), 1]), cb_postenet_acc$coefs_all)
    
    colnames(coef_all) <- paste("X", 1:eval(parse(text="p")), sep = "")
    
    rownames(coef_all) <- c("enet-iCR", "enet-penCR", 
                            # "enet-casebase", "postCB-bias",
                            "enet-casebase-Acc", "postenet-casebase-Acc")
    
    rownames(Res_bias) <- c("enet-iCR", "enet-penCR",
                            # "enet-casebase", "postCB-bias",
                            "enet-casebase-Acc", "postenet-casebase-Acc")
    
    # Models
    models <- list("cox_val_min" = cox_val_min,
                   "penCR_val_min1" = penCR_val_min1,
                   "penCR_val_min2" = penCR_val_min2,
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


