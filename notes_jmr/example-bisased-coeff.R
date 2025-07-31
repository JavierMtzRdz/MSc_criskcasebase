library(casebase)
library(future.apply)
library(glmnet)
library(mtool)
library(parallel)
library(tictoc)
library(tidyverse)
library(foreach)
library(survival)
library(cmprsk)
library(glue)
library(pec)
library(knitr)
library(kableExtra)

# Fitting functions 
source(here::here("src/fitting_functionsV2.R"))

run_simulation <- function(seed) {
    # Set seed
    cli::cli_alert_info("Start simulation w/ seed {seed}")
    
    seed <- as.integer(Sys.time())
    
    # take the last five digits of the initial seed
    the_seed= seed %% 100000
    set.seed(the_seed)
    
    
    # Setup
    n <- 400
    p <- 1000
    num_true <- 24
    beta1 <- c(rep(0, p))
    beta2 <- c(rep(0, p))
    nu_ind <- seq(num_true)
    nu_ind <- seq(num_true)
    beta1[nu_ind] <- c(rep(1, 6), 0.5, -0.5, 0.5, -0.5, 0.5, -0.5, rep(1, 6), rep(0, 6))
    beta2[nu_ind] <- c(rep(0, 6), 0.5, -0.5, 0.5, -0.5, 0.5, -0.5, rep(0, 6), rep(1, 6))
    
    
    # Simulate data
    sim.data <- cause_hazards_sim(n = n, p = p, nblocks = 4, num.true = 24, 
                                  beta1 = beta1, beta2 = beta2, rate_cens = 0.05, 
                                  h1 = 0.55, h2 = 0.35, gamma1 = 1.5, gamma2 = 1.5)
    
    
    # Censoring proportion
    cen.prop <- c(prop.table(table(sim.data$fstatus)), 0, 0, 0, 0)
    
    # Training-test split 
    train.index <- caret::createDataPartition(sim.data$fstatus, p = 0.75, list = FALSE)
    train <- sim.data[train.index,]
    test <- sim.data[-train.index,]
    ######################## Fit cox-regression model ###############################
    ######################### Cause-1 #########################################
    # Censor competing event
    y_train1 <- Surv(time = train$ftime, event = train$fstatus == 1)
    
    x_train1 <- model.matrix(~ . -ftime -fstatus, data = train)[, -1] 
    
    # Censor competing event
    y_test1 <- Surv(time = test$ftime, event = test$fstatus == 1)
    
    x_test1 <- model.matrix(~ . -ftime -fstatus, data = test)[, -1] 
    
    tictoc::tic()
    # Fit cause-specific cox model with glmnet on training set 
    cox_mod1 <- cv.glmnet(x = x_train1, y = y_train1, family = "cox", alpha = 0.7)
    
    # Fit on validation set 
    cox_val_min1 <- glmnet(x = x_test1, y = y_test1, family = "cox", alpha = 0.7, 
                           lambda = cox_mod1$lambda.min)
    
    time_cox1 <- as.numeric(tictoc::toc()$toc - time_casebase$tic)
    
    cli::cli_alert_success("Fit cox-regression model cause 1 finished")
    
    cc_min1 <- coef(cox_val_min1)
    
    res_cox_min1 <- varsel_perc(cc_min1, beta1) %>% bind_cols(tibble(Time = time_cox1))
    
    # MSE
    mse_cox <- mean((cc_min1[1:24] - beta1[1:24])^2)
    
    # True coefs
    cc_true <- cc_min1[1:24]
    ########################## Cause 2 #####################################
    tictoc::tic()
    # Censor competing event
    y_train2 <- Surv(time = train$ftime, event = train$fstatus == 2)
    
    x_train2 <- model.matrix(~ . -ftime -fstatus, data = train)[, -1] 
    
    # Censor competing event
    y_test2 <- Surv(time = test$ftime, event = test$fstatus == 2)
    
    x_test2 <- model.matrix(~ . -ftime -fstatus, data = test)[, -1] 
    
    # Fit cause-specific cox model with glmnet on training set 
    cox_mod2 <- cv.glmnet(x = x_train2, y = y_train2, family = "cox", alpha = 0.7)
    
    # Fit on validation set 
    cox_val_min2 <- glmnet(x = x_test2, y = y_test2, family = "cox", alpha = 0.7, 
                           lambda = cox_mod2$lambda.min)
    
    time_cox2 <- as.numeric(tictoc::toc()$toc - time_casebase$tic)
    
    cli::cli_alert_success("Fit cox-regression model cause 2 finished")
    cc_min2 <- coef(cox_val_min2)
    
    res_cox_min2 <- varsel_perc(cc_min2, beta2) %>% bind_cols(tibble(Time = time_cox2))
    ########################## Fit casebase model #################################
    tictoc::tic()
    # Train case-base model 
    cv.lambda <- mtool.multinom.cv(train, seed = 1, alpha = 0.7, nfold = 5, train_ratio = 20)
    
    cv.lambda$lambda.1se
    
    # Test set 
    surv_obj_val <- with(test, Surv(ftime, as.numeric(fstatus), type = "mstate"))
    
    # Covariance matrix
    cov_val <- cbind(test[, c(grepl("X", colnames(test)))], time = log(test$ftime))
    
    # Case-base dataset
    cb_data_val <- create_cbDataset(surv_obj_val, as.matrix(cov_val), ratio = 10)
    
    time_casebase_cv <- as.numeric(tictoc::toc()$toc - time_casebase$tic)
    
    # Case-base fits 
    # Lambda.min
    tictoc::tic()
    fit_val_min <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                               lambda = cv.lambda$lambda.min , alpha = 0.7, unpen_cov = 2)
    
    time_casebase_min <- as.numeric(tictoc::toc()$toc - time_casebase$tic)
    
    tictoc::tic()
    fit_val_1se <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                               lambda = cv.lambda$lambda.1se , alpha = 0.7, unpen_cov = 2)
    time_casebase_1se <- as.numeric(tictoc::toc()$toc - time_casebase$tic)
    
    cli::cli_alert_success("Fit Case-base model org finished")
    
    
    coef_val1 <- fit_val_min$coefficients[1:eval(parse(text="p")), 1]
    coef_val2 <- fit_val_min$coefficients[1:eval(parse(text="p")), 2]
    
    coef_val3 <- fit_val_1se$coefficients[1:eval(parse(text="p")), 1]
    coef_val4 <- fit_val_1se$coefficients[1:eval(parse(text="p")), 2]
    
    # True coefs
    cb_true_min <- coef_val1[1:24]
    cb_true_1se <- coef_val3[1:24]
    
    #coef_val1 <- ifelse(coef_val1 < 1e-05, 0, coef_val1)
    #coef_val2 <- ifelse(coef_val2 < 1e-05, 0, coef_val2)
    
    #coef_val3 <- ifelse(coef_val1 < 1e-05, 0, coef_val1)
    #coef_val4 <- ifelse(coef_val2 < 1e-05, 0, coef_val2)
    
    
    res_cb_min1 <- varsel_perc(coef_val1, beta1) %>% bind_cols(tibble(Time = time_casebase_cv + time_casebase_min))
    
    res_cb_min2 <- varsel_perc(coef_val2, beta2) %>% bind_cols(tibble(Time = time_casebase_cv + time_casebase_min))
    
    res_cb_min3 <- varsel_perc(coef_val3, beta1) %>% bind_cols(tibble(Time = time_casebase_cv + time_casebase_1se))
    
    res_cb_min4 <- varsel_perc(coef_val4, beta2) %>% bind_cols(tibble(Time = time_casebase_cv + time_casebase_1se))
    
    
    ########################## Fit casebase > tol model #################################
    tictoc::tic()
    # Train case-base model 
    cv.lambda2 <- mtool.multinom.cv(train, seed = 1, alpha = 0.7, nfold = 5, 
                                    train_ratio = 20,
                                    fit_fun = purrr::partial(mtool::mtool.MNlogistic, 
                                                             tol = 1e-6,
                                                             learning_rate = 1e-3))
    
    time_casebase2_cv <- as.numeric(tictoc::toc()$toc - time_casebase$tic)
    
    # Case-base fits 
    # Lambda.min
    tictoc::tic()
    fit2_val_min <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                               lambda = cv.lambda2 , alpha = 0.7, unpen_cov = 2,
                               fit_fun = purrr::partial(mtool::mtool.MNlogistic, 
                                                        tol = 1e-6,
                                                        learning_rate = 1e-3))
    
    
    time_casebase2_min <- as.numeric(tictoc::toc()$toc - time_casebase$tic)
    
    tictoc::tic()
    fit2_val_1se <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                               lambda = cv.lambda2$lambda.1se , alpha = 0.7, 
                               unpen_cov = 2,
                               fit_fun = purrr::partial(mtool::mtool.MNlogistic, 
                                                        tol = 1e-6,
                                                        learning_rate = 1e-3))
    time_casebase2_1se <- as.numeric(tictoc::toc()$toc - time_casebase$tic)
    
    cli::cli_alert_success("Fit Case-base model Acc. finished")
    
    
    coef2_val1 <- fit2_val_min$coefficients[1:eval(parse(text="p")), 1]
    coef2_val2 <- fit2_val_min$coefficients[1:eval(parse(text="p")), 2]
    
    coef2_val3 <- fit2_val_1se$coefficients[1:eval(parse(text="p")), 1]
    coef2_val4 <- fit2_val_1se$coefficients[1:eval(parse(text="p")), 2]
    
    # True coefs
    cb2_true_min <- coef2_val1[1:24]
    cb2_true_1se <- coef2_val3[1:24]
    
    #coef_val1 <- ifelse(coef_val1 < 1e-05, 0, coef_val1)
    #coef_val2 <- ifelse(coef_val2 < 1e-05, 0, coef_val2)
    
    #coef_val3 <- ifelse(coef_val1 < 1e-05, 0, coef_val1)
    #coef_val4 <- ifelse(coef_val2 < 1e-05, 0, coef_val2)
    
    
    res_cb2_min1 <- varsel_perc(coef2_val1, beta1) %>% bind_cols(tibble(Time = time_casebase2_cv + time_casebase2_min))
    
    res_cb2_min2 <- varsel_perc(coef2_val2, beta2) %>% bind_cols(tibble(Time = time_casebase2_cv + time_casebase2_min))
    
    res_cb2_min3 <- varsel_perc(coef2_val3, beta1) %>% bind_cols(tibble(Time = time_casebase2_cv + time_casebase2_1se))
    
    res_cb2_min4 <- varsel_perc(coef2_val4, beta2) %>% bind_cols(tibble(Time = time_casebase2_cv + time_casebase2_1se))
    
    
    # MSE
    mse_cb2 <- mean((fit2_val_min$coefficients[1:eval(parse(text="p")), 1][1:24] - beta1[1:24])^2)
    
    ###################################################################################
    Res <- rbind(res_cb_min1, res_cb_min2, res_cb_min3, res_cb_min4, 
                 res_cb2_min1, res_cb2_min2, res_cb2_min3, res_cb2_min4, 
                 res_cox_min1, res_cox_min2,   cen.prop)
    
    rownames(Res) <- c("casebase.lambda.min_cause1", "casebase.lambda.min_cause2","casebase.lambda.1se_cause1", "casebase.lambda.1se_cause2",
                       "cox.lambda.min_cause1", "cox.lambda.min_cause2", "cens.prop")
    
    Res
    
    mse_res <- rbind(mse_cox, mse_cb, mse_cb2)
    
    cli::cli_alert_success("Simulation w/ seed {seed} finished.")
    
    return(list(Res = Res, mse_res = mse_res, cc_true = cc_true, 
                cb_true_min = cb_true_min, 
                cb_true_1se = cb_true_1se,
                cb2_true_min = cb2_true_min, 
                cb2_true_1se = cb2_true_1se))
    
}




Res_all <- purrr::map2_dfr(results_list, seeds, ~{
    as.data.frame(.x$Res) %>%
        mutate(model = rownames(.x$Res), seed = .y)
})

mse_all <- purrr::map2_dfr(results_list, seeds, ~{
    as.data.frame(.x$mse_res) %>%
        rownames_to_column("model") %>%
        mutate(seed = .y)
})

Res_all %>%
    kable(format = "html", caption = "📊 Sensitivity and Specificity Across Seeds") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
    scroll_box(width = "100%", height = "400px")

mse_all %>%
    kable(format = "html", caption = "📈 MSE Across Models and Seeds") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive")) %>%
    scroll_box(width = "100%", height = "300px")

coef_truth_long <- purrr::map2_dfr(results_list, seeds, ~{
    tibble(
        seed = .y,
        index = 1:24,
        glmnet = as.numeric(.x$cc_true),
        mtool_min = as.numeric(.x$cb_true_min),
        mtool_1se = as.numeric(.x$cb_true_1se)
    )
}) %>%
    pivot_longer(cols = c(glmnet, mtool_min, mtool_1se),
                 names_to = "method", values_to = "coefficient")

coef_truth_long %>%
    mutate(across(everything(), as.character)) %>%  
    kable(format = "html", caption = "🔍 Comparison of True Coefficients Across Methods and Seeds") %>%
    kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), full_width = F) %>%
    scroll_box(width = "100%", height = "400px")

print("==== Sensitivity & Specificity Across Runs ====")

print(Res_all)

print("==== MSE for Each Model Across Runs ====")

print(mse_all)

print("==== True Coefficients from Each Method ====")
print(coef_truth_long)
