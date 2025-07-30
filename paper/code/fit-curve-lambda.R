##################################################################
##             Project: Selection Curves Simulation             ##
##################################################################
##
## Description:    Texto
##                 Texto
##
## Author:         Javier Mtz.-Rdz.  
##
## Creation date:  2025-07-10
##
## Email:          javier.mr@stat.ubc.ca
##
## ---------------------------
## Notes:          
## ---------------------------

# Setup ----
## Packages to use ----

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
               pec,
               purrr,
               glue)

# Fitting functions 
source(here("paper",
            "code", "fitting_functionsV2.R"))

if(!exists("start_sims")) start_sims <- 1

n_lambda <- 50


for (i in start_sims:sims) {
    
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
    
    # Fit indepedent cox-regression model
    # Cause-1
    # Censor competing event
    y_train <- Surv(time = train$ftime, event = train$fstatus == 1)
    
    x_train <- model.matrix(~ . -ftime -fstatus, data = train)[, -1]
    
    # Censor competing event
    y_test <- Surv(time = test$ftime, event = test$fstatus == 1)
    
    x_test <- model.matrix(~ . -ftime -fstatus, data = test)[, -1]
    
    # Fit cause-specific cox model with glmnet on training set
    cox_mod <- glmnet(x = x_train, y = y_train, family = "cox", alpha = 0.7,
                      lambda.min.ratio = 0,
                      nlambda = n_lambda)


    table_result_cox <- map_df(cox_mod$lambda,
                               \(x){as_tibble(varsel_perc(coef(cox_mod, s = x),
                                                      beta1),
                                          lambda = x) %>%
                                       mutate(mse_bias = mse_bias(coef(cox_mod, s = x), beta1))}) %>%
        mutate(model = "enet-CR",
               sim = i)

    # casebase

    mtool_cv <- purrr::partial(mtool::mtool.MNlogisticAcc,
                               niter_inner_mtplyr = 1.5,
                               maxit = 100,
                               momentum_gamma = 0,
                               tolerance = 1e-8,
                               learning_rate = 1e-4
    )

    tic("cv.lambdaAcc")
    cv.lambdaAcc <- mtool.multinom(train, seed = 1, alpha = 0.7,
                                   fit_fun = mtool_cv,
                                   ncores = 2,
                                   # train_ratio = 5,
                                   # lambda_max = 0.5,
                                   # ws = F,
                                   grid_size = n_lambda
    )
    toc()

    # plot_cv.multinom(cv.lambdaAcc)


    table_result_cb <-map_df(cv.lambdaAcc$lambdagrid,
                                \(x){as_tibble(varsel_perc(cv.lambdaAcc$cov_coeffs[[as.character(x)]][,1],
                                                           beta1),
                                               lambda = x) %>%
                                        mutate(mse_bias = mse_bias(cv.lambdaAcc$cov_coeffs[[as.character(x)]][,1], beta1))}) %>%
        mutate(model = "enet-casebase-Acc",
               sim = i)

    table_result_cb %>%
        ggplot(aes(x = 1 - Specificity, y = Sensitivity)) +
        geom_path() +
        coord_equal()
    
    
    # SCAD -----
    mtool_cv <- purrr::partial(mtool::mtool.MNlogisticAcc,
                               niter_inner_mtplyr = 1.5,
                               maxit = 100,
                               momentum_gamma = 0,
                               tolerance = 1e-8,
                               learning_rate = 1e-4
    )
    tic("cv.lambdaSCAD")
    cv.lambdaSCAD <- mtool.multinom(train, seed = 1, 
                                    fit_fun = mtool_cv,
                                    grid_size = n_lambda,
                                    regularization = "SCAD",
                                    ncores = 2,
                                    train_ratio = 5
    )
    toc()
    
    table_result_SCAD <- map_df(cv.lambdaSCAD$lambdagrid,
                                \(x){as_tibble(varsel_perc(cv.lambdaSCAD$cov_coeffs[[as.character(x)]][,1],
                                                       beta1),
                                           lambda = x) %>% 
                                        mutate(mse_bias = mse_bias(cv.lambdaSCAD$cov_coeffs[[as.character(x)]][,1], beta1))}) %>% 
        mutate(model = "SCAD-casebase",
               sim = i)
    
    table_result_SCAD %>%
        ggplot(aes(x = 1 - Specificity, y = Sensitivity)) +
        geom_path() +
        coord_equal()
    
    
    
    # 
    # surv_obj_val <- with(train, Surv(ftime, as.numeric(fstatus), type = "mstate"))
    # 
    # # Covariance matrix
    # cov_val <- cbind(train[, c(grepl("X", colnames(train)))], time = log(train$ftime))
    # 
    # # Case-base dataset
    # cb_data_val <- create_cbDataset(surv_obj_val, as.matrix(cov_val))
    # 
    # # Case-base fits 
    # # Lambda.min
    # lambda <- 0.2
    # cli::cli_alert_info("Lambda.min: {lambda}")
    # fit_val_min_SCAD <- fit_cbmodel(cb_data_val, regularization = 'SCAD',
    #                                lambda = lambda, #alpha = 0.7,
    #                                unpen_cov = 2,
    #                                fit_fun = mtool_cv)
    # 
    # varsel_perc(fit_val_min_SCAD$coefficients[1:eval(parse(text="p")), 1], beta1)
    
    
    
    
    
    
    # table_select <- readRDS(here("paper",
    #              "results",
    #              "lambda_paths",
    #              glue("models_selec-{setting}_iter-{i}_p-{p}_k-{num_true}.rds"))) %>%
    #     filter(model != "SCAD-casebase")
    # 
    # cli::cli_alert_success("SCAD path ran!")
    # 
    # table_select <- as_tibble(bind_rows(table_select,
    #                                 table_result_SCAD))
    
    table_select <- as_tibble(rbind(table_result_cox,
                                    table_result_cb,
                                    table_result_SCAD))
    
    print(table_select %>%
        ggplot(aes(x = 1 - Specificity, y = Sensitivity,
                   color = model)) +
        geom_path() +
        coord_equal())
    
    
    saveRDS(table_select,
            here("paper",
                 "results",
                 "lambda_paths",
                 glue("models_selec-{setting}_iter-{i}_p-{p}_k-{num_true}.rds")))
    
}
