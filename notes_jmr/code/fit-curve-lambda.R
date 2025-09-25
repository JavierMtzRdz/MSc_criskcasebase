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
               glue,
               cbSCRIP)

# Fitting functions 
source(here("paper",
            "code", "fitting_functionsV2.R"))
source(here("paper",
            "code", "fitting_functionsV3.R"))

n_lambda <- 50

# for (i in start_sims:sims) {
    
n <- 400          # Number of samples

args <- commandArgs(trailingOnly = TRUE)

p <- as.numeric(args[1])

num_true <- as.integer(args[2])

setting <- as.numeric(args[3])

i <- as.numeric(args[4])

rm(args)

data <- gen_data(n = n, 
                 p = p, 
                 num_true = num_true,
                 setting = setting,
                 iter = i)

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
    tic("cox_mod")
    cox_mod <- glmnet(x = x_train, y = y_train, family = "cox", alpha = 0.5,
                      lambda.min.ratio = ifelse(length(y_test) < ncol(x_test), 
                                                0, 0.001),
                      nlambda = 50)
    toc()
    
    
    cox_mod0 <- glmnet(x = x_train, y = y_train, family = "cox", alpha = 0.5,
                      lambda = 0)
    

    
    table_result_cox <- map_df(cox_mod$lambda,
                               \(x){as_tibble(varsel_perc(coef(cox_mod, s = x),
                                                      beta1)) %>%
                                       mutate(mse_bias = mse_bias(coef(cox_mod, s = x), beta1),
                                              lambda = x,
                                              params = list(coef(cox_mod, s = x)))}) %>% 
        bind_rows(as_tibble(varsel_perc(coef(cox_mod0, s = 0),
                                        beta1)) %>%
                      mutate(mse_bias = mse_bias(coef(cox_mod0, s = 0), beta1),
                             lambda = 0,
                             params = list(coef(cox_mod0, s = 0)))) %>% 
        mutate(model = "enet-CR",
               sim = i,
               model_size = TP + FP)
    
    
    table_result_cox %>%
        ggplot(aes(x = 1 - Specificity, y = Sensitivity)) +
        geom_path() +
        coord_equal()

    
    # casebase ----
    if(p <500){
        fit_fun <- purrr::partial(cbSCRIP::MNlogisticAcc,
                                  niter_inner_mtplyr = 1.5,
                                  maxit = 250,
                                  c_factor =500000,
                                  v_factor = 100000,
                                  tolerance = 1e-3,
                                  save_history = F,
                                  verbose = F)
    } else if(p <700){
    fit_fun <- purrr::partial(cbSCRIP::MNlogisticAcc,
                               niter_inner_mtplyr = 1,
                               maxit = 250,
                               c_factor = 100000,
                               v_factor = 500,
                               tolerance = 1e-3,
                               save_history = F,
                               verbose = F)
    }else if(p > 700){
        fit_fun <- purrr::partial(cbSCRIP::MNlogisticAcc,
                                   niter_inner_mtplyr = 0.5,
                                   maxit = 250,
                                   c_factor = 2000,
                                   v_factor = 5000,
                                   tolerance = 1e-3,
                                   save_history = F,
                                   verbose = F)
    }

    
    # fit_fun <- purrr::partial(cbSCRIP::MNlogisticAcc,
    #                           niter_inner_mtplyr = 1.5,
    #                           maxit = 500,
    #                           c_factor =1000000,
    #                           v_factor = 500000,
    #                           tolerance = 1e-3,
    #                           save_history = F,
    #                           verbose = F)
    # 
    # fit_fun <- purrr::partial(cbSCRIP::MNlogistic,
    #                           niter_inner_mtplyr = 2,
    #                           maxit = 200,
    #                           tolerance = 1e-4,
    #                           learning_rate = 1e-4,
    #                           verbose = F,
    #                           save_history = F)   
    
    tic("cv.lambdaAcc")
    set.seed(123)
    cv.lambdaAcc <- cbSCRIP(Surv(ftime, fstatus) ~ .,
                            train,
                            nlambda = 50,
                            alpha = 0.5,
                            warm_start = F,
                            fit_fun = fit_fun,
                            ratio = 50)
    toc()
    
    table_result_cb <- map_df(seq_along(cv.lambdaAcc$lambdagrid),
                              \(x){
                                  as_tibble(varsel_perc(cv.lambdaAcc$coefficients[[x]][,1],
                                                        beta1)) %>%
                                      mutate(mse_bias = mse_bias(cv.lambdaAcc$coefficients[[x]][,1], beta1),
                                             lambda = cv.lambdaAcc$lambdagrid[x],
                                             params = list(cv.lambdaAcc$coefficients[[x]][,1]))}) %>%
        mutate(model = "enet-casebase-Acc",
               sim = i,
               model_size = TP + FP)
    
    table_select <- as_tibble(rbind(table_result_cox,
                                    table_result_cb
                                    # table_result_SCAD
    ))
    
    (plot <- table_select %>%
            ggplot(aes(x = 1 - Specificity, y = Sensitivity,
                       color = model,
                       linetype = model)) +
            geom_path() +
            coord_equal())
    
 # X_t <- model.matrix(~., 
    #                     data = data.frame(cbind(cv.lambdaAcc$cb_data$covariates, 
    #                                             time = log(cv.lambdaAcc$cb_data$time))))
    # 
    # X_t <- cbind(X_t[,-1], X_t[,1])
    # 
    # p.fac <- rep(1, ncol(X_t))
    # # p.fac[(1388+1-pen_last):1388] <- 0
    # 
    # cv.lambdaAcc$models_info[[1]]$convergence_pass
    # cv.lambdaAcc$models_info[[1]]$coefficients_sparse
    # loss <- map_vec(cv.lambdaAcc$models_info[[1]]$coefficients_history,
    #                 ~calculate_penalized_multinomial_loss(
    #                     .x,
    #                     alpha = 0.5,
    #                     lambda = 0.01,
    #                     Y = cv.lambdaAcc$cb_data$event,
    #                     X = X_t,
    #                     offset = cv.lambdaAcc$cb_data$offset,
    #                     penalty_weights = p.fac))
    # 
    # if(min(loss) < min_loss) min_loss <- min(loss)
    # 
    # tibble(iter = 1:length(loss),
    #        loss = loss) %>% 
    #     mutate(loss_diff = loss - min_loss) %>% 
    #     ggplot(aes(x = iter, y = loss_diff)) +
    #     geom_line() +
    #     scale_y_log10() +
    #     labs(y = "F(x)-F(x*)")
    
    
    
    plot(cv.lambdaAcc) +
        geom_vline(xintercept = cv.lambdaAcc$lambdagrid,
                   alpha = 0.5)
    
    # table_result_cb <- map_df(cv.lambdaAcc$lambdagrid,
    #                             \(x){
    #                                 as_tibble(varsel_perc(cv.lambdaAcc$coefficients[[as.character(x)]][,1],
    #                                                        beta1)) %>%
    #                                     mutate(mse_bias = mse_bias(cv.lambdaAcc$coefficients[[as.character(x)]][,1], beta1),
    #                                            lambda = x,
    #                                            params = list(cv.lambdaAcc$coefficients[[as.character(x)]][,1]))}) %>%
    #     mutate(model = "enet-casebase-Acc",
    #            sim = i,
    #            model_size = TP + FP)

    table_result_cb %>%
        ggplot(aes(x = model_size, y = Sensitivity)) +
        geom_path() 
    
    table_result_cb %>%
        ggplot(aes(x = 1- Specificity, y = Sensitivity)) +
        geom_path() +
        coord_equal()
    
    
    # SCAD -----
  
    # tic("cv.lambdaSCAD")
    # cv.lambdaSCAD <- cbSCRIP(Surv(ftime, fstatus) ~ .,
    #                         train,
    #                         lambda.min.ratio = 0.01,
    #                         nlambda = 50,
    #                         warm_start = T,
    #                         regularization = "SCAD",
    #                         fit_fun = fit_fun)
    # toc()
    # cv.lambdaSCAD$lambdagrid
    # plot(cv.lambdaSCAD) +
    #     geom_vline(xintercept = cv.lambdaSCAD$lambdagrid)
    # # plot_cv.multinom(cv.lambdaSCAD)
    # 
    # table_result_SCAD <- map_df(cv.lambdaSCAD$lambdagrid,
    #                             \(x){as_tibble(varsel_perc(cv.lambdaSCAD$coefficients[[as.character(x)]][,1],
    #                                                    beta1)) %>%
    #                                     mutate(mse_bias = mse_bias(cv.lambdaSCAD$coefficients[[as.character(x)]][,1], beta1),
    #                                            lambda = x,
    #                                            params = list(cv.lambdaSCAD$coefficients[[as.character(x)]][,1]))}) %>%
    #     mutate(model = "SCAD-casebase",
    #            sim = i,
    #            model_size = TP + FP)
    # 
    # table_result_SCAD %>%
    #     ggplot(aes(x = 1 - Specificity, y = Sensitivity)) +
    #     geom_path() +
    #     coord_equal()
    # 
    # 
    # surv_obj_val <- with(train, Surv(ftime, as.numeric(fstatus), type = "mstate"))
    # 
    # # Covariance matrix
    # cov_val <- cbind(train[, c(grepl("X", colnames(train)))], time = log(train$ftime))
    # 
    # # Case-base dataset
    # cb_data_val <- create_cbDataset(surv_obj_val, as.matrix(cov_val),
    #                                 ratio = 20)

    
    
    
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
                                    table_result_cb
                                    # table_result_SCAD
                                    ))
    
    (plot <- table_select %>%
        ggplot(aes(x = 1 - Specificity, y = Sensitivity,
                   color = model,
                   linetype = model)) +
        geom_path() +
        coord_equal())
    
    ggsave(here("paper",
                "results",
                glue("models_selec-{setting}_iter-{i}_p-{p}_k-{num_true}.png")),
           plot = plot,
           bg = "transparent",
           width = 200,                 # Ancho de la gráfica
           height = 120,
           units = "mm",
           dpi = 300)
    
    
    saveRDS(table_select,
            here("paper",
                 "results",
                 "lambda_paths",
                 glue("models_selec-{setting}_iter-{i}_p-{p}_k-{num_true}.rds")))

