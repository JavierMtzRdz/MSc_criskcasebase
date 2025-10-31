
library(casebase)
library(future.apply)
library(glmnet)
#library(mtool)
library(cbSCRIP)
library(parallel)
library(tictoc)
library(tidyverse)
library(foreach)
library(survival)
library(cmprsk)
library(glue)
library(pec)
library(CoxBoost)
library(riskRegression)

data <- cbSCRIP:::gen_data(n = 400, p = 120, num_true = 20, setting = 1, iter = 123)

beta1 <- data$beta1
beta2 <- data$beta2
train <- data$train
test <- data$test
p <- ncol(train) -2

tic()
fit_cv <- cv_cbSCRIP(Surv(ftime, fstatus) ~ .,
                          train,
                          nlambda = 50,
                          alpha = 0.5,
                          ratio = 50,
                     fit_fun = MNlogisticSAGA_Native,
                     lr_adj = 1,
                     warm_start = T)
toc()

(res_cb_min2SE_1 <-varsel_perc(fit_cv$fit.min$coefficients[1:p, 1], 
                              beta1))

(bias_min2SE_1 <-mse_bias(fit_cv$fit.min$coefficients[1:p, 1], 
                               beta1))

fit_cv$fit.min$convergence_pass

plot(fit_cv)


fit_min <- cbSCRIP(Surv(ftime, fstatus) ~ .,
                    train,
                    lambda = 0,
                    alpha = 0.5,
                    ratio = 50,
                    coeffs = "original")
fit_min$convergence_pass
fit_min$coefficients

tic()
fit_path <- cbSCRIP(Surv(ftime, fstatus) ~ .,
                          train,
                    coeffs = "original",
                    lr_adj = 1,
                    tolerance = 1e-4,
                    maxit = 1000,
                    fit_fun = cbSCRIP::MNlogisticSAGA_Native)
toc()

plot(fit_path)

map_dbl(fit_path$models_info, ~.$convergence_pass)

map_dbl(fit_path$non_zero, ~.)

map_dfr(1:length(fit_path$coefficients), 
        ~varsel_perc(fit_path$models_info[[.]]$coefficients[1:p, 1], 
                     beta1) |> 
            mutate(lambda = fit_path$lambdagrid[.])) 

map_dfr(1:length(fit_path$coefficients), 
        ~varsel_perc(fit_path$models_info[[.]]$coefficients[1:p, 1], 
                     beta1) |> 
            mutate(lambda = fit_path$lambdagrid[.])) |> 
    ggplot(aes(x = 1-Specificity, y = Sensitivity)) +
    geom_line()

map_dfr(1:length(fit_path$coefficients), 
        ~varsel_perc(fit_path$models_info[[.]]$coefficients[1:p, 1], 
                     beta1) |> 
            mutate(lambda = fit_path$lambdagrid[.])) |> 
    ggplot(aes(x = TP + FP, y = Sensitivity)) +
    geom_line()


