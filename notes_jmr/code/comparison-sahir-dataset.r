set.seed(12345)
library(casebase)
library(glmnet)
library(survival)

data(bmtcrr)
head(bmtcrr)


y <- Surv(time = bmtcrr$ftime, 
          event = as.numeric(bmtcrr$Status == 1))

x <- model.matrix(Status ~ . -ftime,
                  data = bmtcrr)

cox_enet_mod <- cv.glmnet(x = x, y = y, family = "cox",
                          # family = "binomial",
                          nfolds = 10,
                          alpha = 0.5)

coef(cox_enet_mod, s = cox_enet_mod$lambda.min)

plot(cox_enet_mod)

cc_enet_min <- coef(cox_enet_mod, s = cox_enet_mod$lambda.min)

select_vars_enet <- cc_enet_min@Dimnames[[1]][-1][cc_enet_min@i]

selected_coefs_enet <- cc_enet_min@x

names(selected_coefs_enet) <- select_vars_enet

selected_coefs_enet



bmtcrr_mtx <- model.matrix( ~ . -1,
                  data = bmtcrr)[,-1]

mtool_fit_fun <- purrr::partial(cbSCRIP::MNlogistic,
                                niter_inner_mtplyr = 1,
                                maxit = 200,
                                tolerance = 1e-5,
                                learning_rate = 1e-1,
                                verbose = F,
                                save_history = F)   

cv_multinom_enet <- cv_cbSCRIP(
    Surv(ftime, Status) ~ .,
    data = bmtcrr_mtx,
    alpha = 0.5,
    nfold = 5,
    nlambda = 50,
    warm_start = T,
    fit_fun = mtool_fit_fun,
    ratio = 50)

plot(cv_multinom_enet)

cv_multinom_enet$fit.min$coefficients

multinom_enet <- cbSCRIP(
    Surv(ftime, Status) ~ .,
    data = bmtcrr_mtx,
    alpha = 0.5,
    # lambda = .01,
    warm_start = F,
    fit_fun = mtool_fit_fun,
    ratio = 50)

plot(multinom_enet)

multinom_enet <- cbSCRIP(
    Surv(ftime, Status) ~ .,
    data = bmtcrr_mtx,
    alpha = 0.5,
    lambda = 0,
    warm_start = F,
    fit_fun = mtool_fit_fun,
    ratio = 50)
multinom_enet$coefficients

model1 <- fitSmoothHazard(Status ~ ftime + Sex + D + Phase + Source + Age, 
                          data = bmtcrr, 
                          ratio = 100,
                          time = "ftime")
summary(model1)

plot(cv_multinom_enet)

cv_multinom_enet$fit.min$coefficients

