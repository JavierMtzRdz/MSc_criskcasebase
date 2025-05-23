########################### n = 400, p = 1000, Tp = 20 Effect only on cause 1 ######################
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

# Fitting functions 
source("../fitting_functions.R")

# Set seed
seed <- as.integer(Sys.time())

# take the last five digits of the initial seed
the_seed= seed %% 100000
set.seed(the_seed)


# Setup
n <- 400
p <- 1000
num_true <- 20
beta1 <- c(rep(0, p))
beta2 <- c(rep(0, p))
nu_ind <- seq(num_true)
beta1[nu_ind] <- c(rep(1, 20))
beta2[nu_ind] <- c(rep(0, 20))

# Simulate data
sim.data <- cause_hazards_sim(n = n, p = p, nblocks = 4, 
                              beta1 = beta1, beta2 = beta2, rate_cens = 0.25, 
                              h1 = 0.55, h2 = 0.35, gamma1 = 1.5, gamma2 = 1.5, exchangeable = TRUE)


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

# Fit cause-specific cox model with glmnet on training set 
cox_mod1 <- cv.glmnet(x = x_train1, y = y_train1, family = "cox", alpha = 0.7)

# Fit on validation set 
cox_val_min1 <- glmnet(x = x_test1, y = y_test1, family = "cox", alpha = 0.7, 
                       lambda = cox_mod1$lambda.min)

cc_min1 <- coef(cox_val_min1)

res_cox_min1 <- varsel_perc(cc_min1, beta1)
########################## Cause 2 #####################################
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

cc_min2 <- coef(cox_val_min2)

res_cox_min2 <- varsel_perc(cc_min2, beta2)

########################## Fit PenCR model ##################################
penCR = cv.glmnet.CR(train, family="cox", alpha= 0.7, standardize= TRUE,
                     nlambda = 30, t.BS = median(train$ftime), seed = 115, causeOfInt = 1,
                     nfold = 10)

lambda_penCR1 <- penCR$glmnet.fits$models$`Cause 1`$glmnet.res$lambda[penCR$min.index[1]]
lambda_penCR2 <- penCR$glmnet.fits$models$`Cause 2`$glmnet.res$lambda[penCR$min.index[2]]

# Fit on validation set 
res_pencr_min1 <- glmnet(x = x_test1, y = y_test1, family = "cox", alpha = 0.7, 
                         lambda = lambda_penCR1)

res_pencr_min2 <- glmnet(x = x_test2, y = y_test2, family = "cox", alpha = 0.7, 
                         lambda = lambda_penCR2)

cc_pencr_min1 <- coef(res_pencr_min1)
cc_pencr_min2 <- coef(res_pencr_min2)


res_pencr_min1 <- varsel_perc(cc_pencr_min1, beta1)
res_pencr_min2 <- varsel_perc(cc_pencr_min2, beta2)

########################## Fit casebase model #################################
# Train case-base model 
cv.lambda <- mtool.multinom.cv(train, seed = 1, alpha = 0.7, nfold = 5)


cv.lambda

# Test set 
surv_obj_val <- with(test, Surv(ftime, as.numeric(fstatus), type = "mstate"))

# Covariance matrix
cov_val <- cbind(test[, c(grepl("X", colnames(test)))], time = log(test$ftime))

# Case-base dataset
cb_data_val <- create_cbDataset(surv_obj_val, as.matrix(cov_val), ratio = 10)

# Case-base fits 
# Lambda.min
fit_val_min <- fit_cbmodel(cb_data_val, regularization = 'elastic-net',
                           lambda = cv.lambda$lambda.min , alpha = 0.7, unpen_cov = 2)


res_cb_min1 <- varsel_perc(fit_val_min$coefficients[1:eval(parse(text="p")), 1], beta1)

res_cb_min2 <- varsel_perc(fit_val_min$coefficients[1:eval(parse(text="p")), 2], beta2)


###################################################################################
Res <- rbind(res_cb_min1, res_cb_min2, res_cox_min1, res_cox_min2,  res_pencr_min1, res_pencr_min2, cen.prop)

rownames(Res) <- c("casebase.lambda.min_cause1", "casebase.lambda.min_cause2", "cox.lambda.min_cause1",
                   "cox.lambda.min_cause2", "pencr.lambda.mincause1", "pencr.lambda.mincause2", "cens.prop")

Res

write.csv(Res, file = glue("setting10{runif(1)}.csv"))
