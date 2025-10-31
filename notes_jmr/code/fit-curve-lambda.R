##################################################################
##             Project: Selection Curves Simulation             ##
##################################################################
##
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
# Load packages
pacman::p_load(
    casebase,
    glmnet,
    parallel,
    tictoc,
    dplyr,
    tidyr,
    foreach,
    survival,
    cmprsk,
    here,
    caret,
    pec,
    purrr,
    glue,
    cbSCRIP
)

# Global parameters
n_lambda <- 50


## Sim parameters


# Simulation-specific parameters
n <- 400 # Number of samples

# Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
p <- as.numeric(args[1])
num_true <- as.integer(args[2])
setting <- as.numeric(args[3])
i <- as.numeric(args[4])
rm(args)

# Generate data
data <- gen_data(
    n = n,
    p = p,
    num_true = num_true,
    setting = setting,
    iter = i
)

# Unpack data
beta1 <- data$beta1
beta2 <- data$beta2
train <- data$train
test <- data$test


## 3.  Cause-Specific Cox model ----

# Prepare data for cause-1 (censor competing event)
y_train <- Surv(time = train$ftime, event = train$fstatus == 1)
x_train <- model.matrix(~ . - ftime - fstatus, data = train)[, -1]

y_test <- Surv(time = test$ftime, event = test$fstatus == 1)
x_test <- model.matrix(~ . - ftime - fstatus, data = test)[, -1]

# Fit cause-specific cox model
tic("cox_mod")
cox_mod <- glmnet(
    x = x_train,
    y = y_train,
    family = "cox",
    alpha = 0.5,
    lambda.min.ratio = ifelse(length(y_test) < ncol(x_test), 0, 0.001),
    nlambda = 50
)
toc()

# Process Cox model results
table_result_cox <- map_df(cox_mod$lambda, \(x) {
    as_tibble(varsel_perc(coef(cox_mod, s = x), beta1)) %>%
        mutate(
            mse_bias = mse_bias(coef(cox_mod, s = x), beta1),
            lambda = x,
            params = list(coef(cox_mod, s = x))
        )
}) %>%
    mutate(
        model = "enet-CR",
        sim = i,
        model_size = TP + FP
    )

# Quick plot (ROC)
table_result_cox %>%
    ggplot(aes(x = 1 - Specificity, y = Sensitivity)) +
    geom_path() +
    coord_equal()


## 4. cbSCRIP - enet) ----

# Fit the cbSCRIP model
tic("cv.lambdaAcc")
set.seed(123)
cv.lambdaAcc <- cbSCRIP(
    Surv(ftime, fstatus) ~ .,
    train,
    coeffs = "original"
)
toc()

# Process cbSCRIP results
table_result_cb <- map_df(seq_along(cv.lambdaAcc$lambdagrid), \(x) {
    as_tibble(varsel_perc(cv.lambdaAcc$coefficients[[x]][1:length(beta1), 1],
                          beta1)) %>%
        mutate(
            mse_bias = mse_bias(cv.lambdaAcc$coefficients[[x]][1:length(beta1), 1], beta1),
            lambda = cv.lambdaAcc$lambdagrid[x],
            params = list(cv.lambdaAcc$coefficients[[x]][1:length(beta1), 1])
        )
}) %>%
    mutate(
        model = "enet-casebase-Acc",
        sim = i,
        model_size = TP + FP
    )

# Plot coefficient paths
plot(cv.lambdaAcc, plot_intercept = F) +
    geom_vline(xintercept = cv.lambdaAcc$lambdagrid,
               alpha = 0.5)

# Quick plots (Sensitivity vs. Model Size / Specificity)
table_result_cb %>%
    ggplot(aes(x = model_size, y = Sensitivity)) +
    geom_path()

table_result_cb %>%
    ggplot(aes(x = 1 - Specificity, y = Sensitivity)) +
    geom_path() +
    coord_equal()


# --- Commented-out diagnostics & alternative processing ---
# X_t <- model.matrix(~.,
#   data = data.frame(cbind(cv.lambdaAcc$cb_data$covariates,
#                           time = log(cv.lambdaAcc$cb_data$time)))
# )
#
# X_t <- cbind(X_t[, -1], X_t[, 1])
#
# p.fac <- rep(1, ncol(X_t))
# # p.fac[(1388+1-pen_last):1388] <- 0
#
# cv.lambdaAcc$models_info[[1]]$convergence_pass
# cv.lambdaAcc$models_info[[1]]$coefficients_sparse
# loss <- map_vec(cv.lambdaAcc$models_info[[1]]$coefficients_history,
#   ~ calculate_penalized_multinomial_loss(
#     .x,
#     alpha = 0.5,
#     lambda = 0.01,
#     Y = cv.lambdaAcc$cb_data$event,
#     X = X_t,
#     offset = cv.lambdaAcc$cb_data$offset,
#     penalty_weights = p.fac
#   )
# )
#
# if (min(loss) < min_loss) min_loss <- min(loss)
#
# tibble(iter = 1:length(loss),
#        loss = loss) %>%
#   mutate(loss_diff = loss - min_loss) %>%
#   ggplot(aes(x = iter, y = loss_diff)) +
#   geom_line() +
#   scale_y_log10() +
#   labs(y = "F(x)-F(x*)")
#
# cv.lambdaAcc$coefficients[[50]] |> View()
#
# table_result_cb <- map_df(cv.lambdaAcc$lambdagrid,
#   \(x) {
#     as_tibble(varsel_perc(cv.lambdaAcc$coefficients[[as.character(x)]][, 1],
#                           beta1)) %>%
#       mutate(
#         mse_bias = mse_bias(cv.lambdaAcc$coefficients[[as.character(x)]][, 1], beta1),
#         lambda = x,
#         params = list(cv.lambdaAcc$coefficients[[as.character(x)]][, 1])
#       )
#   }
# ) %>%
#   mutate(
#     model = "enet-casebase-Acc",
#     sim = i,
#     model_size = TP + FP
#   )
# ---

#================================================================
## 5. MODEL 3: CASE-BASE (cbSCRIP - SCAD) [COMMENTED] ----
#================================================================

# tic("cv.lambdaSCAD")
# cv.lambdaSCAD <- cbSCRIP(Surv(ftime, fstatus) ~ .,
#   train,
#   lambda.min.ratio = 0.01,
#   nlambda = 50,
#   warm_start = T,
#   regularization = "SCAD",
#   fit_fun = fit_fun
# )
# toc()
# cv.lambdaSCAD$lambdagrid
# plot(cv.lambdaSCAD) +
#   geom_vline(xintercept = cv.lambdaSCAD$lambdagrid)
# # plot_cv.multinom(cv.lambdaSCAD)
#
# table_result_SCAD <- map_df(cv.lambdaSCAD$lambdagrid,
#   \(x) {
#     as_tibble(varsel_perc(cv.lambdaSCAD$coefficients[[as.character(x)]][, 1],
#                           beta1)) %>%
#       mutate(
#         mse_bias = mse_bias(cv.lambdaSCAD$coefficients[[as.character(x)]][, 1], beta1),
#         lambda = x,
#         params = list(cv.lambdaSCAD$coefficients[[as.character(x)]][, 1])
#       )
#   }
# ) %>%
#   mutate(
#     model = "SCAD-casebase",
#     sim = i,
#     model_size = TP + FP
#   )
#
# table_result_SCAD %>%
#   ggplot(aes(x = 1 - Specificity, y = Sensitivity)) +
#   geom_path() +
#   coord_equal()
#
#
# surv_obj_val <- with(train, Surv(ftime, as.numeric(fstatus), type = "mstate"))
#
# # Covariance matrix
# cov_val <- cbind(train[, c(grepl("X", colnames(train)))], time = log(train$ftime))
#
# # Case-base dataset
# cb_data_val <- create_cbDataset(surv_obj_val, as.matrix(cov_val),
#   ratio = 20
# )

#================================================================
## 6. COMBINE, PLOT, AND SAVE RESULTS ----
#================================================================

# --- Commented-out logic for loading/merging SCAD results ---
# table_select <- readRDS(here("paper",
#   "results",
#   "lambda_paths",
#   glue("models_selec-{setting}_iter-{i}_p-{p}_k-{num_true}.rds")
# )) %>%
#   filter(model != "SCAD-casebase")
#
# cli::cli_alert_success("SCAD path ran!")
#
# table_select <- as_tibble(bind_rows(
#   table_select,
#   table_result_SCAD
# ))
# ---

# Combine the results from the active models
table_select <- as_tibble(rbind(
    table_result_cox,
    table_result_cb
    # table_result_SCAD
))

# Create final ROC plot
(plot <- table_select %>%
        ggplot(aes(
            x = 1 - Specificity,
            y = Sensitivity,
            color = model,
            linetype = model
        )) +
        geom_path() +
        coord_equal())

# Save plot
ggsave(
    here(
        "paper",
        "results",
        glue("models_selec-{setting}_iter-{i}_p-{p}_k-{num_true}.png")
    ),
    plot = plot,
    bg = "transparent",
    width = 200,
    height = 120,
    units = "mm",
    dpi = 300
)

# Save results object
saveRDS(
    table_select,
    here(
        "paper",
        "results",
        "lambda_paths",
        glue("models_selec-{setting}_iter-{i}_p-{p}_k-{num_true}.rds")
    )
)