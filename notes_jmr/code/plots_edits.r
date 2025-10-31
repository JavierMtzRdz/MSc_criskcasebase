
#generate title

settings_tbl <- tibble::tribble(
    ~setting,    ~desc,
    "setting1",  "Single effects on end-point of interest",
    "setting2", "Single effects on both end-points",
    "setting3", "Opposing effects",
    "setting4", "Mixture of effects",
    "setting5", "Non-proportional hazards"
)

cpl_palette <- c("cbSCRIP" = "#277DA1", 
                 "Aalen-Johansen" = "#43AA8B",
                 "Boosted Fine-Gray" = "#F9C74F",
                 "iCR" = "#f94144")

## Set theme ------
# mytidyfunctions::set_mytheme(text = element_text(family = "Times New Roman"))

# Brier score plot ----

N <- 400; p <- 120; k <- 20

brier_table <- map_df(1:5,
                      ~ {run_id <- sprintf("setting%d_N%d_p%d_k%d", ., N, p, k)
                      res <- readRDS(here::here("paper","preds", 
                                                paste0("results_", run_id, ".rds")))
                      res$brier_table |> 
                          mutate(setting = .)})

(plot <- brier_table |> 
    tibble() |> 
    ggplot(aes(times, Brier, colour = model,
               linewidth = ifelse(model == "cbSCRIP", 0.5, 0.1))) +
    geom_line() +
    facet_wrap(. ~ title, nrow = 2, scale = "free_x") +
    scale_colour_manual(values = cpl_palette) +
    scale_linewidth_continuous(range = c(0.5, 0.8),
                               guide = "none") +
    labs(x = "Follow-up time (years)",
         color = "Models",
         y = "Brier Score for Cause 1 Predictions") +
    theme(legend.position = c(0.8, 0.2)))


out_name   <- file.path(res_dir, sprintf("brier_grid_N%d_p%d_k%d.png", N, p, k))

ggsave(out_name, plot,
       width = 20, height = 12, units = "cm", dpi = 300)


N <- 400; p <- 500; k <- 20

brier_table <- map_df(1:5,
                      ~ {run_id <- sprintf("setting%d_N%d_p%d_k%d", ., N, p, k)
                      res <- readRDS(here::here("paper","preds", 
                                                paste0("results_", run_id, ".rds")))
                      res$brier_table |> 
                          mutate(setting = .)})

(plot <- brier_table |> 
        tibble() |> 
        ggplot(aes(times, Brier, colour = model,
                   linewidth = ifelse(model == "cbSCRIP", 0.5, 0.1))) +
        geom_line() +
        facet_wrap(. ~ title, nrow = 2, scale = "free_x") +
        scale_colour_manual(values = cpl_palette) +
        scale_linewidth_continuous(range = c(0.5, 0.8),
                                   guide = "none") +
        labs(x = "Follow-up time (years)",
             color = "Models",
             y = "Brier Score for Cause 1 Predictions") +
        theme(legend.position = c(0.8, 0.2)))


out_name   <- file.path(res_dir, sprintf("brier_grid_N%d_p%d_k%d.png", N, p, k))

ggsave(out_name, plot,
       width = 20, height = 12, units = "cm", dpi = 300)


N <- 400; p <- 500; k <- 84

brier_table <- map_df(1:5,
                      ~ {run_id <- sprintf("setting%d_N%d_p%d_k%d", ., N, p, k)
                      res <- readRDS(here::here("paper","preds", 
                                                paste0("results_", run_id, ".rds")))
                      res$brier_table |> 
                          mutate(setting = .)})

(plot <- brier_table |> 
        tibble() |> 
        ggplot(aes(times, Brier, colour = model,
                   linewidth = ifelse(model == "cbSCRIP", 0.5, 0.1))) +
        geom_line() +
        facet_wrap(. ~ title, nrow = 2, scale = "free_x") +
        scale_colour_manual(values = cpl_palette) +
        scale_linewidth_continuous(range = c(0.5, 0.8),
                                   guide = "none") +
        labs(x = "Follow-up time (years)",
             color = "Models",
             y = "Brier Score for Cause 1 Predictions") +
        theme(legend.position = c(0.8, 0.2)))


out_name   <- file.path(res_dir, sprintf("brier_grid_N%d_p%d_k%d.png", N, p, k))

ggsave(out_name, plot,
       width = 20, height = 12, units = "cm", dpi = 300)



N <- 400; p <- 1000; k <- 20

brier_table <- map_df(1:5,
                      ~ {run_id <- sprintf("setting%d_N%d_p%d_k%d", ., N, p, k)
                      res <- readRDS(here::here("paper","preds", 
                                                paste0("results_", run_id, ".rds")))
                      res$brier_table |> 
                          mutate(setting = .)})

(plot <- brier_table |> 
        tibble() |> 
        ggplot(aes(times, Brier, colour = model,
                   linewidth = ifelse(model == "cbSCRIP", 0.5, 0.1))) +
        geom_line() +
        facet_wrap(. ~ title, nrow = 2, scale = "free_x") +
        scale_colour_manual(values = cpl_palette) +
        scale_linewidth_continuous(range = c(0.5, 0.8),
                                   guide = "none") +
        labs(x = "Follow-up time (years)",
             color = "Models",
             y = "Brier Score for Cause 1 Predictions") +
        theme(legend.position = c(0.8, 0.2)))


out_name   <- file.path(res_dir, sprintf("brier_grid_N%d_p%d_k%d.png", N, p, k))

ggsave(out_name, plot,
       width = 20, height = 12, units = "cm", dpi = 300)




N <- 400; p <- 1000; k <- 168

brier_table <- map_df(1:5,
                      ~ {run_id <- sprintf("setting%d_N%d_p%d_k%d", ., N, p, k)
                      res <- readRDS(here::here("paper","preds", 
                                                paste0("results_", run_id, ".rds")))
                      res$brier_table |> 
                          mutate(setting = .)})

(plot <- brier_table |> 
        tibble() |> 
        ggplot(aes(times, Brier, colour = model,
                   linewidth = ifelse(model == "cbSCRIP", 0.5, 0.1))) +
        geom_line() +
        facet_wrap(. ~ title, nrow = 2, scale = "free_x") +
        scale_colour_manual(values = cpl_palette) +
        scale_linewidth_continuous(range = c(0.5, 0.8),
                                   guide = "none") +
        labs(x = "Follow-up time (years)",
             color = "Models",
             y = "Brier Score for Cause 1 Predictions") +
        theme(legend.position = c(0.8, 0.2)))


out_name   <- file.path(res_dir, sprintf("brier_grid_N%d_p%d_k%d.png", N, p, k))

ggsave(out_name, plot,
       width = 20, height = 12, units = "cm", dpi = 300)



# CIF score plot ----

N <- 400; p <- 120; k <- 20

cif_table <- map_df(1:5,
                      ~ {cif_data_fn(N, p, k, .) |> 
                          mutate(setting = .)})

(plot <- cif_table |> 
        tibble() |> 
        mutate(title = settings_tbl$desc[setting]) |> 
        ggplot(aes(Time, Risk, colour = Method,
                   linewidth = ifelse(Method == "cbSCRIP", 0.5, 0.1))) +
        geom_line() +
        facet_wrap(. ~ title, nrow = 2, scale = "free_x") +
        scale_colour_manual(values = cpl_palette) +
        scale_linewidth_continuous(range = c(0.5, 0.8),
                                   guide = "none") +
        labs(x = "Follow-up time (years)",
             color = "Models",
             y = "Absolute Risk") +
        theme(legend.position = c(0.8, 0.2)))


out_name   <- file.path(res_dir, sprintf("cif_grid_N%d_p%d_k%d.png", N, p, k))

ggsave(out_name, plot,
       width = 20, height = 12, units = "cm", dpi = 300)

N <- 400; p <- 500; k <- 20

cif_table <- map_df(1:5,
                    ~ {cif_data_fn(N, p, k, .) |> 
                            mutate(setting = .)})

(plot <- cif_table |> 
        tibble() |> 
        mutate(title = settings_tbl$desc[setting]) |> 
        ggplot(aes(Time, Risk, colour = Method,
                   linewidth = ifelse(Method == "cbSCRIP", 0.5, 0.1))) +
        geom_line() +
        facet_wrap(. ~ title, nrow = 2, scale = "free_x") +
        scale_colour_manual(values = cpl_palette) +
        scale_linewidth_continuous(range = c(0.5, 0.8),
                                   guide = "none") +
        labs(x = "Follow-up time (years)",
             color = "Models",
             y = "Absolute Risk") +
        theme(legend.position = c(0.8, 0.2)))


out_name   <- file.path(res_dir, sprintf("cif_grid_N%d_p%d_k%d.png", N, p, k))

ggsave(out_name, plot,
       width = 20, height = 12, units = "cm", dpi = 300)


N <- 400; p <- 500; k <- 84

cif_table <- map_df(1:5,
                    ~ {cif_data_fn(N, p, k, .) |> 
                            mutate(setting = .)})

(plot <- cif_table |> 
        tibble() |> 
        mutate(title = settings_tbl$desc[setting]) |> 
        ggplot(aes(Time, Risk, colour = Method,
                   linewidth = ifelse(Method == "cbSCRIP", 0.5, 0.1))) +
        geom_line() +
        facet_wrap(. ~ title, nrow = 2, scale = "free_x") +
        scale_colour_manual(values = cpl_palette) +
        scale_linewidth_continuous(range = c(0.5, 0.8),
                                   guide = "none") +
        labs(x = "Follow-up time (years)",
             color = "Models",
             y = "Absolute Risk") +
        theme(legend.position = c(0.8, 0.2)))


out_name   <- file.path(res_dir, sprintf("cif_grid_N%d_p%d_k%d.png", N, p, k))

ggsave(out_name, plot,
       width = 20, height = 12, units = "cm", dpi = 300)


N <- 400; p <- 1000; k <- 20

cif_table <- map_df(1:5,
                    ~ {cif_data_fn(N, p, k, .) |> 
                            mutate(setting = .)})

(plot <- cif_table |> 
        tibble() |> 
        mutate(title = settings_tbl$desc[setting]) |> 
        ggplot(aes(Time, Risk, colour = Method,
                   linewidth = ifelse(Method == "cbSCRIP", 0.5, 0.1))) +
        geom_line() +
        facet_wrap(. ~ title, nrow = 2, scale = "free_x") +
        scale_colour_manual(values = cpl_palette) +
        scale_linewidth_continuous(range = c(0.5, 0.8),
                                   guide = "none") +
        labs(x = "Follow-up time (years)",
             color = "Models",
             y = "Absolute Risk") +
        theme(legend.position = c(0.8, 0.2)))


out_name   <- file.path(res_dir, sprintf("cif_grid_N%d_p%d_k%d.png", N, p, k))

ggsave(out_name, plot,
       width = 20, height = 12, units = "cm", dpi = 300)



N <- 400; p <- 1000; k <- 168

cif_table <- map_df(1:5,
                    ~ {cif_data_fn(N, p, k, .) |> 
                            mutate(setting = .)})

(plot <- cif_table |> 
        tibble() |> 
        mutate(title = settings_tbl$desc[setting]) |> 
        ggplot(aes(Time, Risk, colour = Method,
                   linewidth = ifelse(Method == "cbSCRIP", 0.5, 0.1))) +
        geom_line() +
        facet_wrap(. ~ title, nrow = 2, scale = "free_x") +
        scale_colour_manual(values = cpl_palette) +
        scale_linewidth_continuous(range = c(0.5, 0.8),
                                   guide = "none") +
        labs(x = "Follow-up time (years)",
             color = "Models",
             y = "Absolute Risk") +
        theme(legend.position = c(0.8, 0.2)))


out_name   <- file.path(res_dir, sprintf("cif_grid_N%d_p%d_k%d.png", N, p, k))

ggsave(out_name, plot,
       width = 20, height = 12, units = "cm", dpi = 300)
