#Read data

# run & save (as before)
for (set in 1:5) brier_plot_fn(N, p, k, set)

# load a saved run and tweak the plot label without recomputing
run_id <- sprintf("setting%d_N%d_p%d_k%d", 2, N, p, k)
res <- readRDS(here::here("paper","preds", paste0("results_", run_id, ".rds")))

# edit labels and re-save the figure
res$brier_table$title <- "New label here"
p2 <- ggplot(res$brier_table, aes(times, Brier, colour = model)) +
    geom_line(linewidth = 0.5) +
    facet_grid(. ~ title) +
    scale_colour_manual(values = c("#808080", "#D55E00", "#CC79A7", "#56B4E9", "#E69F00")) +
    theme_bw()

ggsave(here::here("paper","figs", res$meta$fname), p2,
       width = 20, height = 12, units = "cm", dpi = 300)
