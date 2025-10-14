################# Plot CIF #########################################

settings_tbl <- tibble::tribble(
    ~setting,    ~desc,
    "setting1",  "Single effects on end-point of interest",
    "setting2", "Single effects on both end-points",
    "setting3", "Opposing effects",
    "setting4", "Mixture of effects",
    "setting5", "Non-proportional hazards"
)

#----------------------------
# Params
img_dir <- here::here("paper","cif_figs")
res_dir <- here::here("paper","figs")


N <- 400; p <- 120; k <- 20

for (set in 1:5){
    cif_plot_fn(N, p, k, set)
}

# grid of figures
files <- file.path(img_dir, sprintf("cif_setting%d_N%d_p%d_k%d.png", 1:5, N, p, k))
out_name   <- file.path(res_dir, sprintf("cif_grid_N%d_p%d_k%d.png", N, p, k))

plot_grid(files,out_name)
dev.off()   
#---------------------
N <- 400; p <- 500; k <- 20

for (set in 1:5){
    cif_plot_fn(N, p, k, set)
}

# grid of figures
files <- file.path(img_dir, sprintf("cif_setting%d_N%d_p%d_k%d.png", 1:5, N, p, k))
out_name   <- file.path(res_dir, sprintf("cif_grid_N%d_p%d_k%d.png", N, p, k))

plot_grid(files,out_name)
dev.off()  

#---------------------
N <- 400; p <- 500; k <- 84

for (set in 1:5){
    cif_plot_fn(N, p, k, set)
}

# grid of figures
files <- file.path(img_dir, sprintf("cif_setting%d_N%d_p%d_k%d.png", 1:5, N, p, k))
out_name   <- file.path(res_dir, sprintf("cif_grid_N%d_p%d_k%d.png", N, p, k))

plot_grid(files,out_name)
dev.off()  

#---------------------
N <- 400; p <- 1000; k <- 20

for (set in 1:5){
    cif_plot_fn(N, p, k, set)
}

# grid of figures
files <- file.path(img_dir, sprintf("cif_setting%d_N%d_p%d_k%d.png", 1:5, N, p, k))
out_name   <- file.path(res_dir, sprintf("cif_grid_N%d_p%d_k%d.png", N, p, k))

plot_grid(files,out_name)
dev.off()  

#---------------------
N <- 400; p <- 1000; k <- 168

for (set in 1:5){
    cif_plot_fn(N, p, k, set)
}

# grid of figures
files <- file.path(img_dir, sprintf("cif_setting%d_N%d_p%d_k%d.png", 1:5, N, p, k))
out_name   <- file.path(res_dir, sprintf("cif_grid_N%d_p%d_k%d.png", N, p, k))

plot_grid(files,out_name)
dev.off()  
