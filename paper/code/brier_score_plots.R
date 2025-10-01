source("paper/code/brier_score_fn.R")

#Generate data
i=1
#threshold for beta
thr.b <- 1e-6

#generate title

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

#---------------------
N <- 400; p <- 120; k <- 20

for (set in 1:5){
    brier_plot_fn(N, p, k, set)
}

plot_grid_fn(N, p, k)
dev.off()  

#---------------------
N <- 400; p <- 500; k <- 20

for (set in 1:5){
brier_plot_fn(N, p, k, set)
}

plot_grid_fn(N, p, k)
dev.off()    

#----------------------------
N <- 400; p <- 500; k <- 84

for (set in 1:5){
    brier_plot_fn(N, p, k, set)
}

plot_grid_fn(N, p, k)
dev.off() 

#----------------------------
N <- 400; p <- 1000; k <- 20

for (set in 1:5){
    brier_plot_fn(N, p, k, set)
}

plot_grid_fn(N, p, k)
dev.off()  

#----------------------------------
N <- 400; p <- 1000; k <- 168

for (set in 1:5){
    brier_plot_fn(N, p, k, set)
}

plot_grid_fn(N, p, k)
dev.off()  
