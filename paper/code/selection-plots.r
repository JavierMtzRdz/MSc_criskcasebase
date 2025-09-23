#################################################################
##                       Selection plots                       ##
#################################################################
##
## Description:    
##                 
##
## Author:         Javier Mtz.-Rdz.  
##
## Creation date:  2025-09-23
##
## Email:          javier.mr@stat.ubc.ca
##
## ---------------------------
## Notes:          
## ---------------------------

# Setup ----
## Packages to use ----
pacman::p_load(tidyverse, janitor, writexl, 
                            readxl, scales, mytidyfunctions, 
                            presupuestoR)

## Set theme ------
set_mytheme(text = element_text(family = "Times New Roman"))

# Setup ----
## Packages to use ----

#' To install mytidyfunctions, you need 
#' remotes::install_github("JavierMtzRdz/mytidyfunctions")
if (!require("pacman")) install.packages("pacman")
if (!require("mytidyfunctions")) remotes::install_github("JavierMtzRdz/mytidyfunctions")


pacman::p_load(tidyverse, janitor, writexl, 
               readxl, scales, mytidyfunctions,
               patchwork, here, 
               mtool, bench, kableExtra)

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
library(survminer)


source(here(#"MSc_criskcasebase",
    "paper","code", "fitting_functionsV2.R"))

## Load fonts ----
extrafont::loadfonts(quiet = TRUE)

## Set theme ------
mytidyfunctions::set_mytheme(text = element_text(family = "Times New Roman"))


# Set-up process

save <- F


# Load coefficients

## Load related files
name_files <- list.files(here("paper", "results",
                              "lambda_paths"), 
                         pattern = 'models_',
    recursive = T, full.names = T)

name_model <-  tibble(name_files_comp = name_files,
                      name = str_extract(name_files, "(?<=models_)(.*?)(?=.rds)")) %>% 
    mutate(name = str_replace(name, "selec", "setting"))

## Tidy up dataframe
data_proj <- name_model %>%
    select(name_files_comp, name) %>% 
    separate(name, sep = "_", 
             into = c("setting", "sim", "p", "k"),
             remove = F) %>% 
    mutate(setting = as.numeric(str_remove(setting, "setting-")),
           sim = as.numeric(str_remove(sim, "iter-")),
           p = as.numeric(str_remove(p, "p-")),
           k = as.numeric(str_remove(k, "k-")),
           dim = glue::glue("p = {p}, k = {k}", .sep = "=")) %>% 
    filter(sim <= 15) %>% 
    arrange(p, k)


## Read individual coefficients
if (F){
    
    # select_mod <- map_df(1:nrow(data_proj), 
    #                      function(i){data_proj[i,-1] %>% 
    #                              select(-sim) %>% 
    #                              bind_cols(suppressMessages(readRDS(data_proj$name_files_comp[i])) %>% 
    #                                            mutate(model_size = TP + FP) #%>% 
    #                                        # mutate(model = case_match(model,
    #                                        #                           "iCR.bias" ~ "enet-iCR",
    #                                        #                           "penCRbias" ~ "enet-penCR",
    #                                        #                           "CB-bias" ~ "enet-casebase",
    #                                        #                           "CB-bias-Acc" ~ "enet-casebase-Acc",
    #                                        #                           "postCB-bias" ~ "postenet-casebase",
    #                                        #                           "postCB-bias-Acc" ~ "postenet-casebase-Acc"))
    #                              )})
    # # coefs_mod <- coefs_mod %>% 
    # #     filter(!str_detect(model, "bias"))
    # 
    # saveRDS(select_mod, here("paper", "results", "select_mod.rds"))
}

select_mod <- readRDS(here("paper", "results", "select_mod.rds")) %>% 
    filter(!str_detect(model, "SCAD"))


# Extrapolate for missing model size

sum_selec <- select_mod %>% 
    # filter(sim <= 2) %>%
    group_by(model, setting, dim, sim, p, k, model_size) %>% 
    summarise(n_models = n(),
              across(Sensitivity:MCC, mean)) %>% 
    group_by(model, setting, dim, sim, p, k) %>% 
    group_split() %>%
    map_df(\(y){
        y %>% 
            full_join(tibble(model_size = 1:max(y$p)), 
                      by = join_by(model_size)) %>% 
            arrange(model_size) %>% 
            mutate(across(c(Sensitivity, Specificity, MCC), 
                          \(x){zoo::na.approx(x,
                                              model_size,
                                              na.rm = F)}),
                   across(c(model, setting, dim, sim, p, k),
                          \(x){first(na.omit(x))})) 
    })


sum_selec %>% 
    # filter(sim <= 2,
    #        setting == 1) %>% 
    # filter(p == 1000, k == 168) %>% 
    ungroup() %>% 
    group_by(model_size, p, k, setting, model, dim) %>%
    summarise(Sensitivity = mean(Sensitivity, na.rm = T),
              Specificity = mean(Specificity, na.rm = T)) %>%
    arrange(p, model_size) %>%  
    ggplot(aes(x = model_size, y = Sensitivity,
               color = model)) +
    geom_path() +
    facet_grid(paste0("Setting ", setting) ~fct_inorder(dim),
               scale = "free_x") +
    ylim(-0.01,1.01) +
    labs(x = "Model size",
         color = "Model ")

ggsave(here("paper", "figs", "selection-tpr.png"),
       bg = "transparent",
       width = 200,     
       height = 120,
       units = "mm",
       dpi = 300)


sum_selec %>% 
    # filter(sim <= 2,
    #        setting == 1) %>% 
    # filter(p == 1000, k == 168) %>% 
    ungroup() %>% 
    group_by(model_size, p, k, setting, model, dim) %>% 
    summarise(MCC = mean(MCC, na.rm = T)) %>% 
    arrange(p, model_size) %>%  
    ggplot(aes(x = model_size, y = MCC,
               color = model)) +
    geom_line() +
    facet_grid(paste0("Setting ", setting) ~fct_inorder(dim),
               scale = "free_x") +
    ylim(-0.1,1.01) +
    labs(x = "Model size",
         color = "Model ")

ggsave(here("paper", "figs", "selection-mcc.png"),
       bg = "transparent",
       width = 200,     
       height = 120,
       units = "mm",
       dpi = 300)


# selection curve
sum_selec %>% 
    group_by(model_size, p, k, setting, model, dim) %>% 
    summarise(Sensitivity = mean(Sensitivity),
              Specificity = mean(Specificity)) %>% 
    arrange(p, model_size) %>% 
    ggplot(aes(x = 1 - Specificity, y = Sensitivity,
               color = model)) +
    geom_path() +
    facet_grid(paste0("Setting ", setting) ~fct_inorder(dim)) +
    coord_equal() +
    ylim(-0.01,1.01) +
    labs(color = "Model ")

ggsave(here("paper", "figs", "selection-curve.png"),
       bg = "transparent",
       width = 200,     
       height = 120,
       units = "mm",
       dpi = 300)
