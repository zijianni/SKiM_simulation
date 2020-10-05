###################################
# Summarize results
library(tidyverse)
setwd("~/Google Drive/Hallu/codes/ckgroup/SKIM/")
load("sim_out_10012020.RData")

Res_df <- data.frame(n_A=names(Res)[-length(Res)], 
                     power_AtoB=0, power_BtoC=0, power_AtoC=0,
                     FDR_AtoB=0, FDR_BtoC=0, FDR_AtoC=0)

for(n_a in names(Res)[-length(Res)]){
    res_tmp <- Res[[n_a]]
    
    Res_df$power_BtoC[Res_df$n_A==n_a] <- lapply(
        res_tmp, function(x) x$p_F_CtoB$power
    ) %>% unlist %>% mean
    Res_df$FDR_BtoC[Res_df$n_A==n_a] <- lapply(
        res_tmp, function(x) x$p_F_CtoB$FDR
        ) %>% unlist %>% mean
    Res_df$power_AtoB[Res_df$n_A==n_a] <- lapply(
        res_tmp, function(x) x$p_F_BtoA$power
        ) %>% unlist %>% mean
    Res_df$FDR_AtoB[Res_df$n_A==n_a] <- lapply(
        res_tmp, function(x) x$p_F_BtoA$FDR
        ) %>% unlist %>% mean
    Res_df$power_AtoC[Res_df$n_A==n_a] <- lapply(
        res_tmp, function(x) x$p_F_CtoA$power
        ) %>% unlist %>% mean
    Res_df$FDR_AtoC[Res_df$n_A==n_a] <- lapply(
        res_tmp, function(x) x$p_F_CtoA$FDR
        ) %>% unlist %>% mean
}

View(Res_df)

write.csv(Res_df,file = "sim_out_10012020.csv")
