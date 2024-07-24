# Clear environment
rm(list = ls())

# Load required libraries
library(tidyverse)
library(progress)
library(rio)
library(patchwork)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

# Set cohort parameter values 
nreps   <- 1000        # Number of repetitions
nsteps  <- 121         # Number of time steps (Months)
cutoffs <- c(0, 1, 2)  # Threshold for PCF cut-off (Baseline = 2m)  

# Relative infectiousness of subclinical disease
RIs <- c(0.0, 0.125, 0.25, 0.375, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75) 

# Transmission pathways (Table) ----
pb <- progress_bar$new(format = "[:bar] :percent", total = nreps * sum(sapply(RIs, function(ri) if (ri == 0.5) length(cutoffs) else 1)) * 6)

for(q in seq_along(RIs)) {
  sub_transm <- RIs[q]
  sub_folder <- paste0("./data/pathways/RI_", sub_transm)
  if(!dir.exists(sub_folder)) dir.create(sub_folder, recursive = TRUE)
  current_cutoffs <- if (sub_transm == 0.5) cutoffs else 2
  
  for(cutoff in current_cutoffs) {
    cutoff_folder <- if (sub_transm == 0.5) paste0(sub_folder, "/", cutoff, "m") else sub_folder
    if(!dir.exists(cutoff_folder)) dir.create(cutoff_folder, recursive = TRUE)
    
    table <- data.frame(matrix(NA, nrow = 9, ncol = 5))
    
    colnames(table) <- c(
      "Proportion",
      "Transmission contribution overall (%, UI)",
      "Transmission from subclinical (%, UI)",
      "Transmission ≤24 months (%, UI)",
      "Transmission from subclinical ≤24 months (%, UI)")
    
    rownames(table) <- c(
      "Progress", "Regress", "Undulate", 
      "Never PCF", "Pre PCF", "Post PCF",
      "Rapid", "Remote", "Total")
    
    TT_summary <- as.data.frame(matrix(0, nreps, 5))
    TT_sums <- as.data.frame(matrix(0, nreps, 5))
    
    NH_sums <- as.data.frame(matrix(0, nreps, 12))
    NH_summary <- as.data.frame(matrix(0, nreps, 12))
    
    for(i in 1:nreps) {
      load(paste("./data/disease/disease_", i, ".Rdata", sep = ""))
      load(paste(sub_folder, "/RTM_", i, ".Rdata", sep = ""))
      
      # Total pathways
      t_sum  <- sum(RTM, na.rm = TRUE)
      s_sum  <- sum(RTM * (disease == 's'), na.rm = TRUE)
      c_sum  <- sum(RTM * (disease == 'c'), na.rm = TRUE)
      
      # Rapid pathways
      subset_RTM <- RTM[1:25, ]
      t2_sum <- sum(subset_RTM, na.rm = TRUE)
      s2_sum <- sum(subset_RTM * (disease[1:25,] == "s"), na.rm = TRUE)
      c2_sum <- sum(subset_RTM * (disease[1:25,] == "c"), na.rm = TRUE)
      
      # Assign values to TT_sums
      TT_sums[i, 1] <- t_sum
      TT_sums[i, 2] <- s_sum
      TT_sums[i, 3] <- t2_sum
      TT_sums[i, 4] <- s2_sum   
      TT_sums[i, 5] <- c_sum
      
      # Natural history
      load(paste(sub_folder, "/progress_", i, ".Rdata", sep = ""))
      load(paste(sub_folder, "/regress_", i, ".Rdata", sep = ""))
      load(paste(sub_folder, "/undulate_", i, ".Rdata", sep = ""))
      
      subset_progress <- progress[1:25, ]
      subset_regress  <- regress[1:25, ]
      subset_undulate <- undulate[1:25, ]
      
      Pt   <- sum(progress, na.rm = TRUE)
      Ps   <- sum(progress * (disease == 's'), na.rm = TRUE)
      Pt2  <- sum(subset_progress, na.rm = TRUE)
      Rt   <- sum(regress, na.rm = TRUE)
      Rs   <- sum(regress * (disease == 's'), na.rm = TRUE)
      Rt2  <- sum(subset_regress, na.rm = TRUE)
      Ut   <- sum(undulate, na.rm = TRUE)
      Us   <- sum(undulate * (disease == 's'), na.rm = TRUE)
      Ut2  <- sum(subset_undulate, na.rm = TRUE)
      Ps2 <- sum(subset_progress * (disease[1:25,] == "s"), na.rm = TRUE)
      Rs2 <- sum(subset_regress * (disease[1:25,] == "s"), na.rm = TRUE)
      Us2 <- sum(subset_undulate * (disease[1:25,] == "s"), na.rm = TRUE)
      
      # Assign values to NH_sums
      NH_sums[i, 1] <- Pt
      NH_sums[i, 2] <- Ps
      NH_sums[i, 3] <- Pt2
      NH_sums[i, 4] <- Rt
      NH_sums[i, 5] <- Rs 
      NH_sums[i, 6] <- Rt2 
      NH_sums[i, 7] <- Ut
      NH_sums[i, 8] <- Us 
      NH_sums[i, 9] <- Ut2
      NH_sums[i, 10] <- Ps2
      NH_sums[i, 11] <- Rs2
      NH_sums[i, 12] <- Us2
      
      rm(disease, RTM, subset_RTM, progress, regress, undulate, subset_progress, subset_regress, subset_undulate)
      pb$tick()
    }
    
    TT_summary[, c(1,2,3,5)] <- TT_sums[, c(1,2,3,5)] / TT_sums[, 1]
    TT_summary[, 4] <- TT_sums[, 4] / TT_sums[, 3]
    ci <- round(100 * (apply(TT_summary, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)), 1)
    median <- round(100 * (apply(TT_summary, 2, median, na.rm = TRUE)), 1)
    save(TT_summary, file = paste(cutoff_folder, "/TT_summary_", sub_transm, "ri_", cutoff, "m.Rdata", sep = ""))
    
    table["Total", c(1,2)] <- 100
    table["Total", "Transmission from subclinical (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[2], ci[1,2], ci[2,2])
    table["Total", "Transmission ≤24 months (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[3], ci[1,3], ci[2,3])
    table["Total", "Transmission from subclinical ≤24 months (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[4], ci[1,4], ci[2,4])
    
    NH_summary[, c(1,4,7)] <-  NH_sums[, c(1,4,7)] / TT_sums[, 1]
    NH_summary[, 2] <-  NH_sums[, 2] / NH_sums[, 1]
    NH_summary[, 5] <-  NH_sums[, 5] / NH_sums[, 4]
    NH_summary[, 8] <-  NH_sums[, 8] / NH_sums[, 7]
    NH_summary[, c(3,6,9)] <- NH_sums[, c(3,6,9)] / NH_sums[, c(1,4,7)]
    NH_summary[, c(10,11,12)] <- NH_sums[, c(10,11,12)] / NH_sums[, c(3,6,9)]  
    ci <- round(100 * (apply(NH_summary, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)), 1)
    median <- round(100 * (apply(NH_summary, 2, median, na.rm = TRUE)), 1)
    
    table["Progress", "Transmission contribution overall (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[1], ci[1,1], ci[2,1])
    table["Progress", "Transmission from subclinical (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[2], ci[1,2], ci[2,2])
    table["Progress", "Transmission ≤24 months (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[3], ci[1,3], ci[2,3])
    table["Regress", "Transmission contribution overall (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[4], ci[1,4], ci[2,4])
    table["Regress", "Transmission from subclinical (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[5], ci[1,5], ci[2,5])
    table["Regress", "Transmission ≤24 months (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[6], ci[1,6], ci[2,6])
    table["Undulate", "Transmission contribution overall (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[7], ci[1,7], ci[2,7])
    table["Undulate", "Transmission from subclinical (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[8], ci[1,8], ci[2,8])
    table["Undulate", "Transmission ≤24 months (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[9], ci[1,9], ci[2,9])
    table["Progress", "Transmission from subclinical ≤24 months (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[10], ci[1,10], ci[2,10])
    table["Regress", "Transmission from subclinical ≤24 months (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[11], ci[1,11], ci[2,11])
    table["Undulate", "Transmission from subclinical ≤24 months (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[12], ci[1,12], ci[2,12])
    
    rm(ci, median, NH_summary, NH_sums)
    
    # Natural history proportions
    NH_tally <- array(0, dim = c(nreps, 4))
    NH_props <- array(0, dim = c(nreps, 3))
    
    for(i in 1:nreps) {
      load(paste("./data/disease/disease_", i, ".Rdata", sep = ""))
      load(paste(sub_folder, "/progress_", i, ".Rdata", sep = ""))
      load(paste(sub_folder, "/regress_", i, ".Rdata", sep = ""))
      load(paste(sub_folder, "/undulate_", i, ".Rdata", sep = ""))
      
      dis  <- 0
      prog <- 0 
      reg  <- 0
      und  <- 0 
      
      for(j in 1:ncol(disease)) {
        if(length(which(disease[, j] %in% c('c','s'))) > 0) dis <- dis + 1
      }
      
      for(j in 1:ncol(progress)) {
        if(sum(progress[, j], na.rm = TRUE) > 0) prog <- prog + 1
      }
      
      for(j in 1:ncol(regress)) {
        if(sum(regress[, j], na.rm = TRUE) > 0) reg <- reg + 1
      }
      
      for (j in 1:ncol(undulate)) {
        if (sum(undulate[, j], na.rm = TRUE) > 0) und <- und + 1
      }
      
      NH_tally[i, 1] <- dis
      NH_tally[i, 2] <- prog
      NH_tally[i, 3] <- reg
      NH_tally[i, 4] <- und
      
      NH_props[i, 1] <- NH_tally[i, 2] / NH_tally[i, 1]
      NH_props[i, 2] <- NH_tally[i, 3] / NH_tally[i, 1]
      NH_props[i, 3] <- NH_tally[i, 4] / NH_tally[i, 1]
      
      rm(disease, progress, undulate, regress)
      pb$tick()
    }
    
    ci <- round(100 * (apply(NH_props, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)), 1)
    median <- round(100 * (apply(NH_props, 2, median, na.rm = TRUE)), 1)
    
    table["Progress", "Proportion"] <- sprintf("%.1f (%.1f - %.1f)", median[1], ci[1,1], ci[2,1])
    table["Regress", "Proportion"] <- sprintf("%.1f (%.1f - %.1f)", median[2], ci[1,2], ci[2,2])
    table["Undulate", "Proportion"] <- sprintf("%.1f (%.1f - %.1f)", median[3], ci[1,3], ci[2,3])
    
    rm(ci, median, NH_props, NH_tally)
    
    # Passive-case finding
    PCF_summary <- as.data.frame(matrix(0, nreps, 12))
    PCF_sums <- as.data.frame(matrix(0, nreps, 12))
    
    for(i in 1:nreps) {
      load(paste("./data/disease/disease_", i, ".Rdata", sep = ""))
      load(paste(cutoff_folder, "/screen_", i, ".Rdata", sep = ""))
      load(paste(cutoff_folder, "/noscreen_", i, ".Rdata", sep = ""))
      load(paste(cutoff_folder, "/prescreen_", i, ".Rdata", sep = ""))  
      
      subset_screen <- screen[1:25, ]
      subset_noscreen <- noscreen[1:25, ]
      subset_prescreen <- prescreen[1:25, ]
      
      post_t <- sum(screen, na.rm = TRUE) - sum(prescreen, na.rm = TRUE)
      post_s <- sum(screen * (disease == 's'), na.rm = TRUE) - sum(prescreen * (disease == 's'), na.rm = TRUE)
      post_t2 <- sum(subset_screen, na.rm = TRUE) - sum(subset_prescreen, na.rm = TRUE)
      never_t <- sum(noscreen, na.rm = TRUE) 
      never_s <- sum(noscreen * (disease == 's'), na.rm = TRUE) 
      never_t2 <- sum(subset_noscreen, na.rm = TRUE) 
      pre_t <- sum(prescreen, na.rm = TRUE) 
      pre_s <- sum(prescreen * (disease == 's'), na.rm = TRUE) 
      pre_t2 <- sum(subset_prescreen, na.rm = TRUE) 
      post_s2 <- sum(subset_screen * (disease[1:25,] == "s"), na.rm = TRUE) - sum(subset_prescreen * (disease[1:25,] == "s"), na.rm = TRUE)
      never_s2 <- sum(subset_noscreen * (disease[1:25,] == "s"), na.rm = TRUE)
      pre_s2 <- sum(subset_prescreen * (disease[1:25,] == "s"), na.rm = TRUE)
      
      # Assign values to PCF_sums
      PCF_sums[i, 1] <- never_t
      PCF_sums[i, 2] <- never_s
      PCF_sums[i, 3] <- never_t2
      PCF_sums[i, 4] <- pre_t
      PCF_sums[i, 5] <- pre_s
      PCF_sums[i, 6] <- pre_t2
      PCF_sums[i, 7] <- post_t
      PCF_sums[i, 8] <- post_s
      PCF_sums[i, 9] <- post_t2
      PCF_sums[i, 10] <- never_s2
      PCF_sums[i, 11] <- pre_s2
      PCF_sums[i, 12] <- post_s2
      
      rm(disease, prescreen, noscreen, screen, subset_prescreen, subset_noscreen, subset_screen)
      pb$tick()
    }
    
    PCF_summary[, c(1,4,7)] <- PCF_sums[, c(1,4,7)] / TT_sums[, 1]
    PCF_summary[, 2] <- PCF_sums[, 2] / PCF_sums[, 1]
    PCF_summary[, 5] <- PCF_sums[, 5] / PCF_sums[, 4]
    PCF_summary[, 8] <- PCF_sums[, 8] / PCF_sums[, 7]
    PCF_summary[, c(3,6,9)] <- PCF_sums[, c(3,6,9)] / PCF_sums[, c(1,4,7)]
    PCF_summary[, c(10,11,12)] <- PCF_sums[, c(10,11,12)] / PCF_sums[, c(3,6,9)]
    ci <- round(100 * (apply(PCF_summary, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)), 1)
    median <- round(100 * (apply(PCF_summary, 2, median, na.rm = TRUE)), 1)
    save(PCF_summary, file = paste(cutoff_folder, "/PCF_summary_", sub_transm, "ri_", cutoff, "m.Rdata", sep = ""))
    
    table["Never PCF", "Transmission contribution overall (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[1], ci[1,1], ci[2,1])
    table["Never PCF", "Transmission from subclinical (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[2], ci[1,2], ci[2,2])
    table["Never PCF", "Transmission ≤24 months (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[3], ci[1,3], ci[2,3])
    table["Pre PCF", "Transmission contribution overall (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[4], ci[1,4], ci[2,4])
    table["Pre PCF", "Transmission from subclinical (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[5], ci[1,5], ci[2,5])
    table["Pre PCF", "Transmission ≤24 months (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[6], ci[1,6], ci[2,6])
    table["Post PCF", "Transmission contribution overall (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[7], ci[1,7], ci[2,7])
    table["Post PCF", "Transmission from subclinical (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[8], ci[1,8], ci[2,8])
    table["Post PCF", "Transmission ≤24 months (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[9], ci[1,9], ci[2,9])
    table["Never PCF", "Transmission from subclinical ≤24 months (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[10], ci[1,10], ci[2,10])
    table["Pre PCF", "Transmission from subclinical ≤24 months (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[11], ci[1,11], ci[2,11])
    table["Post PCF", "Transmission from subclinical ≤24 months (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[12], ci[1,12], ci[2,12])
    
    rm(ci, median, PCF_summary, PCF_sums)
    
    # Passive-case finding proportions
    PCF_tally <- array(0, dim = c(nreps, 4))
    PCF_props <- array(0, dim = c(nreps, 3))
    
    for(i in 1:nreps) {
      load(paste("./data/disease/disease_", i, ".Rdata", sep = ""))
      load(paste(cutoff_folder, "/screen_", i, ".Rdata", sep = ""))
      load(paste(cutoff_folder, "/noscreen_", i, ".Rdata", sep = ""))
      load(paste(cutoff_folder, "/prescreen_", i, ".Rdata", sep = ""))
      
      dis <- 0
      never <- 0 
      pre <- 0
      post <- 0 
      
      for(j in 1:ncol(disease)) {
        if(length(which(disease[, j] %in% c('c', 's'))) > 0) dis <- dis + 1
      }
      
      for(j in 1:ncol(noscreen)) {
        if(sum(noscreen[, j], na.rm = TRUE) > 0) never <- never + 1
      }
      
      for(j in 1:ncol(prescreen)) {
        if(sum(prescreen[, j], na.rm = TRUE) > 0) pre <- pre + 1
      }
      
      for(j in 1:ncol(screen)) {
        if(sum(screen[, j], na.rm = TRUE) > 0) post <- post + 1
      }
      
      PCF_tally[i, 1] <- dis
      PCF_tally[i, 2] <- never
      PCF_tally[i, 3] <- pre
      PCF_tally[i, 4] <- post
      
      PCF_props[i, 1] <- PCF_tally[i, 2] / PCF_tally[i, 1]
      PCF_props[i, 2] <- PCF_tally[i, 3] / PCF_tally[i, 1]
      PCF_props[i, 3] <- PCF_tally[i, 4] / PCF_tally[i, 1]
      
      rm(disease, prescreen, screen, noscreen)
      pb$tick()
    }
    
    ci <- round(100 * (apply(PCF_props, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)), 1)
    median <- round(100 * (apply(PCF_props, 2, median, na.rm = TRUE)), 1)
    
    table["Never PCF", "Proportion"] <- sprintf("%.1f (%.1f - %.1f)", median[1], ci[1,1], ci[2,1])
    table["Pre PCF", "Proportion"] <- sprintf("%.1f (%.1f - %.1f)", median[2], ci[1,2], ci[2,2])
    table["Post PCF", "Proportion"] <- sprintf("%.1f (%.1f - %.1f)", median[3], ci[1,3], ci[2,3])
    
    rm(ci, median, PCF_props, PCF_tally)
    
    # Rapid/remote progression
    RR_sums <- as.data.frame(matrix(0, nreps, 9))
    RR_summary <- as.data.frame(matrix(0, nreps, 9))
    
    for (i in 1:nreps) {
      load(paste("./data/disease/disease_", i, ".Rdata", sep = ""))
      load(paste(sub_folder, "/remote_", i, ".Rdata", sep = ""))
      load(paste(sub_folder, "/rapid_", i, ".Rdata", sep = ""))
      
      subset_rapid <- rapid[1:25, ]
      subset_remote <- remote[1:25, ]
      subset_rapid2 <- rapid[26:121, ]
      
      rapidt <- sum(rapid, na.rm = TRUE)
      rapids <- sum(rapid * (disease == 's'), na.rm = TRUE)
      rapidt2 <- sum(subset_rapid, na.rm = TRUE)
      remotet <- sum(remote, na.rm = TRUE)
      remotes <- sum(remote * (disease == 's'), na.rm = TRUE)
      remotet2 <- sum(subset_remote, na.rm = TRUE)
      rapids2 <- sum(subset_rapid * (disease[1:25,] == "s"), na.rm = TRUE)
      remotes2 <- sum(subset_remote * (disease[1:25,] == "s"), na.rm = TRUE)
      rapidt8 <- sum(subset_rapid2, na.rm = TRUE)
      
      # Assign values to RR_sums
      RR_sums[i, 1] <- rapidt
      RR_sums[i, 2] <- rapids
      RR_sums[i, 3] <- rapidt2
      RR_sums[i, 4] <- remotet
      RR_sums[i, 5] <- remotes
      RR_sums[i, 6] <- remotet2
      RR_sums[i, 7] <- rapids2
      RR_sums[i, 8] <- remotes2
      RR_sums[i, 9] <- rapidt8
      
      rm(disease, rapid, remote, subset_rapid, subset_remote, subset_rapid2)
      pb$tick()
    }
    
    RR_summary[, c(1,4)] <-  RR_sums[, c(1,4)] / TT_sums[, 1]
    RR_summary[, 1] <-  RR_sums[, 1] / TT_sums[, 1]
    RR_summary[, 2] <-  RR_sums[, 2] / RR_sums[, 1]
    RR_summary[, 5] <-  RR_sums[, 5] / RR_sums[, 4]
    RR_summary[, c(3,6)] <- RR_sums[, c(3,6)] / RR_sums[, c(1,4)]
    RR_summary[, c(7,8)] <- RR_sums[, c(7)] / RR_sums[, c(3)]  
    RR_summary[, 9] <-  RR_sums[, 9] / TT_sums[, 1]
    ci <- round(100 * (apply(RR_summary, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)), 1)
    median <- round(100 * (apply(RR_summary, 2, median, na.rm = TRUE)), 1)
    
    table["Rapid", "Transmission contribution overall (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[1], ci[1,1], ci[2,1])
    table["Rapid", "Transmission from subclinical (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[2], ci[1,2], ci[2,2])
    table["Rapid", "Transmission ≤24 months (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[3], ci[1,3], ci[2,3])
    table["Remote", "Transmission contribution overall (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[4], ci[1,4], ci[2,4])
    table["Remote", "Transmission from subclinical (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[5], ci[1,5], ci[2,5])
    table["Remote", "Transmission ≤24 months (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[6], ci[1,6], ci[2,6])
    table["Rapid", "Transmission from subclinical ≤24 months (%, UI)"] <- sprintf("%.1f (%.1f - %.1f)", median[7], ci[1,7], ci[2,7])
    
    rm(ci, median, RR_summary, RR_sums)
    
    # Rapid/remote progression proportion
    RR_tally <- array(0, dim = c(nreps, 4))
    RR_props <- array(0, dim = c(nreps, 3))
    
    for(i in 1:nreps) {
      load(paste("./data/disease/disease_", i, ".Rdata", sep = ""))
      load(paste(sub_folder, "/rapid_", i, ".Rdata", sep = ""))
      load(paste(sub_folder, "/remote_", i, ".Rdata", sep = ""))
      
      dis <- 0
      rap <- 0 
      rem <- 0
      
      for(j in 1:ncol(disease)) {
        if(length(which(disease[, j] %in% c('c', 's'))) > 0) dis <- dis + 1
      }
      
      for(j in 1:ncol(rapid)) {
        if(sum(rapid[, j], na.rm = TRUE) > 0) rap <- rap + 1
      }
      
      for (j in 1:ncol(remote)) {
        if (sum(remote[, j], na.rm = TRUE) > 0) rem <- rem + 1
      }
      
      RR_tally[i, 1] <- dis
      RR_tally[i, 2] <- rap
      RR_tally[i, 3] <- rem
      
      RR_props[i, 1] <- RR_tally[i, 2] / RR_tally[i, 1]
      RR_props[i, 2] <- RR_tally[i, 3] / RR_tally[i, 1]
      RR_props[i, 3] <- RR_tally[i, 4] / RR_tally[i, 1]
      
      rm(disease, rapid, remote)
      pb$tick()
    }
    
    ci <- round(100 * (apply(RR_props, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)), 1)
    median <- round(100 * (apply(RR_props, 2, median, na.rm = TRUE)), 1)
    
    table["Rapid", "Proportion"] <- sprintf("%.1f (%.1f - %.1f)", median[1], ci[1,1], ci[2,1])
    table["Remote", "Proportion"] <- sprintf("%.1f (%.1f - %.1f)", median[2], ci[1,2], ci[2,2])
    
    rm(ci, median, RR_props, RR_tally, TT_summary, TT_sums)
    
    if(!dir.exists("results")) dir.create("results", recursive = TRUE)
    save(table, file = paste("results/table_", sub_transm, "ri_", cutoff, "m.Rdata", sep = ""))
    rm(list=ls()[! ls() %in% c("nreps", "nsteps", "cutoffs", "RIs", "pb", "sub_transm", "sub_folder")])
  }
}

rm(pb, sub_folder, sub_transm)

# Relative infectiousness of subclinical TB ----
RI_transm <- data.frame(matrix(ncol = 7, nrow = length(RIs)))
rownames(RI_transm) <- seq_along(RIs)
colnames(RI_transm) <- c("S_Med", "S_Low", "S_High","C_Med", "C_Low", "C_High", "Sub_RI")

for (q in seq_along(RIs)) {
  sub_transm <- RIs[q]
  sub_folder <- paste0("./data/pathways/RI_", sub_transm)
  file_folder <- if(sub_transm == 0.5) paste0(sub_folder, "/2m") else sub_folder
  load(paste(file_folder, "/TT_summary_", sub_transm, "ri_2m.Rdata", sep = ""))
  
  ci <- apply(TT_summary, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  median <- apply(TT_summary, 2, median, na.rm = TRUE)
  
  RI_transm[q, "S_Med"] <- median[2]
  RI_transm[q, "S_Low"] <- ci[1, 2]
  RI_transm[q, "S_High"] <- ci[2, 2]
  RI_transm[q, "C_Med"] <- median[5]
  RI_transm[q, "C_Low"] <- ci[1, 5]
  RI_transm[q, "C_High"] <- ci[2, 5]
  RI_transm[q, "Sub_RI"] <- sub_transm

  rm(list=ls()[! ls() %in% c("nreps", "nsteps", "RIs", "RI_transm", "PCF_transm")])
}

if(!dir.exists("figures")) dir.create("figures", recursive = TRUE)
fig_RI_transm <- ggplot(RI_transm, aes(x = Sub_RI)) +
  geom_point(aes(y = S_Med, colour = "Subclinical transmission"), size = 1) +
  geom_point(aes(y = C_Med, colour = "Clinical transmission"), size = 1) +
  geom_ribbon(aes(ymin = S_Low, ymax = S_High, fill = "Subclinical transmission"), alpha = 0.3) +
  geom_ribbon(aes(ymin = C_Low, ymax = C_High, fill = "Clinical transmission"), alpha = 0.3) +
  geom_vline(xintercept = 0.5, colour = "red") +
  scale_color_manual(values = c("Subclinical transmission" = "#00BFC4", "Clinical transmission" = "#F8766D")) +
  scale_fill_manual(values = c("Subclinical transmission" = "#00BFC4", "Clinical transmission" = "#F8766D")) +
  labs(x = "Relative transmission rate of subclinical TB", y = "Proportion of all transmission", 
       colour = "State", fill = "State") +
  theme_bw() +
  theme(legend.key.size = unit(1.2, "lines"), legend.position = "bottom",
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), panel.border = element_blank(), 
        axis.line = element_line(colour = 'black', size = 0.5), axis.ticks = element_line(colour = 'black', size = 0.5))
save(RI_transm, file="results/RI_transm.Rdata")
png("figures/fig_RI_transm.png",width=600, height=400)
print(fig_RI_transm)
dev.off()

# Passive-case finding transmission ----
PCF_transm <- data.frame(matrix(ncol = 7, nrow = length(RIs)))
rownames(PCF_transm) <- seq_along(RIs)
colnames(PCF_transm) <- c("E_Med", "E_Low", "E_High","I_Med", "I_Low", "I_High", "Sub_RI")

for (q in seq_along(RIs)) {
  sub_transm <- RIs[q]
  sub_folder <- paste0("./data/pathways/RI_", sub_transm)
  file_folder <- if(sub_transm == 0.5) paste0(sub_folder, "/2m") else sub_folder
  load(paste(file_folder, "/PCF_summary_", sub_transm, "ri_2m.Rdata", sep = ""))
  
  PCF_summary[, 13] <- PCF_summary[, 1] + PCF_summary[, 4]
  ci <- apply(PCF_summary, 2, quantile, probs = c(0.025,0.975), na.rm = TRUE)
  median <- apply(PCF_summary, 2, median, na.rm = TRUE)
  
  PCF_transm[q, "E_Med"] <- median[7]
  PCF_transm[q, "E_Low"] <- ci[1, 7]
  PCF_transm[q, "E_High"] <- ci[2, 7]
  PCF_transm[q, "I_Med"] <- median[13]
  PCF_transm[q, "I_Low"] <- ci[1, 13]
  PCF_transm[q, "I_High"] <- ci[2, 13]
  PCF_transm[q, "Sub_RI"] <- sub_transm 
  
  rm(list=ls()[! ls() %in% c("nreps", "nsteps", "RIs", "RI_transm", "PCF_transm")])
}

fig_PCF_transm <- ggplot(PCF_transm, aes(x = Sub_RI)) +
  geom_point(aes(y = E_Med, colour = "PCF-elligible transmission"), size = 1) +
  geom_point(aes(y = I_Med, colour = "PCF-inelligible transmission"), size = 1) +
  geom_ribbon(aes(ymin = E_Low, ymax = E_High, fill = "PCF-elligible transmission"), alpha = 0.3) +
  geom_ribbon(aes(ymin = I_Low, ymax = I_High, fill = "PCF-inelligible transmission"), alpha = 0.3) +
  geom_vline(xintercept = 0.5, colour = "red") +
  scale_color_manual(values = c("PCF-elligible transmission" = "#629DFF", "PCF-inelligible transmission" = "#FFA500")) +
  scale_fill_manual(values = c("PCF-elligible transmission" = "#629DFF", "PCF-inelligible transmission" = "#FFA500")) +
  labs(x = "Relative transmission rate of subclinical TB", y = "Proportion of all transmission", 
       colour = "PCF pathway", fill = "PCF pathway") +
  theme_bw() +
  theme(legend.key.size = unit(1.2, "lines"), legend.position = "bottom",
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), panel.border = element_blank(), 
        axis.line = element_line(colour = 'black', size = 0.5), axis.ticks = element_line(colour = 'black', size = 0.5))
save(PCF_transm, file="results/PCF_transm.Rdata")
png("figures/fig_PCF_transm.png",width=600, height=400)
print(fig_RI_transm)
dev.off()

# Variable PCF thresholds ----
m0 <- import("./results/table_0.5ri_0m.Rdata", trust = TRUE)[4:6, 2:3] %>% 
  mutate(across(everything(), ~ as.numeric(sub(" .*", "", .)) / 100)) %>% 
  setNames(c("total", "propsctb")) %>% 
  rownames_to_column(var = "pcf") %>% 
  mutate(pcf = case_when(pcf == "Never PCF" ~ "never",
                         pcf == "Pre PCF" ~ "pre",
                         pcf == "Post PCF" ~ "post")) %>% 
  mutate(sctb = total * propsctb, ctb = total * (1 - propsctb)) %>% 
  mutate(type = "0m") %>% 
  select(type, pcf, sctb, ctb)
  
m1 <- import("./results/table_0.5ri_1m.Rdata", trust = TRUE)[4:6, 2:3] %>% 
  mutate(across(everything(), ~ as.numeric(sub(" .*", "", .)) / 100)) %>% 
  setNames(c("total", "propsctb")) %>% 
  rownames_to_column(var = "pcf") %>% 
  mutate(pcf = case_when(pcf == "Never PCF" ~ "never",
                         pcf == "Pre PCF" ~ "pre",
                         pcf == "Post PCF" ~ "post")) %>% 
  mutate(sctb = total * propsctb, ctb = total * (1 - propsctb)) %>% 
  mutate(type = "1m") %>% 
  select(type, pcf, sctb, ctb)

m2 <- import("./results/table_0.5ri_2m.Rdata", trust = TRUE)[4:6, 2:3] %>% 
  mutate(across(everything(), ~ as.numeric(sub(" .*", "", .)) / 100)) %>% 
  setNames(c("total", "propsctb")) %>% 
  rownames_to_column(var = "pcf") %>% 
  mutate(pcf = case_when(pcf == "Never PCF" ~ "never",
                         pcf == "Pre PCF" ~ "pre",
                         pcf == "Post PCF" ~ "post")) %>% 
  mutate(sctb = total * propsctb, ctb = total * (1 - propsctb)) %>% 
  mutate(type = "2m") %>% 
  select(type, pcf, sctb, ctb)

PCFsa <- import("./results/table_0.5ri_2m.Rdata", trust = TRUE)[9, 2:3] %>% 
  mutate(across(everything(), ~ as.numeric(sub(" .*", "", .)) / 100)) %>% 
  setNames(c("total", "propsctb")) %>% 
  rownames_to_column(var = "pcf") %>% 
  mutate(sctb = total * propsctb, ctb = total * (1 - propsctb)) %>% 
  mutate(pcf = "total", type = "all") %>% 
  select(type, pcf, sctb, ctb) %>% 
  rbind(m0, m1, m2) %>% 
  pivot_longer(cols = -c(type, pcf), names_to = "tb", values_to = "val") %>% 
  mutate(tb = factor(tb, levels = c("ctb", "sctb"), labels = c("Clinical", "Subclinical"))) %>% 
  mutate(pcf = factor(pcf, levels = c("never", "pre", "post", "total"),
                      labels = c("Never PCF Eligible", "Pre-PCF Eligible", "PCF Eligible", "Total"))) %>% 
  mutate(type = factor(type, levels = c("0m", "1m", "2m", "all"), 
                       labels = c("0 Months", "1 Month", "2 Months (Baseline)", "Total"))) %>% 
  arrange(type, pcf)

f1 <- ggplot(filter(PCFsa, pcf != "Total")) +
  facet_wrap(~type, scales = 'free_x', nrow = 1) + 
  geom_col(mapping = aes(x = pcf, y = val, fill = tb), position = "stack") + 
  labs(y = 'Percentage of attributable transmission', fill = 'State') +
  scale_y_continuous(labels = scales::label_percent(scale = 100, suffix = '%')) +
  coord_cartesian(ylim = c(0,1)) + 
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        legend.position = "none") +
  theme(strip.text.x = element_text(angle = 0),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm"))

f2 <- ggplot(filter(PCFsa, pcf == "Total")) +
  facet_wrap(~type, scales = 'free_x', nrow = 1) + 
  geom_col(mapping = aes(x = pcf, y = val, fill = tb), position = "stack") + 
  labs(y = 'Percentage of attributable transmission', fill = 'State') +
  scale_y_continuous(labels = scales::label_percent(scale = 100, suffix = '%')) +
  coord_cartesian(ylim = c(0,1)) + 
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y.left = element_blank()) +
  theme(strip.text.x = element_text(angle = 0),
        plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "mm"))

f3 <- (f1 + f2) +
  plot_layout(widths = c(8, 1), ncol = 2)
png("figures/fig_PCF_thresholds.png",width=600, height=400)
print(f3)
dev.off()

rm(list=ls()[! ls() %in% c("nreps", "nsteps", "RIs", "PCFsa", "f3")])

# Monthly PAF by disease state ----
mT <- as.data.frame(matrix(0, nrow = nsteps, ncol = nreps))
mC <- as.data.frame(matrix(0, nrow = nsteps, ncol = nreps))
mS <- as.data.frame(matrix(0, nrow = nsteps, ncol = nreps))

RI <- 0.5
pb <- progress_bar$new(format = "[:bar] :percent", total = nreps)

for (i in 1:nreps) {
  load(paste("./data/pathways/RI_", RI,"/RTM_", i, ".Rdata", sep = ""))
  
  for (j in 1:nsteps){
    mT[j, i] <- as.numeric(sum(RTM[j, ]))
    mC[j, i] <- as.numeric(sum(RTM[j, RTM[j, ] == 1]))
    mS[j, i] <- as.numeric(sum(RTM[j, RTM[j, ] < 1]))
  }
  
  rm("RTM")
  pb$tick()
}

mC <- mC %>%
  rowwise() %>%
  mutate(medC = median(c_across(everything()))) %>%
  ungroup() %>% 
  select(medC)

mS <- mS %>%
  rowwise() %>%
  mutate(medS = median(c_across(everything()))) %>%
  ungroup() %>% 
  select(medS)

RTM <- mT %>%
  rowwise() %>%
  mutate(medT = median(c_across(everything()))) %>%
  ungroup() %>%
  select(medT) %>%
  cbind(mC, mS) %>% 
  mutate(time = row_number()) %>% 
  mutate(propS = medS / (medS + medC), propC = medC / (medS + medC)) %>% 
  mutate(propS = ifelse(is.na(propS), 0, propS), propC = ifelse(is.na(propC), 0, propC)) %>% 
  pivot_longer(cols = starts_with("prop"), names_to = "state", values_to = "pct") %>%
  select(time, state, pct) %>%
  mutate(state = ifelse(state == "propS", "Subclinical", "Clinical"))

fig_RTM <- ggplot(data = RTM, aes(x = time, y = pct, fill = state)) +
  geom_area(position = "stack") +
  labs(x = "ms following infection", y = "Proportion of attributable transmission", fill = "Disease state") +
  theme_bw() +
  theme(legend.key.size = unit(1.2, "lines"), legend.position = "bottom", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), panel.border = element_blank(), 
        axis.line = element_line(colour = 'black', size = 0.5), axis.ticks = element_line(colour = 'black', size = 0.5))
png("figures/fig_RTM.png",width=600, height=400)
print(fig_RTM)
dev.off()

RTMd <- mT %>%
  rowwise() %>%
  mutate(medT = median(c_across(everything()))) %>%
  ungroup() %>%
  select(medT) %>%
  cbind(mC, mS) %>% 
  mutate(time = row_number()) %>% 
  mutate(propS = medS / (sum(medT)), propC = medC / (sum(medT))) %>% 
  mutate(propS = ifelse(is.na(propS), 0, propS), propC = ifelse(is.na(propC), 0, propC)) %>% 
  pivot_longer(cols = starts_with("prop"), names_to = "state", values_to = "pct") %>%
  select(time, state, pct) %>%
  mutate(state = ifelse(state == "propS", "Subclinical", "Clinical"))

fig_RTMd <- ggplot(data = RTMd, aes(x = time, y = pct, fill = state)) +
  geom_area(position = "stack") +
  labs(x = "ms following infection", y = "Proportion of attributable transmission", fill = "Disease state") +
  theme_bw() +
  theme(legend.key.size = unit(1.2, "lines"), legend.position = "bottom", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), panel.border = element_blank(), 
        axis.line = element_line(colour = 'black', size = 0.5), axis.ticks = element_line(colour = 'black', size = 0.5))
png("figures/fig_RTMd.png",width=600, height=400)
print(fig_RTMd)
dev.off()

rm(list=ls()[! ls() %in% c("nreps", "nsteps", "RIs", "PCFsa", "f3", "RTM", "RTMd")])

# Monthly PAF by PCF eligibility ----
mPre <- as.data.frame(matrix(0, nrow = nsteps, ncol = nreps))
mPost <- as.data.frame(matrix(0, nrow = nsteps, ncol = nreps))
mNever <- as.data.frame(matrix(0, nrow = nsteps, ncol = nreps))
mTotal <- as.data.frame(matrix(0, nrow = nsteps, ncol = nreps))

RI <- 0.5
pb <- progress_bar$new(format = "[:bar] :percent", total = nreps)

for (i in 1:nreps) {
  load(paste("./data/pathways/RI_", RI,"/2m/screen_", i, ".Rdata", sep = ""))
  load(paste("./data/pathways/RI_", RI,"/2m/noscreen_", i, ".Rdata", sep = ""))
  load(paste("./data/pathways/RI_", RI,"/2m/prescreen_", i, ".Rdata", sep = ""))
  
  
  for (j in 1:nsteps){
    mPre[j, i] <- (sum(prescreen[j, ]))
    mPost[j, i] <- (sum(screen[j, ])) - (sum(prescreen[j, ]))
    mNever[j, i] <- (sum(noscreen[j, ]))
    mTotal[j, i] <- ((mPre[j, i]) + (mPost[j, i]) + (mNever[j, i]))
  }
  
  rm("screen", "noscreen", "prescreen")
  pb$tick()
}

mPre <- mPre %>%
  rowwise() %>%
  mutate(medPre = median(c_across(everything()))) %>%
  ungroup() %>% 
  select(medPre)

mPost <- mPost %>%
  rowwise() %>%
  mutate(medPost = median(c_across(everything()))) %>%
  ungroup() %>% 
  select(medPost)

mNever <- mNever %>%
  rowwise() %>%
  mutate(medNever = median(c_across(everything()))) %>%
  ungroup() %>% 
  select(medNever)

PCF <- mTotal %>%
  rowwise() %>%
  mutate(medTotal = median(c_across(everything()))) %>%
  ungroup() %>%
  select(medTotal) %>%
  cbind(mPre, mPost, mNever) %>% 
  mutate(time = row_number()) %>% 
  mutate(propPre = medPre / medTotal, 
         propPost = medPost / medTotal, 
         propNever = medNever / medTotal) %>% 
  mutate(propPre = ifelse(is.na(propPre), 0, propPre),
         propPost = ifelse(is.na(propPost), 0, propPost),
         propNever = ifelse(is.na(propNever), 0, propNever)) %>% 
  pivot_longer(cols = starts_with("prop"), names_to = "pcf", values_to = "pct") %>%
  select(time, pcf, pct) %>%
  mutate(pcf = case_when(pcf == "propPre" ~ "Pre-PCF eligible",
                         pcf == "propPost" ~ "Post-PCF eligible",
                         pcf == "propNever" ~ "Never PCF eligible")) %>% 
  mutate(pcf = factor(pcf, levels = c("Never PCF eligible", "Pre-PCF eligible", "Post-PCF eligible")))

fig_PCF <- ggplot(data = PCF, aes(x = time, y = pct, fill = pcf)) +
  geom_area(position = "stack") +
  labs(x = "Months following infection", y = "Proportion of attributable transmission", 
       fill = "PCF eligibility") +
  theme_bw() +
  theme(legend.key.size = unit(1.2, "lines"), legend.position = "bottom", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), panel.border = element_blank(), 
        axis.line = element_line(colour = 'black', size = 0.5), axis.ticks = element_line(colour = 'black', size = 0.5))

png("figures/fig_PCF.png",width=600, height=400)
print(fig_PCF)
dev.off()

PCFd <- mTotal %>%
  rowwise() %>%
  mutate(medTotal = median(c_across(everything()))) %>%
  ungroup() %>%
  select(medTotal) %>%
  cbind(mPre, mPost, mNever) %>% 
  mutate(time = row_number()) %>% 
  mutate(propPre = medPre / sum(medTotal), 
         propPost = medPost / sum(medTotal), 
         propNever = medNever / sum(medTotal)) %>% 
  mutate(propPre = ifelse(is.na(propPre), 0, propPre),
         propPost = ifelse(is.na(propPost), 0, propPost),
         propNever = ifelse(is.na(propNever), 0, propNever)) %>% 
  pivot_longer(cols = starts_with("prop"), names_to = "pcf", values_to = "pct") %>%
  select(time, pcf, pct) %>%
  mutate(pcf = case_when(pcf == "propPre" ~ "Pre-PCF eligible",
                         pcf == "propPost" ~ "Post-PCF eligible",
                         pcf == "propNever" ~ "Never PCF eligible")) %>% 
  mutate(pcf = factor(pcf, levels = c("Never PCF eligible", "Pre-PCF eligible", "Post-PCF eligible")))

fig_PCFd <- ggplot(data = PCFd, aes(x = time, y = pct, fill = pcf)) +
  geom_area(position = "stack") +
  labs(x = "Months following infection", y = "Proportion of attributable transmission", 
       fill = "PCF eligibility") +
  theme_bw() +
  theme(legend.key.size = unit(1.2, "lines"), legend.position = "bottom", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), panel.border = element_blank(), 
        axis.line = element_line(colour = 'black', size = 0.5), axis.ticks = element_line(colour = 'black', size = 0.5))

png("figures/fig_PCFd.png",width=600, height=400)
print(fig_PCFd)
dev.off()

rm(list=ls()[! ls() %in% c("nreps", "nsteps", "RIs", "PCFsa", "f3", "RTM", "RTMd", "PCF", "PCFd")])

# Monthly PAF by NH pathway ----
mPro <- as.data.frame(matrix(0, nrow = nsteps, ncol = nreps))
mReg <- as.data.frame(matrix(0, nrow = nsteps, ncol = nreps))
mUnd <- as.data.frame(matrix(0, nrow = nsteps, ncol = nreps))
mTotal <- as.data.frame(matrix(0, nrow = nsteps, ncol = nreps))

RI <- 0.5
pb <- progress_bar$new(format = "[:bar] :percent", total = nreps)

for (i in 1:nreps) {
  load(paste("./data/pathways/RI_", RI,"/progress_", i, ".Rdata", sep = ""))
  load(paste("./data/pathways/RI_", RI,"/regress_", i, ".Rdata", sep = ""))
  load(paste("./data/pathways/RI_", RI,"/undulate_", i, ".Rdata", sep = ""))
  
  for (j in 1:nsteps){
    mPro[j, i] <- (sum(progress[j, ]))
    mReg[j, i] <- (sum(regress[j, ]))
    mUnd[j, i] <- (sum(undulate[j, ]))
    mTotal[j, i] <- ((mPro[j, i]) + (mReg[j, i]) + (mUnd[j, i]))
  }
  
  rm("progress", "regress", "undulate")
  pb$tick()
}

mPro <- mPro %>%
  rowwise() %>%
  mutate(medPro = median(c_across(everything()))) %>%
  ungroup() %>% 
  select(medPro)

mReg <- mReg %>%
  rowwise() %>%
  mutate(medReg = median(c_across(everything()))) %>%
  ungroup() %>% 
  select(medReg)

mUnd <- mUnd %>%
  rowwise() %>%
  mutate(medUnd = median(c_across(everything()))) %>%
  ungroup() %>% 
  select(medUnd)

NH <- mTotal %>%
  rowwise() %>%
  mutate(medTotal = median(c_across(everything()))) %>%
  ungroup() %>%
  select(medTotal) %>%
  cbind(mPro, mReg, mUnd) %>% 
  mutate(time = row_number()) %>% 
  mutate(propPro = medPro / medTotal, 
         propReg = medReg / medTotal, 
         propUnd = medUnd / medTotal) %>% 
  mutate(propPro = ifelse(is.na(propPro), 0, propPro),
         propReg = ifelse(is.na(propReg), 0, propReg),
         propUnd = ifelse(is.na(propUnd), 0, propUnd)) %>% 
  pivot_longer(cols = starts_with("prop"), names_to = "pcf", values_to = "pct") %>%
  select(time, pcf, pct) %>%
  mutate(pcf = case_when(pcf == "propPro" ~ "Progression",
                         pcf == "propReg" ~ "Regression",
                         pcf == "propUnd" ~ "Undulation")) %>% 
  mutate(pcf = factor(pcf, levels = c("Progression", "Regression", "Undulation")))

fig_NH <- ggplot(data = NH, aes(x = time, y = pct, fill = pcf)) +
  geom_area(position = "stack") +
  labs(x = "Months following infection", y = "Proportion of attributable transmission", 
       fill = "Natural history pathways") +
  theme_bw() +
  theme(legend.key.size = unit(1.2, "lines"), legend.position = "bottom", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), panel.border = element_blank(), 
        axis.line = element_line(colour = 'black', size = 0.5), axis.ticks = element_line(colour = 'black', size = 0.5))

png("figures/fig_NH.png",width=600, height=400)
print(fig_NH)
dev.off()

NHd <- mTotal %>%
  rowwise() %>%
  mutate(medTotal = median(c_across(everything()))) %>%
  ungroup() %>%
  select(medTotal) %>%
  cbind(mPro, mReg, mUnd) %>% 
  mutate(time = row_number()) %>% 
  mutate(propPro = medPro / sum(medTotal), 
         propReg = medReg / sum(medTotal), 
         propUnd = medUnd / sum(medTotal)) %>% 
  mutate(propPro = ifelse(is.na(propPro), 0, propPro),
         propReg = ifelse(is.na(propReg), 0, propReg),
         propUnd = ifelse(is.na(propUnd), 0, propUnd)) %>% 
  pivot_longer(cols = starts_with("prop"), names_to = "pcf", values_to = "pct") %>%
  select(time, pcf, pct) %>%
  mutate(pcf = case_when(pcf == "propPro" ~ "Progression",
                         pcf == "propReg" ~ "Regression",
                         pcf == "propUnd" ~ "Undulation")) %>% 
  mutate(pcf = factor(pcf, levels = c("Progression", "Regression", "Undulation")))

fig_NHd <- ggplot(data = NHd, aes(x = time, y = pct, fill = pcf)) +
  geom_area(position = "stack") +
  labs(x = "Months following infection", y = "Proportion of attributable transmission", 
       fill = "Natural history pathways") +
  theme_bw() +
  theme(legend.key.size = unit(1.2, "lines"), legend.position = "bottom", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), panel.background = element_rect(fill = "white"), panel.border = element_blank(), 
        axis.line = element_line(colour = 'black', size = 0.5), axis.ticks = element_line(colour = 'black', size = 0.5))

png("figures/fig_NHd.png",width=600, height=400)
print(fig_NHd)
dev.off()

rm(list=ls()[! ls() %in% c("nreps", "nsteps", "RIs", "PCFsa", "f3", "RTM", "RTMd", "PCF", "PCFd", "NH", "NHd")])

# Output data for tables
table_050 <- import("./results/table_0.5ri_2m.Rdata", trust = TRUE)
write.csv(table_050,file="./table_050.csv")

table_025 <- import("./results/table_0.25ri_2m.Rdata", trust = TRUE)
write.csv(table_025,file="./table_025.csv")

table_100 <- import("./results/table_1ri_2m.Rdata", trust = TRUE)
write.csv(table_100,file="./table_100.csv")
