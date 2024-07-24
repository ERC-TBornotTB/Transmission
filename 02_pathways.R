# Clear environment
rm(list = ls())

# Load required libraries
library(tidyverse)
library(progress)

# Specify random number generation
set.seed(123)

# Set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

# Set cohort parameter values 
nreps   <- 1000        # Number of repetitions
nsteps  <- 121         # Number of time steps (months)
cutoffs <- c(0, 1, 2)  # Threshold for PCF cut-off (Baseline = 2m)  

# Relative infectiousness of subclinical disease
RIs <- c(0.0, 0.125, 0.25, 0.375, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75) 

if(!dir.exists("data/pathways")) dir.create(file.path("data/pathways"))

# Relative transmission ----
pb <- progress_bar$new(format = "[:bar] :percent", total = (nreps * length(RIs)))

for(q in seq_along(RIs)) {
  sub_transm <- RIs[q]
  sub_folder <- paste0("./data/pathways/RI_", sub_transm)
  dir.create(sub_folder, showWarnings = FALSE)  
  
  for(i in 1:nreps) {
    RTM <- as.data.frame(array(0, dim = c(nsteps, 0)))
    load(paste("./data/disease/disease_", i, ".Rdata", sep = ""))
    
    for(j in 1:ncol(disease)) {
      s_value <- NA
      
      for(k in 1:nrow(disease)) {
        if(disease[k, j] == "s") {
          if(is.na(s_value)) {
            s_value <- runif(1, min = 0.4 * sub_transm, max = 1.6 * sub_transm)
          }
          RTM[k, j] <- s_value
          } else if (disease[k, j] == "c") {
            RTM[k, j] <- 1
            s_value <- NA
          } else {
            RTM[k, j] <- 0
            s_value <- NA
          }
      }
    }
    
    save(RTM, file = paste(sub_folder, "/RTM_", i, ".Rdata", sep = ""))
    rm("disease", "RTM")
    pb$tick()
  }
}
cat(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "- Completed: Relative transmission\n"))

# Passive-case finding ----
pb <- progress_bar$new(format = "[:bar] :percent", total = nreps * (length(RIs) + length(cutoffs)) - nreps)

for(q in seq_along(RIs)) {
  sub_transm <- RIs[q]
  sub_folder <- paste0("./data/pathways/RI_", sub_transm)
  if(!dir.exists(sub_folder)) dir.create(sub_folder, recursive = TRUE)
  
  if(sub_transm == 0.5) {
    for(cutoff in cutoffs) {
      for(i in 1:nreps) {
        load(paste("./data/disease/disease_", i, ".Rdata", sep = ""))
        load(paste(sub_folder, "/RTM_", i, ".Rdata", sep = ""))
        screen <- as.data.frame(matrix(0, nrow = nsteps, ncol = ncol(disease)))
        noscreen <- as.data.frame(matrix(0, nrow = nsteps, ncol = ncol(disease)))
        consec_c <- 0
        is.screened <- FALSE 
        
        for(j in 1:ncol(disease)) {
          for(k in 1:nrow(disease)) {
            if(disease[k, j] == 'c') {
              consec_c <- consec_c + 1
              if(consec_c > cutoff) {
                is.screened <- TRUE
              }
            } else {
              consec_c <- 0
            }
          }
          
          if (is.screened == TRUE) {
            screen[, j] <- RTM[, j]
            noscreen[, j] <- 0
          } else {
            noscreen[, j] <- RTM[, j]
            screen[, j] <- 0 
          }
          is.screened <- FALSE 
        }
        
        cutoff_folder <- paste0(sub_folder, "/", cutoff, "m")
        if(!dir.exists(cutoff_folder)) dir.create(cutoff_folder, recursive = TRUE)
        save(screen, file = paste(cutoff_folder, "/screen_", i, ".Rdata", sep = ""))
        save(noscreen, file = paste(cutoff_folder, "/noscreen_", i, ".Rdata", sep = ""))
        
        rm(list = c("screen", "noscreen", "disease", "RTM", "consec_c"))
        gc()
        pb$tick()
      }
    }
  } else {
    for(i in 1:nreps) {
      load(paste("./data/disease/disease_", i, ".Rdata", sep = ""))
      load(paste(sub_folder, "/RTM_", i, ".Rdata", sep = ""))
      screen <- as.data.frame(matrix(0, nrow = nsteps, ncol = ncol(disease)))
      noscreen <- as.data.frame(matrix(0, nrow = nsteps, ncol = ncol(disease)))
      consec_c <- 0
      is.screened <- FALSE 
      
      for(j in 1:ncol(disease)) {
        for(k in 1:nrow(disease)) {
          if(disease[k, j] == 'c') {
            consec_c <- consec_c + 1
            if(consec_c > 2) {  
              is.screened <- TRUE
            }
          } else {
            consec_c <- 0
          }
        }
        
        if (is.screened == TRUE) {
          screen[, j] <- RTM[, j]
          noscreen[, j] <- 0
        } else {
          noscreen[, j] <- RTM[, j]
          screen[, j] <- 0 
        }
        is.screened <- FALSE 
      }
      
      save(screen, file = paste(sub_folder, "/screen_", i, ".Rdata", sep = ""))
      save(noscreen, file = paste(sub_folder, "/noscreen_", i, ".Rdata", sep = ""))
      
      rm(list = c("screen", "noscreen", "disease", "RTM", "consec_c"))
      gc()
      pb$tick()
    }
  }
}
cat(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "- Completed: Passive-case finding\n"))

# Pre-screening ----
pb <- progress_bar$new(format = "[:bar] :percent", total = nreps * (length(RIs) + length(cutoffs)) - nreps)

for(q in seq_along(RIs)) {
  sub_transm <- RIs[q]
  sub_folder <- paste0("./data/pathways/RI_", sub_transm)
  if(!dir.exists(sub_folder)) dir.create(sub_folder, recursive = TRUE)
  
  if(sub_transm == 0.5) {
    for(cutoff in cutoffs) {
      for(i in 1:nreps) {
        load(paste("./data/disease/disease_", i, ".Rdata", sep = ""))
        load(paste(paste0(sub_folder, "/", cutoff, "m"), "/screen_", i, ".Rdata", sep = ""))
        
        prescreen <- as.data.frame(matrix(0, nrow = nrow(screen), ncol = ncol(screen)))
        
        for(j in 1:ncol(disease)) {
          consec_c <- 0
          col_length <- 0
          post_screen <- 0
          screen_point <- 0
          is.screened <- FALSE
          screen_point_set <- FALSE
          
          for(k in 1:nrow(disease)) {
            if(disease[k , j] == 'c') {
              consec_c <- consec_c + 1
              col_length <- col_length + 1
              
              if (consec_c > cutoff && !screen_point_set) {
                is.screened <- TRUE
                screen_point <- col_length - 1
                screen_point_set <- TRUE
              }
              
            } else {
              consec_c <- 0
              col_length <- col_length + 1
            }
          }
          
          if(is.screened) {
            post_screen <- post_screen + 1
            prescreen[1:screen_point, j] <- screen[1:screen_point, j]
            prescreen <- prescreen[1:nsteps, ]
          } else {
            prescreen[1:nsteps, j] <- 0
          }
        }
        
        cutoff_folder <- paste0(sub_folder, "/", cutoff, "m")
        if(!dir.exists(cutoff_folder)) dir.create(cutoff_folder, recursive = TRUE)
        save(prescreen, file = paste0(cutoff_folder, "/prescreen_", i, ".Rdata"))
        
        rm(list = c("prescreen", "disease", "consec_c", "col_length", "post_screen", "is.screened", "screen_point", "screen_point_set", "screen"))
        gc()
        pb$tick()
      }
    }
  } else {
    for(i in 1:nreps) {
      load(paste("./data/disease/disease_", i, ".Rdata", sep = ""))
      load(paste(sub_folder, "/screen_", i, ".Rdata", sep = ""))
      
      prescreen <- as.data.frame(matrix(0, nrow = nrow(screen), ncol = ncol(screen)))
      
      for(j in 1:ncol(disease)) {
        consec_c <- 0
        col_length <- 0
        post_screen <- 0
        screen_point <- 0
        is.screened <- FALSE
        screen_point_set <- FALSE
        
        for(k in 1:nrow(disease)) {
          if(disease[k , j] == 'c') {
            consec_c <- consec_c + 1
            col_length <- col_length + 1
            
            if (consec_c > 2 && !screen_point_set) {
              is.screened <- TRUE
              screen_point <- col_length - 1
              screen_point_set <- TRUE
            }
            
          } else {
            consec_c <- 0
            col_length <- col_length + 1
          }
        }
        
        if(is.screened) {
          post_screen <- post_screen + 1
          prescreen[1:screen_point, j] <- screen[1:screen_point, j]
          prescreen <- prescreen[1:nsteps, ]
        } else {
          prescreen[1:nsteps, j] <- 0
        }
      }
      
      save(prescreen, file = paste0(sub_folder, "/prescreen_", i, ".Rdata"))
      
      rm(list = c("prescreen", "disease", "consec_c", "col_length", "post_screen", "is.screened", "screen_point", "screen_point_set", "screen"))
      gc()
      pb$tick()
    }
  }
}
cat(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "- Completed: Pre-screening\n"))

# Natural history ----
pb <- progress_bar$new(format = "[:bar] :percent", total = (nreps * length(RIs)))

for(q in seq_along(RIs)) {
  sub_transm <- RIs[q]
  sub_folder <- paste0("./data/pathways/RI_", sub_transm)
  if(!dir.exists(sub_folder)) dir.create(sub_folder, recursive = TRUE)
  
  for(i in 1:nreps) {
    load(paste("./data/disease/disease_", i, ".Rdata", sep = ""))
    load(paste(sub_folder, "/RTM_", i, ".Rdata", sep = ""))
    progress <- as.data.frame(matrix(0, nrow = nsteps, ncol = ncol(disease)))
    regress <- as.data.frame(matrix(0, nrow = nsteps, ncol = ncol(disease)))
    undulate <- as.data.frame(matrix(0, nrow = nsteps, ncol = ncol(disease)))
    
    for(j in 2:ncol(disease)) {
      P <- TRUE
      R <- FALSE
      U <- FALSE
      
      for(k in 2:nrow(disease)) {
        if((disease[k, j] == "s" & disease[k - 1, j] == "c") |
           (disease[k, j] == "m" & disease[k - 1, j] == "s") |
           (disease[k, j] == "r" & disease[k - 1, j] == "m")) {
          R <- TRUE
          }
        
        if(R == TRUE & 
           ((disease[k, j] == "c" & disease[k - 1, j] == "s") |
            (disease[k, j] == "s" & disease[k - 1, j] == "m") |
            (disease[k, j] == "m" & disease[k - 1, j] == "r"))) {
          U <- TRUE
        }
        
        if(R == FALSE | U == FALSE) {
          P <- TRUE
        }
      }
      
      if(P == TRUE & R == FALSE) {
        progress[, j] <- RTM[, j]
      } 
      
      if(R == TRUE & U == FALSE) {
        regress[, j] <- RTM[, j]
      }  
      
      if(U == TRUE) {
        undulate[, j] <- RTM[, j]
      } 
    }
    
    save(progress, file = paste(sub_folder,"/progress_", i, ".Rdata", sep = ""))
    save(regress, file = paste(sub_folder,"/regress_", i, ".Rdata", sep = ""))
    save(undulate, file = paste(sub_folder,"/undulate_", i, ".Rdata", sep = ""))
    rm(list = c("progress", "regress", "undulate", "disease", "RTM"))
    gc()
    pb$tick()
  }
}
cat(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "- Completed: Natural history\n"))

# Rapid progression ----
pb <- progress_bar$new(format = "[:bar] :percent", total = (nreps * length(RIs)))

for(q in seq_along(RIs)) {
  sub_transm <- RIs[q]
  sub_folder <- paste0("./data/pathways/RI_", sub_transm)
  if(!dir.exists(sub_folder)) dir.create(sub_folder, recursive = TRUE)
  
  for(i in 1:nreps) {
    load(paste("./data/disease/disease_", i, ".Rdata", sep = ""))
    load(paste(sub_folder, "/RTM_", i, ".Rdata", sep = ""))
    rapid <- as.data.frame(matrix(0, nrow = nsteps, ncol = ncol(disease)))
    
    for(j in 1:ncol(disease)) {
      if(nrow(disease) == 0) {
        rapid[, j] <- 0
      }
      else if(any(disease[1:25, j] %in% c('s', 'c'))) {
        rapid[, j] <- RTM[, j]
      }
    }
    
    save(rapid, file = paste0(sub_folder, "/rapid_", i, ".Rdata"))
    rm(list = c("rapid", "disease", "RTM" ))
    gc()
    pb$tick()
  }
}
cat(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "- Completed: Rapid progression\n"))

# Remote progression ----
pb <- progress_bar$new(format = "[:bar] :percent", total = (nreps * length(RIs)))

for(q in seq_along(RIs)) {
  sub_transm <- RIs[q]
  sub_folder <- paste0("./data/pathways/RI_", sub_transm)
  if(!dir.exists(sub_folder)) dir.create(sub_folder, recursive = TRUE)
  
  for(i in 1:nreps) {
    load(paste("./data/disease/disease_", i, ".Rdata", sep = ""))
    load(paste(sub_folder,"/RTM_", i, ".Rdata", sep = ""))
    remote <- as.data.frame(matrix(0, nrow = nsteps, ncol = ncol(disease)))
    
    for(j in 1:ncol(disease)) { 
      if(nrow(disease) == 0) {
        remote[, j] <- 0
      } 
      else if(!any(disease[1:25, j] %in% c('s', 'c'))) {
        remote[, j] <- RTM[, j]
      } else {
        remote[, j] <- 0
      }
    }
    
    save(remote, file = paste0(sub_folder, "/remote_", i, ".Rdata"))
    rm(list = c("remote", "disease", "RTM"))
    gc()
    pb$tick()
  }
}
cat(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "- Completed: Remote progression\n"))
