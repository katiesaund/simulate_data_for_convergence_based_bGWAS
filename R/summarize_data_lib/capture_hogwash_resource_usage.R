# To be run from data directory where hogwash results are saved
library(tidyverse)
args <- commandArgs(trailingOnly = TRUE)
num_tree <- as.numeric(args[1])
num_pheno <- as.numeric(args[2])
relative_path <-  args[3]
source(paste0(relative_path, "/simulate_data_for_convergence_based_bGWAS/R/summarize_data_lib/greatlakes_resource_usage.R"))

num_test <- 3 # phyc, sync, and continuous
num_signal <- 2 # BM and WN
num_row <- num_tree * num_pheno * num_test * num_signal
signals <- c("BM", "WN")

resource_df <- as.data.frame(matrix(NA, nrow = num_row, ncol = 6))
colnames(resource_df) <- c("Time (Hours)",
                           "Memory (GB)", 
                           "tree_ID", 
                           "pheno_ID", 
                           "test", 
                           "phylo_signal")
row_num <- 1
for (i in 1:num_tree) {
  for (j in 1:num_pheno) {
    for (k in 1:num_signal) {
      continuous_out_file_name <- paste0("continuous_pheno_", signals[k], "_tree_", i, "_pheno_", j, ".out")
      path_and_job_id <- readLines(con = continuous_out_file_name, n = 2)
      job_id <- path_and_job_id[2]
      time_and_mem <- GetResourceUsage(jobId = job_id) # Returns vector of length 2
      time <- time_and_mem[1]
      mem <- time_and_mem[2]
      resource_df[row_num, ] <- c(time, mem, i, j, "Continuous", signals[k])
      
      phyc_out_file_name <- paste0("phyc_discrete_pheno_", signals[k], "_tree_", i, "_pheno_", j, ".out")
      path_and_job_id <- readLines(con = phyc_out_file_name, n = 2)
      job_id <- path_and_job_id[2]
      time_and_mem <- GetResourceUsage(jobId = job_id) # Returns vector of length 2
      time <- time_and_mem[1]
      mem <- time_and_mem[2]
      resource_df[row_num + 1, ] <- c(time, mem, i, j, "PhyC", signals[k])
      
      sync_out_file_name <- paste0("synchronous_discrete_pheno_", signals[k], "_tree_", i, "_pheno_", j, ".out")
      path_and_job_id <- readLines(con = phyc_out_file_name, n = 2)
      job_id <- path_and_job_id[2]
      time_and_mem <- GetResourceUsage(jobId = job_id) # Returns vector of length 2
      time <- time_and_mem[1]
      mem <- time_and_mem[2]
      resource_df[row_num + 2, ] <- c(time, mem, i, j, "Synchronous", signals[k])
      row_num <- row_num + 3
    }
  }
}
write.csv(resource_df, 
          file = "hogwash_resource_usage.csv", 
          quote = FALSE, 
          row.names = FALSE)
