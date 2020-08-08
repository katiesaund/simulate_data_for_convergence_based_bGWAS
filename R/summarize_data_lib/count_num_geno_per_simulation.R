# Count the number of genotypes input into each hogwash run, or: how many
# genotypes did we actually simulate? We do this by reading in each simulated
# genotype matrix and counting the number of entries. 
# Record summary statistics in the output file. 
data_type <- c("discrete", "continuous")
signal <- c("BM", "WN")
args <- commandArgs(trailingOnly = TRUE)
num_tree <- args[1]
trees <- c(1:num_tree)

geno_count_mat <- matrix(0, nrow = 16, ncol = 3)
colnames(geno_count_mat) <- c("data_type", "signal", "number_of_genotypes")
row_count <- 0
for (i in 1:length(data_type)) {
  for (j in 1:length(signal)) {
    for (k in 1:length(trees)) {
        fname <- paste0("simulated_genotype_for_", 
                        data_type[i], 
                        "_pheno_", 
                        signal[j], 
                        "_tree_", 
                        trees[k], 
                        "_pheno_1.tsv")
        geno_mat <- read.table(fname,
                               header = TRUE, 
                               sep = "\t",
                               stringsAsFactors = FALSE ,
                               row.names = 1)
        num_geno <- ncol(geno_mat)
        row_count <- row_count + 1
        geno_count_mat[row_count, ] <- c(data_type[i], signal[j], num_geno)
    }
  }
}

summary(as.numeric(geno_count_mat[geno_count_mat[, 1] == "discrete", 3]))
summary(as.numeric(geno_count_mat[geno_count_mat[, 1] == "continuous", 3]))

summary_mat <- rbind(summary(as.numeric(geno_count_mat[geno_count_mat[, 1] == "discrete", 3])),
                     summary(as.numeric(geno_count_mat[geno_count_mat[, 1] == "continuous", 3])))
row.names(summary_mat) <- c("discrete", "continuous")
write.csv(summary_mat, "number_of_genotypes_simulated.csv")
