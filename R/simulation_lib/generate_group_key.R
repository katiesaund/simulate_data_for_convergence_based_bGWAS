# Generate a group key mapping variants into genes. 

# Read in genotype used for Figure 4 in the paper
# Figure 4: phyc, BM, tree 1, pheno 1
genotype <- read.table("../data/simulated_genotype_for_discrete_pheno_BM_tree_1_pheno_1.tsv",
                       sep = "\t",
                       row.names = 1,
                       header = TRUE,
                       stringsAsFactors = FALSE,
                       check.names = FALSE)
genotype <- as.matrix(genotype)

# Initialize
group_names <- loci_names <- NULL
set.seed(1)

# Create a list of simulated loci names to use in the grouped key
# Start by including all of the genotypes and then repeating some of them at 
# least a few times (using the sample function). The repetition will simulate
# snps being put into multiple groups (overlapping genes, multiple pathways...)
loci_names <- c(colnames(genotype), 
                sample(colnames(genotype),
                       size = floor(0.2 * ncol(genotype)), 
                       replace = TRUE))
# Delete some of the loci as not all loci in my real data map a group
loci_to_drop_index <- sample(1:length(loci_names), 
                             size = floor(0.1 * length(loci_names)),
                             replace = FALSE)
loci_names <- loci_names[!1:length(loci_names) %in% loci_to_drop_index]

# Initialize the grouping key
grouping_key <- matrix(NA, nrow = length(loci_names), ncol = 2)
colnames(grouping_key) <- c("locus", "group")

# Create a list of names for the group -- it will be longer than needed. 
# Group names will be included in this list between 1 - 10 times (numbers I see
# in my real data). 
for (i in 1:length(loci_names)) {
  num_rep <- sample(1:10, 1)
  temp_name <- paste0("group_", rep(i, num_rep))
  group_names <- c(group_names, temp_name)
}

# Fill in grouping key & save
grouping_key[, 1] <- loci_names
grouping_key[, 2] <- group_names[1:length(loci_names)]
row.names(grouping_key) <- 

write.table(grouping_key, 
            sep = "\t", 
            col.names = TRUE, 
            quote = FALSE,
            file = "../data/group_key_for_phyc_BM_tree_1_pheno_1.tsv")
# hist(table(grouping_key[, 1])) # 1, 2, 3, or 4
# hist(table(grouping_key[, 2])) # Range of 1 - 10
