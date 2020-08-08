library(hogwash)
args <- commandArgs(trailingOnly = TRUE) # Grab arguments from the .sbat file
phenotype <- read.table(args[1],
                    sep = "\t",
                    row.names = 1,
                    header = TRUE,
                    stringsAsFactors = FALSE,
                    check.names = FALSE)
phenotype <- as.matrix(phenotype)

genotype <- read.table(args[2],
                    sep = "\t",
                    row.names = 1,
                    header = TRUE,
                    stringsAsFactors = FALSE,
                    check.names = FALSE)
genotype <- as.matrix(genotype)

tr  <- ape::read.tree(args[3])
f_name <- args[4]
out_dir <- args[5]
perm_num <- as.numeric(args[6])
fdr_value <- as.numeric(args[7])
bootstrap_threshold <- as.numeric(args[8])

test_type <- as.character(args[9])

group_method <- as.character(args[10])

key <- NULL
if (!is.na(args[11])) {
  key <- read.table(args[11],
                    sep = "\t",
                    row.names = 1,
                    header = TRUE,
                    stringsAsFactors = FALSE,
                    check.names = FALSE)
  key <- as.matrix(key)
}

hogwash(pheno = phenotype,
        geno = genotype,
        tree = tr,
        file_name = f_name,
        dir = out_dir,
        perm = perm_num,
        fdr = fdr_value,
        bootstrap = bootstrap_threshold,
        test = test_type,
        grouping_method = group_method, 
        group_genotype_key = key)
