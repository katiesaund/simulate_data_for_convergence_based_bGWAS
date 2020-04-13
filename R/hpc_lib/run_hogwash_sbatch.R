# devtools::install_github("katiesaund/hogwash", ref = "master", dependencies = FALSE)
library(hogwash)
args <- commandArgs(trailingOnly = TRUE) # Grab arguments from the PBS script
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

key <- NULL
if (!is.na(args[10])) {
  key <- read.table(args[10],
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
        group_genotype_key = key)