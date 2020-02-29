# devtools::install_github("katiesaund/hogwash", ref = "master", dependencies = FALSE)
library(hogwash)
args <- commandArgs(trailingOnly = TRUE) # Grab arguments from the PBS script
pheno <- read.table(args[1],
                    sep = "\t",
                    row.names = 1,
                    header = TRUE,
                    stringsAsFactors = FALSE,
                    check.names = FALSE)
pheno <- as.matrix(pheno)

geno <- read.table(args[2],
                    sep = "\t",
                    row.names = 1,
                    header = TRUE,
                    stringsAsFactors = FALSE,
                    check.names = FALSE)
geno <- as.matrix(geno)

tree  <- ape::read.tree(args[3])
file_name <- args[4]
dir <- args[5]
perm <- as.numeric(args[6])
fdr <- as.numeric(args[7])
bootstrap <- as.numeric(args[8])

group_genotype_key <- NULL
if (!is.na(args[9])) {
  group_genotype_key <- read.table(args[9],
                     sep = "\t",
                     row.names = 1,
                     header = TRUE,
                     stringsAsFactors = FALSE,
                     check.names = FALSE)
  group_genotype_key <- as.matrix(group_genotype_key)
}

hogwash(pheno,
        geno,
        tree,
        file_name,
        dir,
        perm,
        fdr,
        bootstrap,
        group_genotype_key)
