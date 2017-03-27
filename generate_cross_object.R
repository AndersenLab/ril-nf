library(qtl)
library(linkagemapping)
library(tidyverse)
data("N2xCB4856cross")
# Generate phenotype info
pheno <- readr::read_tsv("output_strains.tsv", col_names = "strain") %>%
  dplyr::mutate(pheno = 1) %>% 
  dplyr::select(pheno, strain) %>% t() 
pheno <- cbind(c("dummy_pheno", "strain"), pheno)
colnames(pheno) <- NULL
write.table(pheno, file = "cross_obj_pheno.tsv", row.names = F, col.names = F,  quote = F, sep = "\t", na = "-")

cross_obj <- read.cross(format = "csvsr",
                        genfile = "cross_obj_geno.tsv",
                        phefile = "cross_obj_pheno.tsv",
                        genotypes = c("N","H","C"),
                        alleles = c("N","C"),
                        sep = "\t")

# Fix X Chromosome Name
names(cross_obj$geno)[names(cross_obj$geno) == "Xchr"] <- "X"

# Remove dummy pheno
cross_obj$pheno$dummy_pheno <- NULL

# Get rid of funky row names
row.names(cross_obj$pheno) <- NULL

# Merge in set data
cross_obj$pheno <- dplyr::left_join(cross_obj$pheno, N2xCB4856cross$pheno) %>%
  dplyr::mutate(set = ifelse(is.na(set), 3, set))

# Replace 3 with 2
for(x in names(cross_obj$geno)) { 
  cross_obj$geno[[x]]$data[cross_obj$geno[[x]]$data == 3] <- 2
}

save(cross_obj, file = "cross_obj.Rdata" )
