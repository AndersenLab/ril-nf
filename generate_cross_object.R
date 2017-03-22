


# Generate phenotype info
pheno <- readr::read_tsv("output_strains.tsv", col_names = "strain") %>%
dplyr::mutate(pheno = 1) %>% 
dplyr::select(pheno, strain) %>% t() 
pheno <- cbind(c("dummy_pheno", "strain"), pheno)
colnames(pheno) <- NULL
write.table(pheno, file = "cross_obj_pheno.tsv", row.names = F, col.names = F,  quote = F, sep = "\t", na = "-")

# Set pseudo CM distance
crobj <- readr::read_tsv("/Users/dancook/Dropbox/Andersenlab/LabFolders/Dan/Andersen-Lab-RILs/_data/cross_obj_geno.tsv")
crobj[,3] <- crobj[,3]/100000000
names(crobj)[is.na(names(crobj))] <- ""
write.table(crobj, file = "cross_obj_geno_f.tsv", row.names = F, col.names = T,  quote = F, sep = "\t", na = "-")


cross_obj <- read.cross(format = "csvsr",
             genfile = "cross_obj_geno_f.tsv",
             phefile = "cross_obj_pheno.tsv",
             genotypes = c("N","H","C"),
             alleles = c("N2","CB4856"),
             sep = "\t")
