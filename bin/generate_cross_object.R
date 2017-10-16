library(qtl)
library(linkagemapping)
library(tidyverse)
library(snow)
data("N2xCB4856cross")
# Generate phenotype info
pheno <- readr::read_tsv("cross_obj_strains.tsv", col_names = "strain") %>%
  dplyr::mutate(pheno = 1) %>%
  dplyr::select(pheno, strain) %>% t()

pheno <- cbind(c("dummy_pheno", "strain"), pheno)
colnames(pheno) <- NULL

write.table(pheno, file = "cross_obj_pheno.tsv", row.names = F, col.names = F,  quote = F, sep = "\t", na = "-")

cross_obj <- read.cross(format = "csvsr",
                        genfile = "cross_obj_geno.tsv",
                        phefile = "cross_obj_pheno.tsv",
                        genotypes = c("N","C"),
                        alleles = c("N","C"),
                        sep = "\t")

# Fix sorting of strains

# Fix X Chromosome Name
names(cross_obj$geno)[names(cross_obj$geno) == "Xchr"] <- "X"
names(cross_obj$geno$X$map) <- gsub("Xchr_", "X_", names(cross_obj$geno$X$map))
colnames(cross_obj$geno$X$data) <- gsub("Xchr_","X_", colnames(cross_obj$geno$X$data))

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



# check for segregation distortion
gt <- geno.table(cross_obj)

# There is a skew on I, as predicted by peel/zeel

# compare genotypes between individuals
cg <- comparegeno(cross_obj)
png(file = "comparegeno.png")
hist(cg, breaks =100, xlab="Proportion of identical genotypes")
dev.off()

save(cg, file = "comparegeno.Rda")

png(file = "rug.png")
plot(0,type='n',axes=FALSE,ann=FALSE)
rug(cg)
dev.off()

identicals <- which(cg > .98, arr.ind = TRUE)

identicals_list <- apply(identicals, 2, function(x) {
  cross_obj$pheno$strain[x]
})

save(identicals_list, file = "identicals_list.Rda")

# There are 30 individuals that are 98% identical. These are all ECA strains. Which makes sense because they aren't advanced intercross lines.

# check marker order
pairrf <- est.rf(cross_obj)

save(pairrf, file = "estrf.Rda")

# no warnings, alleles likely not switched

checkAlleles(pairrf)
# no issues with alleles

# plot the recombination pattern
rfplot <- plot.rf(pairrf, alternate.chrid = TRUE)

# recombination looks good
CM_map <- est.map(pairrf, error.prob = 0.001, n.cluster = 6)

save(CM_map, file = "CM_map.Rda")

CM_cross_obj <- replace.map(pairrf, CM_map)

save(CM_cross_obj, file = "CM_cross_obj.Rda")
