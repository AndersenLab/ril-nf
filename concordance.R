# Determine concordance of fastqs and among identicle strains.
library("tidyverse")

ind_calls <- dir(pattern="*.tsv")
fname <- readr::read_tsv(ind_calls[[1]])[[1]][[1]]

df <- readr::read_tsv(ind_calls[1], col_names = c("fname", "chrom_pos", "gt")) %>%
      dplyr::rename_(.dots = setNames("gt", fname)) %>%
      dplyr::select(-fname)


for (x in ind_calls[2:length(ind_calls)]) { 
  fname <- readr::read_tsv(x)[[1]][[1]]
  x <- readr::read_tsv(x,  col_names = c("fname", "chrom_pos", "gt")) %>%
       dplyr::rename_(.dots = setNames("gt", fname)) %>%
       dplyr::select(-fname)
  df <- dplyr::left_join(df, x, by = "chrom_pos")  
  df
}

df <- as.matrix(dplyr::select(df, -chrom_pos))

similarity <- apply(df,2,function(x)colSums(x!=df, na.rm = TRUE))
df_complete <- apply(df,2,function(x)colSums(is.na(x)==is.na(df), na.rm = TRUE))

diag(similarity) <- -1
diag(df_complete) <- -1

df <- tbl_df(similarity) %>% 
  dplyr::mutate(f1 = colnames(similarity) ) %>%
  tidyr::gather(f2, diff, -f1)

df_complete <- tbl_df(df_complete) %>% 
  dplyr::mutate(f1 = colnames(df_complete) ) %>%
  tidyr::gather(f2, sites, -f1)

df <- dplyr::left_join(df, df_complete, by = c("f1", "f2")) %>%
      dplyr::filter(sites != -1) %>%
      dplyr::mutate(concordance = 1 - (diff/sites)) %>%
      dplyr::mutate(SM1 = stringr::str_extract(f1, "([^\\.]+)")) %>%
      dplyr::mutate(SM2 = stringr::str_extract(f2, "([^\\.]+)")) %>%
      dplyr::mutate(same_SM = SM1 == SM2) %>%
      dplyr::select(SM1, SM2, same_SM, everything()) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(ff = as.character(list(sort(c(f1, f2))))) %>%
      dplyr::group_by(ff) %>%
      dplyr::distinct(.keep_all = TRUE) %>%
      dplyr::ungroup() %>%
      dplyr::select(-ff)

ggplot(df) +
  geom_histogram(aes(x=concordance, fill = same_SM), binwidth = 0.000025) +
  scale_fill_manual(values = c("#808080", "#0080FF")) +
  labs(x = "Concordance", y = "Number of Comparisons") +
  theme(axis.title = ggplot2::element_text(size=14, face="bold", color="black", vjust=5))

ggsave("fq_concordance.svg")
ggsave("fq_concordance.png")

# Save text files
readr::write_tsv(df, "fq_concordance.tsv")
