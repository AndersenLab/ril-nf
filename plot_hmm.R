library(tidyverse)

df <- readr::read_tsv("gt_hmm_fill.tsv") %>%
  dplyr::group_by(sample) %>%
  dplyr::filter(chrom != "MtDNA") %>%
  dplyr::mutate(gt = ifelse(gt == 1, "CB4856", "N2"))

df$index <- dplyr::group_indices(df)

strain_index <- df$sample
names(strain_index) <- df$index + 0.5

ggplot(df) +
  geom_rect(aes(xmin = start, xmax = end, ymin = index, ymax = index + 1, fill = gt)) +
  scale_fill_manual(values = c("#0080FF", "#FF8000")) +
  facet_grid(.~chrom, scales="free", space="free") +
  theme_bw() +
  scale_x_continuous(labels = function(x) { x/1e6 }, expand = c(0,0)) +
  scale_y_continuous(breaks = unique(df$index) + 0.5, labels = function(x) { strain_index[as.character(x)] }, expand = c(0,0)) + 
  theme(strip.background = element_blank(),
        legend.position = "None")

ggsave("gt_hmm.svg")
ggsave("gt_hmm.png")