library(tidyverse)
library(cowplot)

fq_coverage <- readr::read_tsv("fq_coverage.tsv")
SM_coverage <- readr::read_tsv("SM_coverage.tsv")

pub_theme <- ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(size=14, face="bold", color="black"),
                 axis.text.y = ggplot2::element_text(size=14, face="bold", color="black"),
                 axis.title.x = ggplot2::element_text(size=14, face="bold", color="black", vjust=-.3),
                 axis.title.y = ggplot2::element_text(size=14, face="bold", color="black", vjust=2),
                 strip.text.x = ggplot2::element_text(size=16, face="bold", color="black"),
                 strip.text.y = ggplot2::element_text(size=16, face="bold", color="black"),
                 axis.ticks= element_line( color = "black", size = 0.25), 
                 legend.position="none",
                 plot.margin = unit(c(1.0,0.5,0.5,1.0),"cm"),
                 plot.title = ggplot2::element_text(size=24, face="bold", vjust = 1),
                 panel.background = ggplot2::element_rect(color = "black",size=0.75),
                 strip.background = ggplot2::element_rect(color = "black", size = 0.75)) 


a1 <- ggplot(fq_coverage) +
  geom_histogram(aes(x = coverage)) +
  geom_vline(aes(xintercept = mean(coverage)), color = "red") +
  labs(x = "Coverage", y = "Count", title = "Fastq Coverage") +
  scale_x_continuous(limits = c(0, 6)) +
  pub_theme

a2 <- ggplot(SM_coverage) +
  geom_histogram(aes(x = coverage)) +
  geom_vline(aes(xintercept = mean(coverage)), color = "red") +
  labs(x = "Coverage", y = "Count", title = "SM Coverage") +
  scale_x_continuous(limits = c(0, 6)) +
  pub_theme


cowplot::plot_grid(a1, a2, nrow = 2)


ggsave("coverage_comparison.png", width = 8, height = 10)
ggsave("coverage_comparison.svg", width = 8, height = 10)


un <- readr::read_tsv("SM_bam_idxstats.tsv") %>%
  dplyr::group_by(SM) %>%
  dplyr::mutate(unmapped = sum(unmapped_reads) / (sum(mapped_reads) + sum(unmapped_reads)))

ggplot(un) +
  geom_histogram(aes(x = unmapped)) +
  geom_vline(aes(xintercept = mean(unmapped)), color = "red") +
  labs(x = "Coverage", y = "Count", title = "Unmapped Read %") +
  pub_theme

ggsave("unmapped_reads.png")
ggsave("unmapped_reads.svg")


dup <- readr::read_tsv("bam_duplicates.tsv")

ggplot(dup) +
  geom_histogram(aes(x = percent_duplication), binwidth = 0.01) +
  geom_vline(aes(xintercept = mean(percent_duplication)), color = "red") +
  labs(x = "Coverage", y = "Count", title = "% Duplicated") +
  scale_x_continuous(limits = c(0,1)) +
  pub_theme

ggsave("duplicates.png")
ggsave("duplicates.svg")

