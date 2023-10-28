require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)
taxa_counts <- read_tsv("../data/ref_taxa_counts.txt")
taxa_counts$Rank <- factor(
  taxa_counts$Rank,
  levels = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
)

ggplot(taxa_counts %>% filter(Rank != "kingdom"),
             aes(x=Rank, y=Count, color=Rank)) +
  geom_violin(scale = "width", color="gray40", trim = TRUE) +
  geom_point(aes(x=Rank, y=Count), alpha = 0.9, size = 1, position = "jitter") +
  scale_y_continuous(trans = "log2") +
  theme_cowplot(font_size = 20) +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_blank(),
  ) +
  scale_colour_brewer(palette = "Paired") +
  theme(panel.spacing.x = unit(1, "lines")) +
  labs(y = "# of genomes in the taxon", x = "") +
  theme(
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16)
  )

