require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)

taxa_counts <- read_tsv("../data/ref_taxa_counts.txt")
taxa_counts$Rank <- factor(
  taxa_counts$Rank,
  levels = c(
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
  )
)

ggplot(taxa_counts %>% filter(Rank == "class"), aes(x = Count, fill = Rank)) +
  # geom_violin(scale = "width", color="gray40", trim = TRUE) +
  # geom_point(aes(x=Rank, y=Count), alpha = 0.9, size = 5, position = "jitter") +
  stat_bin(color = "black") +
  # stat_density() +
  scale_y_continuous(trans = "log1p") +
  theme_cowplot(font_size = 20) +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_blank(), ) +
  scale_fill_brewer(palette = "Dark2") +
  theme(panel.spacing.x = unit(1, "lines")) +
  labs(x = "# of genomes in the same class", y = "Number of taxon") +
  theme(axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16))

species_counts <- read_tsv("../data/ref_taxa_counts.txt")
ggplot(species_counts, aes(x = Count)) +
  # geom_violin(scale = "width", color="gray40", trim = TRUE) +
  # geom_point(aes(x=Rank, y=Count), alpha = 0.9, size = 5, position = "jitter") +
  stat_bin(fill = "navy", color = "white") +
  # stat_density() +
  scale_y_continuous(trans = "log1p", breaks = c(1, 10, 100, 1000, 10000)) +
  scale_x_continuous(trans = "log1p", breaks = c(1, 10, 100, 1000, 10000)) +
  theme_cowplot(font_size = 27) +
  theme(strip.background = element_rect(fill = "gray")) +
  labs(x = "# of genomes sequenced", y = "# of species")

ggplot(taxa_counts %>% filter(Rank != "kingdom"), aes(x = Rank, y = Count)) +
  geom_boxplot(
    color = "gray40",
    outlier.shape = NA,
    linewidth = 0.75
  ) +
  geom_point(
    aes(x = Rank, y = Count, color = Rank),
    alpha = 0.7,
    size = 0.5,
    position = "jitter"
  ) +
  scale_y_continuous(trans = "log2") +
  theme_cowplot(font_size = 18) +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_blank()) +
  scale_colour_manual(values = palette.colors(n = 7, palette = "Paired")[2:7]) +
  theme(panel.spacing.x = unit(1, "lines")) +
  theme(legend.position = "none") +
  labs(y = "Number of genomes in the taxon", x = "") +
  theme(axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16))