require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)

counts_lca <- read_tsv("../results/comparison-LCAcounts-10kSpecies.txt")

counts_lca <- counts_lca %>% filter(Rank %in% c("superkingdom", "phylum", "class", "order", "family", "genus", "species"))
counts_lca$Rank <- factor(
  counts_lca$Rank,
  levels = c("phylum", "class", "order", "family", "genus", "species")
)

ggplot(
  counts_lca,
  aes(y = Count, color = Method, shape = Method, linetype = Method)
) +
  facet_wrap(vars(Rank), scale = "free") +
  stat_ecdf()+
  scale_colour_brewer(palette = "Dark2") +
  theme_cowplot(font_size = 17) +
  scale_y_continuous(trans = "log1p", breaks = c(1, 10, 100, 1000, 10000, 100000))

ggplot(
  counts_lca,
  aes(x = Number_of_genomes,  y = Count, color = Method, shape = Method, linetype = Method)
) +
  facet_wrap(vars(Rank), scale = "free") +
  geom_point(alpha = 0.2) +
  stat_smooth(se = F, method = "glm", size = 1) +
  scale_colour_brewer(palette = "Dark2") +
  theme_cowplot(font_size = 17) +
  scale_x_continuous(trans = "log1p", breaks = c(1, 10, 100, 1000, 10000)) +
  scale_y_continuous(trans = "log1p", breaks = c(1, 10, 100, 1000, 10000, 100000))

ggplot(
  counts_lca %>% mutate(Number_of_genomes = cut( Number_of_genomes, include.lowest = TRUE, breaks = c(0, 1, 10, 100, 1000, 10000, 100000))),
  aes(x = Number_of_genomes, y = Count, color = Method, shape = Method, linetype = Method)) +
facet_wrap(vars(Rank))+
geom_boxplot() +
geom_point(position = "jitter", alpha = 0.1) +
scale_colour_brewer(palette = "Dark2") +
theme_cowplot(font_size = 17) +
scale_y_continuous(trans = "log1p", breaks = c(1, 10, 100, 1000, 10000, 100000)) +
theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1))

ggplot(
  counts_lca,
  aes(x = Rank, y = Count, color = Method, shape = Method)
) +
  geom_boxplot() +
  scale_colour_brewer(palette = "Dark2") +
  theme_cowplot(font_size = 17) +
  scale_y_continuous(trans = "log1p", breaks = c(1, 10, 100, 1000, 10000, 100000))