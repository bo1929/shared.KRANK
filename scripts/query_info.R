require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)

dists <- read_tsv("../data/dist_to_closest.txt", col_names = c("G1", "G2", "d"))
p1 <- ggplot(dists) + aes(y=d, x=reorder(G1,d)) +
  geom_point(size=0.25) +
  geom_hline(yintercept = c(0.001, 0.025, 0.05, 0.1, 0.2, 0.35), color="red", linetype=2) +
  theme_cowplot(font_size = 24) +
  labs(x="Query genome ID", y="Distance to the closest reference genome") +
  theme(
    axis.text.x = element_text(size=2, angle=90, vjust=0.5),
    axis.ticks.length.x = unit(.01, "cm")
  )
ggsave2("../figures/distance_to_closest.pdf", height = 8, width = 20)

ranks <- read_tsv("../data/query_ranks.tsv", col_names = c("genome", "kingdom",  "phylum", "class", "order", "family", "genus", "species"))
ranks$phylum <- as.factor(ranks$phylum)
ranks$class <- as.factor(ranks$class)
ranks <- ranks %>% group_by(phylum) %>% mutate(count_phylum = -n())

p2 <- ggplot(ranks) + geom_bar(aes(x=reorder(phylum, count_phylum))) +
  theme_cowplot(font_size = 24) +
  labs(x="Phylum ID (NCBI)", y="Number of query genomes") +
  theme(axis.text.x = element_text(angle=90, size=18, vjust = 0.5)) +
  scale_y_continuous(trans = "log1p", breaks = c(1, 10, 50, 100, 200))
  
plot_grid(p1, p2, nrow=2, rel_widths = c(4, 1), labels = c("A", "B"), label_size = 25)
ggsave2("../figures/query_details.pdf", height = 18, width = 20)