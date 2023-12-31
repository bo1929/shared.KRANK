require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)

scores <- read_csv("../results/cscores-10kSpecies_-_KRANK-sizeconst_comparison.csv")
scores
scores$Taxonomic_rank <- factor(
  scores$Taxonomic_rank,
  levels = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
)
scores <- scores %>% filter(Distance_to_closest < 0.35)
# scores$Taxon <- as.factor(scores$Taxon)
options(scipen=10000)
scores$Method[scores$Method == "k=32 w=35 h=12 b=16 l=2 no const. (0.00)"] <- "None"
scores$Method[scores$Method == "k=32 w=35 h=12 b=16 l=2 species-count const. (0.00)"] <- "species"
scores$Method[scores$Method == "k=32 w=35 h=12 b=16 l=2 mixed const. (mer-count and no) (0.00)"] <- "Mixed"
scores$Method[scores$Method == "k=32 w=35 h=12 b=16 l=2 mer-count const. (0.00)"] <- "k-mer"

p1 <- ggplot(
  scores %>% mutate(
    Reference_genome_in_phylum = cut(
      Reference_genome_in_phylum,
      include.lowest = TRUE,
      breaks = c(0, 50, 500, 5000)
    )
  ) %>% 
    group_by(Taxonomic_rank, Reference_genome_in_phylum, Method) %>%
    summarise(Recall = mean(Recall), Precision = mean(Precision), F1 = mean(F1))
) +
  aes(x=reorder(Method,F1), y = F1, color = Taxonomic_rank, shape = Method) +
  facet_wrap(.~Reference_genome_in_phylum) +
  geom_line(aes(group=Taxonomic_rank))+
  geom_point(size = 2.5, alpha = 0.85) +
  labs(shape = "Constraint", colour = "Taxonomic rank", linetype = "Constraint", x = "Constraint", y = "F1") +
  scale_colour_brewer(palette = "Paired") +
  scale_shape_manual(values = c(15,17,19,8)) +
  theme_cowplot(font_size = 18) +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(
    axis.text.y = element_text(),
    axis.text.x = element_text(angle = 20)
  ) + guides(color=guide_legend(nrow=1, byrow=TRUE))
p1

p2 <- ggplot(
  scores %>% 
    filter(Method %in% c("Mixed", "None")) %>%
    filter(Taxonomic_rank != "kingdom") %>%
    mutate(
      Distance_to_closest = cut(
        Distance_to_closest,
        include.lowest = TRUE,
        breaks = c(0, 0.001, 0.025, 0.05, 0.1, 0.2, 0.35)
      )
    ) %>%
    group_by(Taxonomic_rank, Method, Distance_to_closest) %>%
    summarise(Precision = mean(Precision), Recall = mean(Recall))
) + aes(color = Taxonomic_rank, x=Precision, y=Recall, shape=reorder(Method, Recall)) + 
  geom_point(size = 2.5, alpha = 0.85) +
  facet_wrap(vars(Distance_to_closest), scales = "free", nrow = 2) +
  labs(shape = "Constraint", colour = "Taxonomic rank", x = "Precision", y = "Recall") +
  geom_line(aes(group=Taxonomic_rank),color="grey60")+
  scale_colour_manual(values = palette.colors(palette = "Paired", n = 7)[-1]) +
  theme_cowplot(font_size = 18) +
  theme(panel.spacing.x = unit(1.15, "lines"))
p2

prow <- plot_grid(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"),
  ncol = 2,
  rel_heights = c(1, 1.25),
  vjust = 2,
  hjust = 2,
  labels = c("B", "C"),
  label_size = 20,
  label_x = 0.075, label_y = 1
)
prow
legend <- get_legend(
  p1 + theme(legend.box.margin = margin(6, 0, 0, 0)) + 
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.direction = "horizontal",
      legend.box = "vertical"
    ) + guides(color=guide_legend(nrow=1, byrow=TRUE))
)

taxa_counts <- read_tsv("../data/ref_taxa_counts.txt")
taxa_counts$Rank <- factor(
  taxa_counts$Rank,
  levels = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
)

p3 <- ggplot(taxa_counts %>% filter(Rank != "kingdom"),
             aes(x=Rank, y=Count)) +
  geom_violin(aes(color=Rank), linewidth=0.8, scale = "width", color="gray20", trim = TRUE, draw_quantiles = TRUE) +
  geom_point(aes(x=Rank, y=Count, color=Rank), alpha = 0.75, size = 1, position = "jitter") +
  stat_summary(color="black", fun = median) +
  scale_y_continuous(trans = "log2") +
  theme_cowplot(font_size = 18) +
  theme(strip.background = element_rect(fill = "gray")) +
  scale_colour_manual(values = palette.colors(palette = "Paired", n = 7)[-1]) +
  theme(
    panel.spacing.x = unit(1, "lines"),
    axis.title.x = element_blank() 
  ) +
  labs(y = "Number of\ngenomes in taxon", x = "")
p3
ggsave2("../figures/num_genomes_per_taxon.pdf", width = 12, height = 4)

plot_grid(p3 + theme(legend.position = "none"), plot_grid(prow, legend, rel_heights =  c(3, .4), nrow=2), ncol=1, rel_widths = c(1, 5), rel_heights = c(1, 3), labels = c("A", ""), label_size = 20)
ggsave2("../figures/size_const_comparison-10kSpecies.pdf", width = 14, height = 9.5)
