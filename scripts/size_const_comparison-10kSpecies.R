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
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_shape_manual(values = c(15,17,19,8)) +
  theme_cowplot(font_size = 18) +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(
    axis.text.y = element_text(),
    axis.text.x = element_text(angle = 20)
  ) + theme(legend.position = "none")
p1
ggsave2("../figures/size_const_comparison-wrt_group_size.pdf", width = 7, height = 6)

options(scipen=10000)
ggplot(scores %>% mutate(
      Reference_genome_in_phylum = cut(
      Reference_genome_in_phylum,
      include.lowest = TRUE,
      dig.lab = 50,
      breaks = c(0, 10, 100, 1000))
    ) %>% filter(Method %in% c("None", "Mixed")) %>% filter(!is.na(Reference_genome_in_phylum))) +
  aes(y = F1, x = Reference_genome_in_phylum, color = Method, shape = Method) +
  # geom_point(position = position_dodge2(width = 0.5), size = 0.75, alpha = 0.25) +
  stat_summary(aes(y = F1, group = Method), geom = "line", color="darkgrey", linewidth=1) +
  stat_summary(size = 1) + 
  scale_colour_brewer(palette = "Dark2") +
  # scale_x_continuous(trans="log10") +
  theme_cowplot(font_size = 18) +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_text(), axis.text.x = element_text()) + labs(x="Number of genomes in the same phylum")
px <- ggplot(
  scores %>% filter(Reference_genome_in_phylum < 1000) %>% filter(Method %in% c("None", "Mixed")) %>% filter(!is.na(Taxonomic_rank)) %>% group_by(Method, Taxonomic_rank) %>%  summarise(Precision = mean(Precision), Recall = mean(Recall), F1 = mean(F1))
) +
  aes(x = Taxonomic_rank, y = F1, color = Method, shape = Method) +
  # geom_boxplot() +
  # geom_point(size = 5,alpha = 1) +
  geom_line(aes(group = Method, color = Method), linewidth = 2, alpha = 0.5) +
  geom_point(size = 5, alpha = 1)+
  # geom_smooth() + 
  # facet_wrap(c("expt")) +
  # facet_wrap(c("Taxonomic_rank")) +
  labs(shape = "Approach", colour = "Approach", x = "Taxonomic rank", y = "F1") +
  # geom_line(aes(group = Distance_to_closest), color = "grey20") +
  scale_colour_brewer(palette = "Dark2") +
  # scale_x_continuous(trans="log10") +
  theme_cowplot(font_size = 21) +
  theme(strip.background = element_rect(fill = "gray")) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
px

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
  labels = c("A", "B"),
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

ggplot(data=data.frame(x=0:1)) +
  geom_vline(xintercept = 0.7, color = "black", linetype = 2) +
  # geom_hline(yintercept = 0.1, color = "navy", linetype = 2) + geom_hline(yintercept = sqrt(0.1), color = "darkred", linetype = 2) + 
  geom_vline(xintercept = 0.1, color = "black", linetype = 4) +
  # geom_hline(yintercept = 0.25, color = "navy", linetype = 4) + geom_hline(yintercept = sqrt(0.25), color = "darkred", linetype = 4) + 
  stat_function(xlim = c(0, 1), fun = function(x) {x}, linewidth = 1, aes(color = "linear")) +
  stat_function(xlim = c(0, 1), fun = function(x) {sqrt(x)}, linewidth = 1, aes(color = "square root")) +
  theme_cowplot() +
  scale_color_manual(
    values=c("navy", "darkred"),
    breaks=c("linear", "square root")
  ) + 
  labs(y = "r(t)", x = "Total k-mer count ratio", color = "")