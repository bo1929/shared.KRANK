require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)
require(stringr)

scores <-
  read_csv("../results/cscores-10kSpecies_-_KRANK-rankingkmers_comparison.csv")
scores
scores$Taxonomic_rank <- factor(
  scores$Taxonomic_rank,
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
scores <- scores %>% filter(Distance_to_closest < 0.35)

scores$Method[scores$Method == "random k-mer ranking ~ k=32 w=35 h=12 b=16 l=2 mer-count const. (0.00)"] <-
  "random"
scores$Method[scores$Method == "negative species count ranking ~ k=32 w=35 h=12 b=16 l=2 mer-count const. (0.00)"] <-
  "species discrim."
scores$Method[scores$Method == "negative child taxon count ranking ~ k=32 w=35 h=12 b=16 l=2 mer-count const. (0.00)"] <-
  "children discrim."
scores$Method[scores$Method == "positive species count ranking ~ k=32 w=35 h=12 b=16 l=2 mer-count const. (0.00)"] <-
  "species common"
scores$Method[scores$Method == "positive child taxon count ranking ~ k=32 w=35 h=12 b=16 l=2 mer-count const. (0.00)"] <-
  "children common"
scores$Method[scores$Method == "weighted sum of species counts w.r.t. k-mer coverages ~ k=32 w=35 h=12 b=16 l=2 mer-count const. (0.00)"] <-
  "taxon covering"
scores_p1 <- scores %>%
  filter(Method %in% c("random", "species common", "species discrim.")) %>%
  # mutate(Distance_to_closest = cut(Distance_to_closest, include.lowest = TRUE, breaks = c(0, 0.001, 0.025, 0.05, 0.1, 0.2, 0.35))) %>%
  group_by(Method, Taxonomic_rank) %>%
  summarise(
    Precision = mean(Precision),
    Recall = mean(Recall),
    F1 = mean(F1)
  )
scores_p1$expt = "Species-based ranking (R)"
p1 <- ggplot(scores_p1) +
  aes(
    x = reorder(Method, F1),
    y = F1,
    color = Taxonomic_rank,
    shape = Method
  ) +
  geom_point(size = 5, alpha = 0.85) +
  facet_wrap(c("expt")) +
  labs(
    shape = "Ranking",
    colour = "Distance to the closest",
    x = "Ranking",
    y = "F1"
  ) +
  geom_line(aes(group = Taxonomic_rank), color = "grey20") +
  scale_colour_brewer(palette = "Paired") +
  theme_cowplot(font_size = 18) + scale_shape_manual(values = c(16, 15, 4)) +
  scale_x_discrete(
    labels = function(x)
      str_wrap(x, width = 10)
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
p1
scores_p2 = scores %>%
  filter(Method %in% c("random", "children common", "children discrim.")) %>%
  mutate(Distance_to_closest = cut(
    Distance_to_closest,
    include.lowest = TRUE,
    breaks = c(0, 0.001, 0.025, 0.05, 0.1, 0.2, 0.35)
  )) %>%
  group_by(Method, Taxonomic_rank) %>%
  summarise(
    Precision = mean(Precision),
    Recall = mean(Recall),
    F1 = mean(F1)
  )
scores_p2$expt = "Children-based ranking (R')"
p2 <- ggplot(scores_p2) +
  aes(
    x = reorder(Method, F1),
    y = F1,
    color = Taxonomic_rank,
    shape = Method
  ) +
  geom_point(size = 5, alpha = 0.85) +
  facet_wrap(c("expt")) +
  labs(
    shape = "Ranking",
    colour = "Distance to the closest",
    x = "Ranking",
    y = "F1"
  ) +
  geom_line(aes(group = Taxonomic_rank), color = "grey20") +
  scale_colour_brewer(palette = "Paired") +
  theme_cowplot(font_size = 18) + scale_shape_manual(values = c(8, 18, 16)) +
  scale_x_discrete(
    labels = function(x)
      str_wrap(x, width = 10)
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
p2
scores_p3 <- scores %>%
  filter(Method %in% c("random", "taxon covering")) %>%
  mutate(Distance_to_closest = cut(
    Distance_to_closest,
    include.lowest = TRUE,
    breaks = c(0, 0.001, 0.025, 0.05, 0.1, 0.2, 0.35)
  )) %>%
  group_by(Method, Taxonomic_rank) %>%
  summarise(
    Precision = mean(Precision),
    Recall = mean(Recall),
    F1 = mean(F1)
  )
scores_p3$expt <- "Weighted sum (R*)"
p3 <- ggplot(scores_p3) +
  aes(
    x = reorder(Method, F1),
    y = F1,
    color = Taxonomic_rank,
    shape = Method
  ) +
  geom_point(size = 5, alpha = 0.85) +
  facet_wrap(c("expt")) +
  labs(
    shape = "Ranking",
    colour = "Distance to the closest",
    x = "Ranking",
    y = "F1"
  ) +
  geom_line(aes(group = Taxonomic_rank), color = "grey20") +
  scale_colour_brewer(palette = "Paired") +
  theme_cowplot(font_size = 18) + scale_shape_manual(values = c(16, 17)) +
  scale_x_discrete(
    labels = function(x)
      str_wrap(x, width = 10)
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
p3

p4 <- ggplot(
  scores %>%
    filter(Method %in% c("random", "taxon covering")) %>%
    mutate(
      Distance_to_closest = cut(
        Distance_to_closest,
        include.lowest = TRUE,
        breaks = c(0, 0.001, 0.025, 0.05, 0.1, 0.2, 0.35)
      )
    ) %>%
    mutate(
      Reference_genome_in_phylum = cut(
        Reference_genome_in_phylum,
        include.lowest = TRUE,
        breaks = c(0, 5, 10, 50, 100, 500, 5000)
      )
    ) %>%
    group_by(Method, Distance_to_closest, Taxonomic_rank) %>%
    summarise(
      Precision = mean(Precision),
      Recall = mean(Recall),
      F1 = mean(F1)
    )
) +
  aes(
    x = Precision,
    y = Recall,
    color = Taxonomic_rank,
    shape = Method
  ) +
  facet_wrap("Distance_to_closest", scale = "free") +
  geom_point(size = 3.2, alpha = 0.85) +
  labs(
    shape = "Ranking",
    colour = "Distance to the closest",
    x = "Precision",
    y = "Recall"
  ) +
  geom_line(aes(group = Taxonomic_rank), color = "grey60") +
  scale_colour_brewer(palette = "Paired") +
  theme_cowplot(font_size = 18) + # scale_shape_manual(values = c(16, 17)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))

prow <- plot_grid(
  p1 + theme(legend.position = "none") + theme(axis.title.x = element_text(size =
                                                                             0)),
  p2 + theme(legend.position = "none") + theme(axis.title.x = element_text(size =
                                                                             0)),
  p3 + theme(legend.position = "none") + theme(axis.title.x = element_text(size =
                                                                             0)),
  ncol = 3,
  rel_heights = c(1, 1, 1)
)
p_cd <- plot_grid(
  prow,
  p4 + theme(legend.position = "none"),
  labels = c("C", "D"),
  label_size = 20
)
p_cd

scores <- read_csv("../results/cscores-10kSpecies_-_KRANK-sizeconst_comparison.csv")
scores
scores$Taxonomic_rank <- factor(
  scores$Taxonomic_rank,
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
scores <- scores %>% filter(Distance_to_closest < 0.35)
# scores$Taxon <- as.factor(scores$Taxon)
options(scipen = 10000)
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
    summarise(
      Recall = mean(Recall),
      Precision = mean(Precision),
      F1 = mean(F1)
    )
) +
  aes(
    x = reorder(Method, F1),
    y = F1,
    color = Taxonomic_rank,
    shape = Method
  ) +
  facet_wrap(. ~ Reference_genome_in_phylum) +
  geom_line(aes(group = Taxonomic_rank)) +
  geom_point(size = 2.5, alpha = 0.85) +
  labs(
    shape = "Constraint",
    colour = "Taxonomic rank",
    linetype = "Constraint",
    x = " ",
    y = "F1"
  ) +
  scale_colour_brewer(palette = "Paired") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  scale_shape_manual(values = c(15, 17, 19, 8)) +
  theme_cowplot(font_size = 18) +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_text(), axis.text.x = element_text(angle = 20)) + theme(legend.position = "none")
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
) + aes(
  color = Taxonomic_rank,
  x = Precision,
  y = Recall,
  shape = reorder(Method, Recall)
) +
  geom_point(size = 2.5, alpha = 0.85) +
  facet_wrap(vars(Distance_to_closest),
             scales = "free",
             nrow = 2) +
  labs(
    shape = "Constraint",
    colour = "Taxonomic rank",
    x = "Precision",
    y = "Recall"
  ) +
  geom_line(aes(group = Taxonomic_rank), color = "grey60") +
  scale_colour_manual(values = palette.colors(palette = "Paired", n = 7)[-1]) +
  theme_cowplot(font_size = 18) +
  theme(panel.spacing.x = unit(1.15, "lines"))
p2
p_ab <- plot_grid(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"),
  ncol = 2,
  rel_heights = c(1, 1.25),
  vjust = 2,
  hjust = 2,
  labels = c("A", "B"),
  label_size = 20,
  label_x = 0.075,
  label_y = 1
)
p_ab
plot_grid(p_ab, p_cd, nrow = 2)
ggsave2(
  "../figures/kmer_ranking_comparison-10kSpecies.pdf",
  width = 15,
  height = 11
)
