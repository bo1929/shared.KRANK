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
    "superkingdom",
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
  "negative species"
scores$Method[scores$Method == "negative child taxon count ranking ~ k=32 w=35 h=12 b=16 l=2 mer-count const. (0.00)"] <-
  "negative children"
scores$Method[scores$Method == "positive species count ranking ~ k=32 w=35 h=12 b=16 l=2 mer-count const. (0.00)"] <-
  "positive species"
scores$Method[scores$Method == "positive child taxon count ranking ~ k=32 w=35 h=12 b=16 l=2 mer-count const. (0.00)"] <-
  "positive children"
scores$Method[scores$Method == "weighted sum of species counts w.r.t. k-mer coverages ~ k=32 w=35 h=12 b=16 l=2 mer-count const. (0.00)"] <-
  "weighted sum"

p1 <- ggplot(
  scores %>%
    filter(Method %in% c("random", "positive species", "negative species")) %>%
    mutate(Distance_to_closest = cut(Distance_to_closest, include.lowest = TRUE, breaks = c(0, 0.001, 0.025, 0.05, 0.1, 0.2, 0.35))) %>%
    group_by(Method, Distance_to_closest) %>%
    summarise(Precision = mean(Precision), Recall = mean(Recall), F1 = mean(F1))
  ) +
  aes(x = Method, y = F1, color = Distance_to_closest, shape = Method) +
  geom_point(size = 5, alpha = 0.85) +
  labs(title = "Species-based ranking", shape = "Ranking", colour = "Distance to the closest", x = "Ranking", y = "F1") +
  geom_line(aes(group = Distance_to_closest), color = "grey20") +
  scale_colour_brewer(palette = "Paired") + 
  theme_cowplot(font_size = 16) + scale_shape_manual(values = c(4, 15, 16)) +
  theme(
    plot.title = element_text(size =  16),
    axis.text.y = element_text(),
    axis.text.x = element_text(angle = 20, vjust = 0.75),
  ) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
p1

p2 <- ggplot(
  scores %>%
    filter(Method %in% c("random", "positive children", "negative children")) %>%
    mutate(Distance_to_closest = cut(Distance_to_closest, include.lowest = TRUE, breaks = c(0, 0.001, 0.025, 0.05, 0.1, 0.2, 0.35))) %>%
    group_by(Method, Distance_to_closest) %>%
    summarise(Precision = mean(Precision), Recall = mean(Recall), F1 = mean(F1))
) +
  aes(x = Method, y = F1, color = Distance_to_closest, shape = Method) +
  geom_point(size = 5, alpha = 0.85) +
  labs(title = "Children-based ranking", shape = "Ranking", colour = "Distance to the closest", x = "Ranking", y = "F1") +
  geom_line(aes(group = Distance_to_closest), color = "grey20") +
  scale_colour_brewer(palette = "Paired") + 
  theme_cowplot(font_size = 16) + scale_shape_manual(values = c(8, 18, 16)) +
  theme(
    plot.title = element_text(size =  16),
    axis.text.y = element_text(),
    axis.text.x = element_text(angle = 20, vjust = 0.75),
  ) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
p2

p3 <- ggplot(
  scores %>%
    filter(Method %in% c("random", "weighted sum")) %>%
    mutate(Distance_to_closest = cut(Distance_to_closest, include.lowest = TRUE, breaks = c(0, 0.001, 0.025, 0.05, 0.1, 0.2, 0.35))) %>%
    group_by(Method, Distance_to_closest) %>%
    summarise(Precision = mean(Precision), Recall = mean(Recall), F1 = mean(F1))
) +
  aes(x = Method, y = F1, color = Distance_to_closest, shape = Method) +
  geom_point(size = 5, alpha = 0.85) +
  labs(title = "Weighted sum heuristic", shape = "Ranking", colour = "Distance to the closest", x = "Ranking", y = "F1") +
  geom_line(aes(group = Distance_to_closest), color = "grey20") +
  scale_colour_brewer(palette = "Paired") + 
  theme_cowplot(font_size = 16) + scale_shape_manual(values = c(16, 17)) +
  theme(
    plot.title = element_text(size =  16),
    axis.text.y = element_text(),
    axis.text.x = element_text(angle = 20, vjust = 0.75),
  ) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
p3

prow <- plot_grid(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none"),
  p3 + theme(legend.position = "none"),
  ncol = 3,
  rel_heights = c(1, 1, 1),
  vjust = 2,
  hjust = 2
)
prow

legend <- get_legend(
  ggplot(scores %>%
      mutate(Distance_to_closest = cut(Distance_to_closest, include.lowest = TRUE, breaks = c(0, 0.001, 0.025, 0.05, 0.1, 0.2, 0.35))) %>%
      group_by(Method, Distance_to_closest) %>%
      summarise(Precision = mean(Precision), Recall = mean(Recall), F1 = mean(F1))
  ) +
  aes(
      x = reorder(Method, F1),
      y = F1,
      shape = Method,
      color = Distance_to_closest
    ) +
    labs(shape = "Ranking", colour = "Distance to the closest", x = "Ranking", y = "F1") +
    geom_point(size = 5, alpha = 0.85) +
    scale_colour_brewer(palette = "Paired") +
    theme_cowplot(font_size = 15) +
    scale_shape_manual(
      values = c(
        "weighted sum" = 17,
        "positive children" = 18,
        "random" = 16,
        "positive species" = 15,
        "negative species" = 4,
        "negative children" = 8
      )
    ) + theme(legend.box.margin = margin(3, 0, 0, 0)) +
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.direction = "horizontal",
      legend.box = "vertical"
    ) + guides(color = guide_legend(nrow = 1, byrow = TRUE))
)
plot_grid(prow, legend, nrow=2, rel_heights = c(3, 1))