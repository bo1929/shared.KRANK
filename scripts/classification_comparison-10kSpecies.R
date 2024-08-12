require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)

scores <- read_csv("../results/cscores-10kSpecies_-_combined.csv")

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
scores$Taxon <- as.factor(scores$Taxon)

scores$Method[scores$Method == "k=30 w=33 h=14 b=16 l=2 mixed const. ~ 51.2Gb 4.29B×2 k-mers (0.00)"] <- "KRANK-hs v0.3.2 (51.2GB)"
scores$Method[scores$Method == "k=29 w=32 h=13 b=16 l=2 mixed const. ~ 12.8Gb 1.07B×2 k-mers (0.00)"] <- "KRANK-lw v0.3.2 (12.8GB)"
scores$Method[scores$Method == "Kraken-II"] <- "Kraken-II v2.1.3 (46.5GB)"
scores$Method[scores$Method == "CONSULT-II 140Gb (0.00)"] <- "CONSULT-II v0.4.0 (140.7GB)"
scores$Method[scores$Method == "CLARK"] <- "CLARK v1.2.6.1 (149.6GB)"
scores$Method <- factor(
  scores$Method,
  levels = c(
    "CLARK v1.2.6.1 (149.6GB)",
    "CONSULT-II v0.4.0 (140.7GB)",
    "Kraken-II v2.1.3 (46.5GB)",
    "KRANK-hs v0.3.2 (51.2GB)",
    "KRANK-lw v0.3.2 (12.8GB)"
  )
)
p1 <- ggplot(
  scores %>% filter(
    Method %in% c(
      "CLARK v1.2.6.1 (149.6GB)",
      "CONSULT-II v0.4.0 (140.7GB)",
      "Kraken-II v2.1.3 (46.5GB)",
      "KRANK-hs v0.3.2 (51.2GB)",
      "KRANK-lw v0.3.2 (12.8GB)"
    )
  ) %>%
    mutate(
      Distance_to_closest = cut(
        Distance_to_closest,
        include.lowest = TRUE,
        breaks = c(0, 0.001, 0.025, 0.05, 0.1, 0.2, 0.35, 1)
      )
    ) %>%
    group_by(Method, Taxonomic_rank, Distance_to_closest) %>%
    summarise(
      Recall = mean(Recall),
      F1 = mean(F1),
      Precision = mean(Precision)
    ) %>% filter(!is.na(Distance_to_closest))
) +
  aes(
    y = F1,
    x = Taxonomic_rank,
    shape = Method,
    color = Method
  ) +
  facet_wrap(~ factor(Distance_to_closest),
             scales = "free_y",
             nrow = 3) +
  geom_line(aes(group = Method), size = 1, alpha = 0.75) +
  geom_point(size = 2, alpha = 0.75) +
  labs(
    shape = "Tool & Memory",
    colour = "Tool & Memory",
    x = "",
    y = "F1"
  ) +
  scale_color_manual(values = c("#ff7f00", "#33a02c", "#e31a1c" , "#1f78b4", "#a6cee3")) +
  scale_shape_manual(values = c(18, 17, 15, 16, 16)) +
  theme_cowplot(font_size = 20) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.margin = margin(-30, 0, 0, 0)
  ) + guides(color = guide_legend(nrow = 3, byrow = TRUE))
p1
ggsave2(
  "../figures/classification_comparison-defaultsF1-10kSpecies.pdf",
  width = 6,
  height = 8
)

p3 <- ggplot(
  scores %>% filter(
    Method %in% c(
      "CLARK v1.2.6.1 (149.6GB)",
      "CONSULT-II v0.4.0 (140.7GB)",
      "Kraken-II v2.1.3 (46.5GB)",
      "KRANK-hs v0.3.2 (51.2GB)",
      "KRANK-lw v0.3.2 (12.8GB)"
    )
  ) %>%
    mutate(
      Distance_to_closest = cut(
        Distance_to_closest,
        include.lowest = TRUE,
        breaks = c(0, 0.001, 0.025, 0.05, 0.1, 0.2, 0.35, 1)
      )
    ) %>% group_by(Method, Taxonomic_rank, Distance_to_closest) %>%
    summarise(
      Recall = mean(Recall),
      F1 = mean(F1),
      Precision = mean(Precision)
    )
) +
  aes(
    y = Recall,
    x = Precision,
    color = Method,
    shape = Taxonomic_rank
  ) +
  facet_wrap("Distance_to_closest", scale = "fixed") +
  geom_point(alpha = 0.8, size = 4) +
  labs(
    shape = "Taxonomic Rank",
    colour = "Tool & Memory",
    x = "Precision",
    y = "Recall"
  ) +
  scale_shape_manual(values = c(8, 9, 10, 15, 16, 17, 18)) +
  # scale_color_brewer(palette = "Paired") +
  scale_color_manual(values = c("#ff7f00", "#33a02c", "#e31a1c" , "#1f78b4", "#a6cee3")) +
  theme_cowplot(font_size = 17) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.spacing.x = unit(2, "lines"),
    panel.spacing.x = unit(1, "lines")
  ) + guides(shape = guide_legend(nrow = 3, byrow = TRUE)) + guides(color =
                                                                      guide_legend(nrow = 3, byrow = TRUE))
p3
ggsave(
  "../figures/classification_comparison-defaultsPrecisionRecall-10kSpecies.pdf",
  width = 10,
  height = 12
)

scores <- read_csv("../results/cscores-10kSpecies_-_combined.csv")
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
scores$Taxon <- as.factor(scores$Taxon)
scores$Memory <- 0
scores$Memory[scores$Method == "k=30 w=33 h=14 b=16 l=3 mixed + no const. ~ 76.8.2Gb 4.29B×3 k-mers (0.00)"] <- 76.8
scores$Memory[scores$Method == "k=30 w=33 h=14 b=16 l=2 mixed const. ~ 51.2Gb 4.29B×2 k-mers (0.00)"] <- 51.2
scores$Memory[scores$Method == "k=30 w=33 h=14 b=16 l=1 mixed const. ~ 25.6Gb 4.29B×1 k-mers (0.00)"] <- 25.6
scores$Memory[scores$Method == "k=29 w=32 h=13 b=16 l=2 mixed const. ~ 12.8Gb 1.07B×2 k-mers (0.00)"] <- 12.8
scores$Memory[scores$Method == "k=29 w=32 h=13 b=16 l=1 no const. ~ 6.4Gb 1.07B×1 k-mers (0.00)"] <- 6.4
scores$Memory[scores$Method == "k=28 w=31 h=12 b=16 l=2 mixed const. ~ 3.2Gb 268M×2 k-mers (0.00)"] <- 3.2
scores$Memory[scores$Method == "k=28 w=31 h=12 b=16 l=1 no const. ~ 1.6Gb 268M×1 k-mers (0.00)"] <- 1.6
scores$Memory[scores$Method == "Kraken-II"] <- 46.5
scores$Memory[scores$Method == "Kraken-II (4Gb)"] <- 4
scores$Memory[scores$Method == "Kraken-II (16Gb)"] <- 16
scores$Memory[scores$Method == "CONSULT-II 140Gb (0.00)"] <- 140.7
scores$Memory[scores$Method == "CONSULT-II 16Gb (0.00)"] <- 16
scores$Memory[scores$Method == "CONSULT-II 5Gb (0.00)"] <- 5
scores$Memory[scores$Method == "CONSULT-II 3Gb (0.00)"] <- 3
scores$Memory[scores$Method == "CLARK"] <- 149.6
scores$Method[scores$Method == "k=30 w=33 h=14 b=16 l=3 mixed + no const. ~ 76.8.2Gb 4.29B×3 k-mers (0.00)"] <- "KRANK v0.3.2"
scores$Method[scores$Method == "k=30 w=33 h=14 b=16 l=2 mixed const. ~ 51.2Gb 4.29B×2 k-mers (0.00)"] <- "KRANK v0.3.2"
scores$Method[scores$Method == "k=30 w=33 h=14 b=16 l=1 mixed const. ~ 25.6Gb 4.29B×1 k-mers (0.00)"] <- "KRANK v0.3.2"
scores$Method[scores$Method == "k=29 w=32 h=13 b=16 l=2 mixed const. ~ 12.8Gb 1.07B×2 k-mers (0.00)"] <- "KRANK v0.3.2"
scores$Method[scores$Method == "k=29 w=32 h=13 b=16 l=1 no const. ~ 6.4Gb 1.07B×1 k-mers (0.00)"] <- "KRANK v0.3.2"
scores$Method[scores$Method == "k=28 w=31 h=12 b=16 l=2 mixed const. ~ 3.2Gb 268M×2 k-mers (0.00)"] <- "KRANK v0.3.2"
scores$Method[scores$Method == "k=28 w=31 h=12 b=16 l=1 no const. ~ 1.6Gb 268M×1 k-mers (0.00)"] <- "KRANK v0.3.2"
scores$Method[scores$Method == "Kraken-II"] <- "Kraken-II v2.1.3"
scores$Method[scores$Method == "Kraken-II (4Gb)"] <- "Kraken-II v2.1.3"
scores$Method[scores$Method == "Kraken-II (16Gb)"] <- "Kraken-II v2.1.3"
scores$Method[scores$Method == "CONSULT-II 140Gb (0.00)"] <- "CONSULT-II v0.4.0"
scores$Method[scores$Method == "CONSULT-II 16Gb (0.00)"] <- "CONSULT-II v0.4.0"
scores$Method[scores$Method == "CONSULT-II 5Gb (0.00)"] <- "CONSULT-II v0.4.0"
scores$Method[scores$Method == "CONSULT-II 3Gb (0.00)"] <- "CONSULT-II v0.4.0"
scores$Method[scores$Method == "CLARK"] <- "CLARK v1.2.6.1"
scores$Method <- factor(
  scores$Method,
  levels = c(
    "CLARK v1.2.6.1",
    "Kraken-II v2.1.3",
    "CONSULT-II v0.4.0",
    "KRANK v0.3.2"
  )
)

p2 <- ggplot(
  scores %>% filter(Memory > 0) %>%  # filter(Method %in% c("CONSULT-II v0.4.0", "KRANK v0.3.2")) %>%
    mutate(
      Distance_to_closest = cut(
        Distance_to_closest,
        include.lowest = TRUE,
        breaks = c(0.35, 0.2, 0.1, 0.05, 0.025, 0.001, 0.0)
      )
    ) %>%
    group_by(Method, Memory, Distance_to_closest) %>%
    summarise(
      Recall = mean(Recall),
      F1 = mean(F1),
      Precision = mean(Precision)
    ) %>% filter(!is.na(Distance_to_closest))
) +
  aes(
    y = F1,
    x = Memory,
    shape = reorder(Method, F1),
    color = reorder(Method, F1)
  ) +
  facet_wrap(~ factor(Distance_to_closest),
             scales = "free_y",
             nrow = 3) +
  geom_line(aes(group = Method), size = 1) +
  geom_point(size = 5) +
  labs(
    shape = "Tool",
    colour = "Tool",
    x = "Memory (GB)",
    y = "F1"
  ) +
  scale_color_manual(values = c("#e31a1c", "#ff7f00", "#33a02c", "#1f78b4")) +
  scale_shape_manual(values = c(18, 17, 15, 16)) +
  theme_cowplot(font_size = 20) +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.margin = margin(0, 0, 0, 0)
  ) +
  geom_vline(
    xintercept = 51.2,
    colour = "#333",
    linewidth = 1.2,
    alpha = 0.65,
    linetype = 3
  ) +
  geom_vline(
    xintercept = 12.8,
    colour = "#333",
    linewidth = 1.2,
    alpha = 0.65,
    linetype = 3
  ) +
  scale_x_continuous(trans = "log", breaks = c(4, 16, 64)) +
  guides(shape = guide_legend(nrow = 2, byrow = TRUE)) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))
p2
ggsave2(
  "../figures/classification_comparison-varying_memory-10kSpecies.pdf",
  width = 6,
  height = 8
)

plot_grid(p1, p2, labels = c("A", "B"), label_size = 24)
ggsave2(
  "../figures/classification_comparison-main_new-10kSpecies-vline.pdf",
  width = 13,
  height = 8
)

scores <- read_csv("../results/cscores-10kSpecies_-_combined.csv")
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
scores$Taxon <- as.factor(scores$Taxon)
scores$nkmers <- 0
scores$nkmers[scores$Method == "k=30 w=33 h=14 b=16 l=3 mixed + no const. ~ 76.8.2Gb 4.29B×3 k-mers (0.00)"] <- 4.29 *
  3
scores$nkmers[scores$Method == "k=30 w=33 h=14 b=16 l=2 mixed const. ~ 51.2Gb 4.29B×2 k-mers (0.00)"] <- 4.29 *
  2
scores$nkmers[scores$Method == "k=30 w=33 h=14 b=16 l=1 mixed const. ~ 25.6Gb 4.29B×1 k-mers (0.00)"] <- 4.29
scores$nkmers[scores$Method == "k=29 w=32 h=13 b=16 l=2 mixed const. ~ 12.8Gb 1.07B×2 k-mers (0.00)"] <- 1.07 *
  2
scores$nkmers[scores$Method == "k=29 w=32 h=13 b=16 l=1 no const. ~ 6.4Gb 1.07B×1 k-mers (0.00)"] <- 1.07
scores$nkmers[scores$Method == "k=28 w=31 h=12 b=16 l=2 mixed const. ~ 3.2Gb 268M×2 k-mers (0.00)"] <- 0.268 *
  2
scores$nkmers[scores$Method == "k=28 w=31 h=12 b=16 l=1 no const. ~ 1.6Gb 268M×1 k-mers (0.00)"] <- 0.268
scores$nkmers[scores$Method == "CONSULT-II 140Gb (0.00)"] <- 2 ** 30 * 7 *
  2 / 10 ** 9
scores$nkmers[scores$Method == "CONSULT-II 16Gb (0.00)"] <- 2 ** 26 * 14 *
  2 / 10 ** 9
scores$nkmers[scores$Method == "CONSULT-II 5Gb (0.00)"] <- 2 ** 24 * 16 *
  2 / 10 ** 9
scores$nkmers[scores$Method == "CONSULT-II 3Gb (0.00)"] <- 2 ** 24 * 9 *
  2 / 10 ** 9
scores$Method[scores$Method == "k=30 w=33 h=14 b=16 l=3 mixed + no const. ~ 76.8.2Gb 4.29B×3 k-mers (0.00)"] <- "KRANK v0.3.2"
scores$Method[scores$Method == "k=30 w=33 h=14 b=16 l=2 mixed const. ~ 51.2Gb 4.29B×2 k-mers (0.00)"] <- "KRANK v0.3.2"
scores$Method[scores$Method == "k=30 w=33 h=14 b=16 l=1 mixed const. ~ 25.6Gb 4.29B×1 k-mers (0.00)"] <- "KRANK v0.3.2"
scores$Method[scores$Method == "k=29 w=32 h=13 b=16 l=2 mixed const. ~ 12.8Gb 1.07B×2 k-mers (0.00)"] <- "KRANK v0.3.2"
scores$Method[scores$Method == "k=29 w=32 h=13 b=16 l=1 no const. ~ 6.4Gb 1.07B×1 k-mers (0.00
)"] <- "KRANK v0.3.2"
scores$Method[scores$Method == "k=28 w=31 h=12 b=16 l=2 mixed const. ~ 3.2Gb 268M×2 k-mers (0.00)"] <- "KRANK v0.3.2"
scores$Method[scores$Method == "k=28 w=31 h=12 b=16 l=1 no const. ~ 1.6Gb 268M×1 k-mers (0.00)"] <- "KRANK v0.3.2"
scores$Method[scores$Method == "CONSULT-II 140Gb (0.00)"] <- "CONSULT-II v0.4.0"
scores$Method[scores$Method == "CONSULT-II 16Gb (0.00)"] <- "CONSULT-II v0.4.0"
scores$Method[scores$Method == "CONSULT-II 5Gb (0.00)"] <- "CONSULT-II v0.4.0"
scores$Method[scores$Method == "CONSULT-II 3Gb (0.00)"] <- "CONSULT-II v0.4.0"
scores$Method <- factor(scores$Method, levels = c("CONSULT-II v0.4.0", "KRANK v0.3.2"))

p1 <- ggplot(
  scores %>% filter(nkmers > 0) %>%
    mutate(
      Distance_to_closest = cut(
        Distance_to_closest,
        include.lowest = TRUE,
        breaks = c(0, 0.001, 0.025, 0.05, 0.1, 0.2, 0.35, 1)
      )
    ) %>%
    group_by(Method, nkmers, Distance_to_closest) %>%
    summarise(
      Recall = mean(Recall),
      F1 = mean(F1),
      Precision = mean(Precision)
    )
) +
  aes(
    y = F1,
    x = nkmers,
    shape = Method,
    color = Method
  ) +
  facet_wrap(vars(Distance_to_closest),
             scales = "free_y",
             nrow = 1) +
  geom_line(aes(group = Method), size = 2, alpha = 0.9) +
  geom_point(size = 4, alpha = 0.9) +
  labs(
    shape = "Tool",
    colour = "Tool",
    x = "Total number of k-mers in the library (billion)",
    y = "F1"
  ) +
  scale_color_manual(values = c("#33a02c", "#1f78b4")) +
  scale_shape_manual(values = c(18, 17, 15, 16, 16)) +
  theme_cowplot(font_size = 17) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.margin = margin(0, 0, 0, 0)
  ) + guides(color = guide_legend(nrow = 1, byrow = TRUE)) +
  scale_x_continuous(trans = "log", breaks = c(0.1, 0.5, 1, 2, 4, 8))
p1
ggsave2(
  "../figures/classification_comparison-varying_numkmers-10kSpecies.pdf",
  width = 6,
  height = 8
)

null <- read_tsv("../data/query_ranks.tsv", col_names = FALSE)
c(
  sum(null$X2[null$X2 != 0] == tail(names(sort(
    table(null$X2[null$X2 != 0])
  )), 1)) / sum(null$X2 != 0),
  sum(null$X3[null$X3 != 0] == tail(names(sort(
    table(null$X3[null$X3 != 0])
  )), 1)) / sum(null$X3 != 0),
  sum(null$X4[null$X4 != 0] == tail(names(sort(
    table(null$X4[null$X4 != 0])
  )), 1)) / sum(null$X4 != 0),
  sum(null$X5[null$X5 != 0] == tail(names(sort(
    table(null$X5[null$X5 != 0])
  )), 1)) / sum(null$X5 != 0),
  sum(null$X6[null$X6 != 0] == tail(names(sort(
    table(null$X6[null$X6 != 0])
  )), 1)) / sum(null$X6 != 0),
  sum(null$X7[null$X7 != 0] == tail(names(sort(
    table(null$X7[null$X7 != 0])
  )), 1)) / sum(null$X7 != 0),
  sum(null$X8[null$X8 != 0] == tail(names(sort(
    table(null$X8[null$X8 != 0])
  )), 1)) / sum(null$X8 != 0)
)

n <- c(
  sum(null$X2[null$X2 != 0] == tail(names(sort(
    table(null$X2[null$X2 != 0])
  )), 1)) / (sum(null$X2[null$X2 != 0] == tail(names(
    sort(table(null$X2[null$X2 != 0]))
  ), 1)) + 0.5 * sum(null$X2[null$X2 != 0] != tail(names(
    sort(table(null$X2[null$X2 != 0]))
  ), 1))),
  sum(null$X3[null$X3 != 0] == tail(names(sort(
    table(null$X3[null$X3 != 0])
  )), 1)) / (sum(null$X3[null$X3 != 0] == tail(names(
    sort(table(null$X3[null$X3 != 0]))
  ), 1)) + 0.5 * sum(null$X3[null$X3 != 0] != tail(names(
    sort(table(null$X2[null$X3 != 0]))
  ), 1))),
  sum(null$X4[null$X4 != 0] == tail(names(sort(
    table(null$X4[null$X4 != 0])
  )), 1)) / (sum(null$X4[null$X4 != 0] == tail(names(
    sort(table(null$X4[null$X4 != 0]))
  ), 1)) + 0.5 * sum(null$X4[null$X4 != 0] != tail(names(
    sort(table(null$X4[null$X4 != 0]))
  ), 1))),
  sum(null$X5[null$X5 != 0] == tail(names(sort(
    table(null$X5[null$X5 != 0])
  )), 1)) / (sum(null$X5[null$X5 != 0] == tail(names(
    sort(table(null$X5[null$X5 != 0]))
  ), 1)) + 0.5 * sum(null$X5[null$X5 != 0] != tail(names(
    sort(table(null$X5[null$X5 != 0]))
  ), 1))),
  sum(null$X6[null$X6 != 0] == tail(names(sort(
    table(null$X6[null$X6 != 0])
  )), 1)) / (sum(null$X6[null$X6 != 0] == tail(names(
    sort(table(null$X6[null$X6 != 0]))
  ), 1)) + 0.5 * sum(null$X6[null$X6 != 0] != tail(names(
    sort(table(null$X6[null$X6 != 0]))
  ), 1))),
  sum(null$X7[null$X7 != 0] == tail(names(sort(
    table(null$X7[null$X7 != 0])
  )), 1)) / (sum(null$X7[null$X7 != 0] == tail(names(
    sort(table(null$X7[null$X7 != 0]))
  ), 1)) + 0.5 * sum(null$X7[null$X7 != 0] != tail(names(
    sort(table(null$X7[null$X7 != 0]))
  ), 1))),
  sum(null$X8[null$X8 != 0] == tail(names(sort(
    table(null$X8[null$X8 != 0])
  )), 1)) / (sum(null$X8[null$X8 != 0] == tail(names(
    sort(table(null$X8[null$X8 != 0]))
  ), 1)) + 0.5 * sum(null$X8[null$X8 != 0] != tail(names(
    sort(table(null$X8[null$X8 != 0]))
  ), 1)))
)
n <- c(
  sum(null$X2 == tail(names(sort(
    table(null$X2[null$X2 != 0])
  )), 1)) / (sum(null$X2 == tail(names(
    sort(table(null$X2[null$X2 != 0]))
  ), 1)) + 0.5 * sum(null$X2 != tail(names(
    sort(table(null$X2[null$X2 != 0]))
  ), 1))),
  sum(null$X3 == tail(names(sort(
    table(null$X3[null$X3 != 0])
  )), 1)) / (sum(null$X3 == tail(names(
    sort(table(null$X3[null$X3 != 0]))
  ), 1)) + 0.5 * sum(null$X3 != tail(names(
    sort(table(null$X2[null$X3 != 0]))
  ), 1))),
  sum(null$X4 == tail(names(sort(
    table(null$X4[null$X4 != 0])
  )), 1)) / (sum(null$X4 == tail(names(
    sort(table(null$X4[null$X4 != 0]))
  ), 1)) + 0.5 * sum(null$X4 != tail(names(
    sort(table(null$X4[null$X4 != 0]))
  ), 1))),
  sum(null$X5 == tail(names(sort(
    table(null$X5[null$X5 != 0])
  )), 1)) / (sum(null$X5 == tail(names(
    sort(table(null$X5[null$X5 != 0]))
  ), 1)) + 0.5 * sum(null$X5 != tail(names(
    sort(table(null$X5[null$X5 != 0]))
  ), 1))),
  sum(null$X6 == tail(names(sort(
    table(null$X6[null$X6 != 0])
  )), 1)) / (sum(null$X6 == tail(names(
    sort(table(null$X6[null$X6 != 0]))
  ), 1)) + 0.5 * sum(null$X6 != tail(names(
    sort(table(null$X6[null$X6 != 0]))
  ), 1))),
  sum(null$X7 == tail(names(sort(
    table(null$X7[null$X7 != 0])
  )), 1)) / (sum(null$X7 == tail(names(
    sort(table(null$X7[null$X7 != 0]))
  ), 1)) + 0.5 * sum(null$X7 != tail(names(
    sort(table(null$X7[null$X7 != 0]))
  ), 1))),
  sum(null$X8 == tail(names(sort(
    table(null$X8[null$X8 != 0])
  )), 1)) / (sum(null$X8 == tail(names(
    sort(table(null$X8[null$X8 != 0]))
  ), 1)) + 0.5 * sum(null$X8 != tail(names(
    sort(table(null$X8[null$X8 != 0]))
  ), 1)))
)
n

scores <- read_csv("../results/cscores-10kSpecies_-_combined.csv")
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
scores$Taxon <- as.factor(scores$Taxon)

scores$Method[scores$Method == "k=30 w=33 h=14 b=16 l=2 mixed const. ~ 51.2Gb 4.29B×2 k-mers (0.00)"] <- "KRANK-hs v0.3.2 (51.2Gb)"
scores$Method[scores$Method == "k=29 w=32 h=13 b=16 l=2 mixed const. ~ 12.8Gb 1.07B×2 k-mers (0.00)"] <- "KRANK-lw v0.3.2 (12.8Gb)"
scores$Method[scores$Method == "Kraken-II"] <- "Kraken-II v2.1.3 (46.5Gb)"
scores$Method[scores$Method == "CONSULT-II 140Gb (0.00)"] <- "CONSULT-II v0.4.0 (140.7Gb)"
scores$Method[scores$Method == "CLARK"] <- "CLARK v1.2.6.1 (149.6Gb)"
scores$Method <- factor(
  scores$Method,
  levels = c(
    "CLARK v1.2.6.1 (149.6Gb)",
    "CONSULT-II v0.4.0 (140.7Gb)",
    "Kraken-II v2.1.3 (46.5Gb)",
    "KRANK-hs v0.3.2 (51.2Gb)",
    "KRANK-lw v0.3.2 (12.8Gb)"
  )
)

p1n <- ggplot(
  scores %>% filter(
    Method %in% c(
      "CLARK v1.2.6.1 (149.6Gb)",
      "CONSULT-II v0.4.0 (140.7Gb)",
      "Kraken-II v2.1.3 (46.5Gb)",
      "KRANK-hs v0.3.2 (51.2Gb)",
      "KRANK-lw v0.3.2 (12.8Gb)"
    )
  ) %>%
    group_by(Method, Taxonomic_rank) %>%
    summarise(
      Recall = mean(Recall),
      F1 = mean(F1),
      Precision = mean(Precision)
    )
) +
  aes(y = F1, x = Taxonomic_rank, color = Method) +
  # facet_wrap(vars(Distance_to_closest), scales = "fixed", nrow = 3) +
  geom_line(aes(shape = Method, group = Method),
            size = 1,
            alpha = 0.75) +
  geom_point(size = 2, alpha = 0.75) +
  labs(
    shape = "Tool & Memory",
    colour = "Tool & Memory",
    x = "",
    y = "F1"
  ) +
  scale_color_manual(values = c("#ff7f00", "#33a02c", "#e31a1c" , "#1f78b4", "#a6cee3")) +
  scale_shape_manual(values = c(18, 17, 15, 16, 16)) +
  theme_cowplot(font_size = 16) +
  theme(axis.text.x = element_text(
    angle = 90,
    vjust = 0.5,
    hjust = 1
  )) +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.margin = margin(-30, 0, 0, 0)
  ) + guides(color = guide_legend(nrow = 3, byrow = TRUE)) +
  geom_line(
    aes(x = p, y = c, group = 1),
    data = data.frame(
      p = c(
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species"
      ),
      c = n
    ),
    color = "black",
    linetype = 2
  ) +
  geom_point(
    aes(x = p, y = c, group = 1),
    data = data.frame(
      p = c(
        "kingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species"
      ),
      c = n
    ),
    color = "black",
    shape = 17
  )
p1n
ggsave2(
  "../figures/classification_comparison-null_model-10kSpecies.pdf",
  width = 7,
  height = 6
)
