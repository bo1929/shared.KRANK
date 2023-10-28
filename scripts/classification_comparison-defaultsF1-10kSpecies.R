require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)

scores <- read_csv("../results/cscores-10kSpecies_-_combined.csv")

scores$Taxonomic_rank <- factor(
  scores$Taxonomic_rank,
  levels = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
)
scores <- scores %>% filter(Distance_to_closest < 0.35)
# scores <- scores %>% filter(Distance_to_closest > 0.00)
scores$Taxon <- as.factor(scores$Taxon)

scores$Method[scores$Method == "k=30 w=33 h=14 b=16 l=2 mixed const. ~ 51.2Gb 4.29B×2 k-mers (0.00)"]<- "KRANK 51.2Gb"
scores$Method[scores$Method == "k=29 w=32 h=13 b=16 l=2 mixed const. ~ 12.8Gb 1.07B×2 k-mers (0.00)"] <- "KRANK 12.8Gb"
scores$Method[scores$Method == "k=29 w=32 h=13 b=16 l=2 mixed const. ~ 12.8Gb 1.07B×2 k-mers (0.00)"] <- "KRANK 12.8Gb"
scores$Method[scores$Method == "Kraken-II"] <- "Kraken-II 46.5Gb"
scores$Method[scores$Method == "CONSULT-II 140Gb (0.00)"] <- "CONSULT-II 140.7Gb"
scores$Method[scores$Method == "CLARK"] <- "CLARK 149.6Gb"

scores$Method <- factor(
  scores$Method,
  levels = c("CLARK 149.6Gb", "Kraken-II 46.5Gb", "CONSULT-II 140.7Gb", "KRANK 51.2Gb", "KRANK 12.8Gb")
)
ggplot(
  scores %>% filter(
    Method %in% c(
      "CLARK 149.6Gb",
      "CONSULT-II 140.7Gb",
      "Kraken-II 46.5Gb",
      "KRANK 12.8Gb",
      "KRANK 51.2Gb"
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
    summarise(F1 = mean(F1))
) +
  aes(y = F1, x = Taxonomic_rank, shape=Method, color=Method) +
  facet_wrap(vars(Distance_to_closest), scales = "fixed", nrow = 2) +
  geom_line(aes(group=Method),size = 1, alpha = 0.75) +
  geom_point(size = 2, alpha = 0.75) +
  labs(shape = "Tool & Memory", colour = "Tool & Memory", x = "", y = "F1") +
  scale_colour_manual(values=c("#ff7f00", "#e31a1c" , "#33a02c","#1f78b4", "#a6cee3")) +
  scale_shape_manual(values=c(18, 17,15,16,16)) +
  theme_cowplot(font_size = 20) +
  theme(
    axis.text.y = element_text(size = 16),
    axis.text.x = element_text(size = 16, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(   legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "vertical",
    legend.margin=margin(-15,-10,0,-10),
    legend.spacing.x = unit(2, "lines"),
    legend.box.margin = margin(-20,-10,0,-10),
    panel.spacing.x = unit(1, "lines"),
    legend.text = element_text(size = 19),
    legend.title = element_text(size = 20),
  ) + guides(color=guide_legend(nrow=2,byrow=TRUE))
ggsave2("../figures/classification_comparison-defaultsF1-10kSpecies.pdf", width=13, height = 8)