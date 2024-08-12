require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)

scores <- read_csv("../results/cscores-10kSpecies_-_combined.csv")

scores$Taxonomic_rank <- factor(
  scores$Taxonomic_rank,
  levels = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
)
scores <- scores %>% filter(Distance_to_closest < 0.35)
scores$Taxon <- as.factor(scores$Taxon)

scores$Method[scores$Method == "CONSULT-II 5Gb (0.00)"] <- "CONSULT-II (v0.4.0) with 5Gb ~ 268Mx2 k-mers"
scores$Method[scores$Method == "CONSULT-II 3Gb (0.00)"] <- "CONSULT-II (v0.4.0) with 3Gb ~ 150Mx2 k-mers"
scores$Method[scores$Method == "k=28 w=31 h=12 b=16 l=2 mixed const. ~ 3.2Gb 268M×2 k-mers (0.00)"]<- "KRANK (v0.3.2) with 3.2Gb ~ 268M×2 k-mers"

ggplot(
  scores %>% filter(
    Method %in% c(
      "CONSULT-II (v0.4.0) with 5Gb ~ 268Mx2 k-mers",
      "CONSULT-II (v0.4.0) with 3Gb ~ 150Mx2 k-mers",
      "KRANK (v0.3.2) with 3.2Gb ~ 268M×2 k-mers"
    )
  ) %>% 
    mutate(
      Distance_to_closest = cut(
        Distance_to_closest,
        include.lowest = TRUE,
        breaks = c(0, 0.001, 0.025, 0.05, 0.1, 0.2, 0.35, 1)
      )
    )  %>% group_by(Method, Distance_to_closest, Taxonomic_rank) %>%
    summarise(Recall = mean(Recall), F1 = mean(F1), Precision = mean(Precision))
) + facet_wrap("Distance_to_closest", scale="free_y", nrow = 2) + 
  geom_point(aes(x=Taxonomic_rank, y=F1, color=Method, shape=Method), size=3, alpha=0.9) +
  geom_line(aes(x=Taxonomic_rank, y=F1, color=Method, group=Method)) +
  scale_colour_manual(values=c("#00a02c", "#99a02c", "#1078b4")) +
  scale_shape_manual(values=c(17, 15, 16)) +
  theme_cowplot(font_size = 18) +
  theme(legend.title = element_blank()) +
  guides(color=guide_legend(nrow=3,byrow=TRUE)) +
  labs(x="Taxonomic Rank") +
  theme(legend.position = "bottom",
        legend.justification = "center",
        legend.direction = "vertical"
  ) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("../figures/classification_comparison_-_withCONSULT-II.pdf", width = 10, height = 10)