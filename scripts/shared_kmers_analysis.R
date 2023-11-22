require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)
require(scales)
require(reshape2)
require(stringr)   
require(pals)

taxon_dist <- read_csv("../data/dist_wrt_lastcommonrank.csv")

taxon_dist$rank <- factor(
  taxon_dist$rank,
  levels = c("root", "kingdom", "phylum", "class", "order", "family", "genus", "species")
)
# taxon_dist <- melt(taxon_dist, id = c("g1", "g2", "rank"), measure = c("JI", "d")) 

p1 <- ggplot(taxon_dist %>% filter(!rank %in% c("root","kingdom")), aes(x=rank, y=d)) +
  geom_boxplot() +
  theme_cowplot(font_size = 16) +
  scale_color_brewer(palette = "Paired")+
  scale_y_continuous(labels = percent, name="Pairwise distance")+
  stat_summary(color="blue")+
  scale_x_discrete(name="",labels = c("Phylum\nClass","Class\nOrder","Order\nFamily","Family\nGenus","Genus\nSpecies")) +
  theme(
    axis.text.y = element_text(size = 14.5),
    axis.text.x = element_text(size = 14.5),
    axis.title.x = element_blank()
  )
p1 <- ggplot(taxon_dist %>% filter(!rank %in% c("root","kingdom")), aes(color=rank, x=d)) +
  stat_ecdf() +
  theme_cowplot(font_size = 16) +
  scale_y_continuous(labels = percent, name = "ECDF") +
  scale_color_discrete(
    type = palette.colors(palette = "Paired", n = 7)[-1], 
    name="Taxonomic rank",
    labels = c("same phylum, diff. class","same class, diff. order", "same order, diff. family","same family, diff. genus","same genus, diff. species", "same species")) +
  theme(
    axis.text.y = element_text(size = 14.5),
    axis.text.x = element_text(size = 14.5),
  ) +
  geom_vline(xintercept = c(5,10,15,20,25,33)/100, color="grey70", linetype=3)+
  scale_x_continuous(name="Pairwise distance") +
  labs(x = "Pairwise distance")
p1

p2 <- ggplot(taxon_dist %>% filter(rank !="species" & !rank %in% c("root","kingdom")), aes(x=rank, y=JI)) +
  geom_boxplot() +
  theme_cowplot(font_size = 16) +
  scale_color_brewer(palette = "Paired")+
  scale_y_continuous(labels = percent, name="Shared 30-mers (portion)", trans = "log2")+
  stat_summary(color="blue")+
  scale_x_discrete(name="",labels = c("Phylum\nClass","Class\nOrder","Order\nFamily","Family\nGenus","Genus\nSpecies")) +
  theme(
    axis.text.y = element_text(size = 14.5),
    axis.text.x = element_text(size = 14.5),
    axis.title.x = element_blank()
    )
p2 <- ggplot(taxon_dist %>% filter(!rank %in% c("root", "kingdom")), aes(color=rank, x=JI)) +
  stat_ecdf() +
  theme_cowplot(font_size = 16) +
  scale_y_continuous(labels = percent, name = "ECDF") +
  scale_color_discrete(
    type = palette.colors(palette = "Paired", n = 7)[-1], 
    name="Taxonomic rank",
    labels = c("same phylum, diff. class","same class, diff. order", "same order, diff. family","same family, diff. genus","same genus, diff. species", "same species")) +
  theme(
    axis.text.y = element_text(size = 14.5),
    axis.text.x = element_text(size = 14.5),
  ) +
  scale_x_continuous(name="Shared 30-mers (portion)", trans = "log10") +
  labs(x = "Shared 30-mers (portion)")
p2

L=10^6
n=32
k=30
prob_change = function(k,d) 1-(1-d)^k
shared_kmers = function(k,d,n=n,L=L)  L*n - L*n*(prob_change(k,d)+((1-prob_change(k,d))*prob_change(k,d)^(n-1)))
shared_kmer_prop = function(k,d,n=n,L=L) shared_kmers(k,d,n=n,L=L)/(n*L)
p3 <- ggplot(aes(x=x), data=data.frame(x=1:1000))+
  stat_function(aes(color="5%"),fun=function(x) shared_kmer_prop(30,0.025,x+1,10^6))+
  stat_function(aes(color="10%"),fun=function(x) shared_kmer_prop(30,0.05,x+1,10^6))+
  stat_function(aes(color="15%"),fun=function(x) shared_kmer_prop(30,0.075,x+1,10^6))+
  stat_function(aes(color="20%"),fun=function(x) shared_kmer_prop(30,0.1,x+1,10^6))+
  stat_function(aes(color="25%"),fun=function(x) shared_kmer_prop(30,0.125,x+1,10^6))+
  stat_function(aes(color="33%"),fun=function(x) shared_kmer_prop(30,1/6,x+1,10^6))+
  theme_cowplot(font_size = 16) +
  labs(y="Expected shared 30-mers (portion)", x="Number of reference genomes", color="Within group\ndiversity") +
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(labels = percent,trans = "log10") +
  theme(
    axis.text.y = element_text(size = 14.5),
    axis.text.x = element_text(size = 14.5),
    legend.position=c(0.7,0.25)
  ) + scale_color_manual(
    values = c("#d73027", "#f46d43", "#fdae61", "#fee090",  "#74add1", "#4575b4"),
    breaks=c("5%", "10%", "15%", "20%", "25%", "33%")
    )
p3

p <- plot_grid(p1 + theme(legend.position = "none"), p2 + theme(legend.position = "none"), rel_widths = c(1, 1), labels = c('A', 'B'), label_size = 19)
legend_p <- get_legend(p1 +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "horizontal",
    legend.box = "vertical"
  ) + guides(color = guide_legend(nrow = 2, byrow = TRUE))
)
p <- plot_grid(p, legend_p, nrow = 2, rel_heights = c(4, 1))
c1 <- plot_grid(p, p3, ncol = 2, labels = c(' ', 'C'), rel_widths = c(2, 1), label_size = 19)
c1
ggsave2("../figures/shared_kmers_analysis.pdf", width=18, height = 6)

# COMPARISON
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
    # mutate(Distance_to_closest = cut(Distance_to_closest, include.lowest = TRUE, breaks = c(0, 0.001, 0.025, 0.05, 0.1, 0.2, 0.35))) %>%
    group_by(Method, Taxonomic_rank) %>%
    summarise(Precision = mean(Precision), Recall = mean(Recall), F1 = mean(F1))
) +
  aes(x = Method, y = F1, color = Taxonomic_rank, shape = Method) +
  geom_point(size = 5, alpha = 0.85) +
  labs(title = "Species-based ranking", shape = "Ranking", colour = "Taxonomic rank", x = "Ranking", y = "F1") +
  geom_line(aes(group = Taxonomic_rank), color = "grey20") +
  scale_colour_brewer(palette = "Paired") + 
  theme_cowplot(font_size = 16) + scale_shape_manual(values = c(4, 15, 16)) +
  theme(
    plot.title = element_text(size =  16),
    axis.text.y = element_text(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 20, vjust = 0.75),
  ) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))

p2 <- ggplot(
  scores %>%
    filter(Method %in% c("random", "positive children", "negative children")) %>%
    # mutate(Distance_to_closest = cut(Distance_to_closest, include.lowest = TRUE, breaks = c(0, 0.001, 0.025, 0.05, 0.1, 0.2, 0.35))) %>%
    group_by(Method, Taxonomic_rank) %>%
    summarise(Precision = mean(Precision), Recall = mean(Recall), F1 = mean(F1))
) +
  aes(x = Method, y = F1, color = Taxonomic_rank, shape = Method) +
  geom_point(size = 5, alpha = 0.85) +
  labs(title = "Children-based ranking", shape = "Ranking", colour = "Taxonomic rank", x = "Ranking", y = "F1") +
  geom_line(aes(group = Taxonomic_rank), color = "grey20") +
  scale_colour_brewer(palette = "Paired") + 
  theme_cowplot(font_size = 16) + scale_shape_manual(values = c(8, 18, 16)) +
  theme(
    plot.title = element_text(size =  16),
    axis.text.y = element_text(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 20, vjust = 0.75),
  ) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))

p3 <- ggplot(
  scores %>%
    filter(Method %in% c("random", "weighted sum")) %>%
    # mutate(Distance_to_closest = cut(Distance_to_closest, include.lowest = TRUE, breaks = c(0, 0.001, 0.025, 0.05, 0.1, 0.2, 0.35))) %>%
    group_by(Method, Taxonomic_rank) %>%
    summarise(Precision = mean(Precision), Recall = mean(Recall), F1 = mean(F1))
) +
  aes(x = Method, y = F1, color = Taxonomic_rank, shape = Method) +
  geom_point(size = 5, alpha = 0.85) +
  labs(title = "Weighted sum heuristic", shape = "Ranking", colour = "Taxonomic rank", x = "Ranking", y = "F1") +
  geom_line(aes(group = Taxonomic_rank), color = "grey20") +
  scale_colour_brewer(palette = "Paired") + 
  theme_cowplot(font_size = 16) + scale_shape_manual(values = c(16, 17)) +
  theme(
    plot.title = element_text(size =  16),
    axis.text.y = element_text(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 20, vjust = 0.75),
  ) + 
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))

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
           group_by(Method, Taxonomic_rank) %>%
           summarise(Precision = mean(Precision), Recall = mean(Recall), F1 = mean(F1))
  ) +
    aes(
      x = reorder(Method, F1),
      y = F1,
      shape = Method,
      color = Taxonomic_rank
    ) +
    labs(shape = "Ranking", colour = "Taxonomic rank", x = "Ranking", y = "F1") +
    geom_point(size = 5, alpha = 0.85) +
    scale_colour_brewer(palette = "Paired") +
    theme_cowplot(font_size = 17) +
    scale_shape_manual(
      values = c(
        "weighted sum" = 17,
        "positive children" = 18,
        "random" = 16,
        "positive species" = 15,
        "negative species" = 4,
        "negative children" = 8
      )
    ) + theme(legend.box.margin = margin(0, 0, 0, 0)) +
    theme(
      legend.position = "bottom",
      legend.justification = "center",
      legend.direction = "horizontal",
      legend.box = "vertical"
    ) + guides(color = guide_legend(nrow = 1, byrow = TRUE))
)
c2 <- plot_grid(prow, legend, nrow=2, rel_heights = c(3, 0.75))

c3 <- ggplot(
  scores %>%
    filter(Method %in% c("random", "weighted sum")) %>%
    mutate(Distance_to_closest = cut(Distance_to_closest, include.lowest = TRUE, breaks = c(0, 0.001, 0.025, 0.05, 0.1, 0.2, 0.35))) %>%
    group_by(Method, Distance_to_closest, Taxonomic_rank) %>%
    summarise(Precision = mean(Precision), Recall = mean(Recall), F1 = mean(F1))
) +
  aes(x = Precision, y = Recall, color = Taxonomic_rank, shape = Method) +
  facet_wrap("Distance_to_closest", scale="free") +
  geom_point(size = 3.2, alpha = 0.85) +
  labs(shape = "Ranking", colour = "Distance to the closest", x = "Precision", y = "Recall") +
  geom_line(aes(group = Taxonomic_rank), color = "grey60") +
  scale_colour_brewer(palette = "Paired") + 
  theme_cowplot(font_size = 16) + scale_shape_manual(values = c(16, 17)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  theme(legend.position = "none")
c3

plot_grid(c1, plot_grid(c2, c3, ncol=2, rel_widths = c(3, 2), labels = c("D", "E"), label_size = 19), nrow=2, rel_heights = c(1, 1.2))
ggsave2("../figures/shared_kmers_analysis-with_comparison.pdf", width=16.5, height = 12)
