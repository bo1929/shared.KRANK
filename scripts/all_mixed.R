require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)

scores <- read_csv("../results/cscores-10kSpecies_-_combined.csv")

scores$Taxonomic_rank <- factor(
  scores$Taxonomic_rank,
  levels = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
)
scores <- scores %>% filter(Distance_to_closest < 1)
scores$Taxon <- as.factor(scores$Taxon)

ggplot(
  scores,
  aes(x = Distance_to_closest, y = F1, color = Method, linetype = Method, shape = Method)
) +
  facet_wrap(vars(Taxonomic_rank), nrow = 1) +
  geom_point(alpha = 0.15) +
  labs(shape = "Tool", colour = "Tool", linetype = "Tool", x = "Distance to the closest",  y = "F1") +
  stat_smooth(se = F, span = 0.7, method = "glm", method.args = list(family = binomial), size = 1) +
  scale_colour_brewer(palette = "Dark2") +
  scale_x_continuous(limits = c(NA, 0.3)) +
  theme_cowplot(font_size = 17) +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12),) +
  theme(legend.position = "bottom", legend.justification = "center", legend.direction = "vertical") +
  theme(panel.spacing.x = unit(1, "lines"))


ggplot(
  scores %>%
    mutate(
      Distance_to_closest = cut(Distance_to_closest, include.lowest = TRUE, breaks = c(0, 0.001, 0.01, 0.025, 0.05, 0.1, 0.2, 0.35, 1))
      ) %>%
    group_by(Taxonomic_rank, Method, Distance_to_closest) %>%
    summarise(F1 = mean(F1)),
  aes(x = Distance_to_closest, y = F1, color = Method, shape = Method)
) +
  facet_wrap(vars(Taxonomic_rank), nrow = 1) +
  geom_point(size = 3, alpha = 0.85) +
  geom_line(aes(group = Method)) +
  labs(shape = "Tool", colour = "Tool", x = "Distance to the closest", y = "F1") +
  scale_colour_brewer(palette = "Dark2") +
  theme_cowplot(font_size = 17) +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1),
    aspect.ratio = 1.25
  ) +
  theme(legend.position = "bottom", legend.justification = "center", legend.direction = "vertical") +
  theme(panel.spacing.x = unit(1, "lines"))

ggplot(
  scores %>%
    group_by(Taxon, Taxonomic_rank, Method) %>%
    summarise(F1 = mean(F1), Distance_to_closest = mean(Distance_to_closest)),
  aes(x = Method,  y = F1, group = Method, color = Method)
) +
  facet_wrap(vars(Taxonomic_rank), nrow = 2) +
  geom_boxplot(orientation = "x") +
  geom_point(position = "jitter",  alpha = 0.25, size = 1) +
  scale_colour_brewer(palette = "Dark2") +
  theme_cowplot(font_size = 17) +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12), aspect.ratio = 1.25) +
  labs(colour = "Method", x = "Class", y = "F1") +
  theme(panel.spacing.x = unit(1, "lines"), axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  theme(legend.position = "bottom", legend.justification = "center",  legend.direction = "vertical")

ggplot(
  scores,
  aes(x = Reference_genome_in_phylum, y = F1, color = Method, linetype = Method, shape = Method)
) +
  facet_wrap(vars(Taxonomic_rank), nrow = 1) +
  geom_point(alpha = 0.1) +
  labs(shape = "Tool",
    colour = "Tool",
    linetype = "Tool",
    x = "# of reference genomes in phylum",
    y = "F1"
  ) +
  stat_smooth(
    se = F,
    span = 0.7,
    method = "lm",
    formula = y ~ splines::bs(x, 3)
  ) +
  scale_colour_brewer(palette = "Dark2") +
  scale_x_continuous(trans = "log1p", breaks = c(1, 10, 100, 1000)) +
  theme_cowplot(font_size = 17) +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    aspect.ratio = 1.25
  ) +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "vertical"
  ) +
  theme(panel.spacing.x = unit(1, "lines"))


ggplot(
  scores %>% mutate(
    Reference_genome_in_phylum = cut(
      Reference_genome_in_phylum,
      include.lowest = TRUE,
      breaks = c(0, 10, 100, 1000, 10000)
      )
  ),
  aes(x = Reference_genome_in_phylum, y = F1, color = Method,)
) +
  facet_wrap(vars(Taxonomic_rank), nrow = 3) +
  geom_boxplot(aes(x = Reference_genome_in_phylum, y = F1, color = Method)) +
  labs(
    shape = "Tool",
    colour = "Tool",
    linetype = "Tool",
    x = "# of reference genomes in phylum",
    y = "F1"
  ) +
  scale_colour_brewer(palette = "Dark2") +
  theme_cowplot(font_size = 17) +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 90),) +
  theme(
    legend.position = "left",
    legend.justification = "center",
    legend.direction = "vertical"
  ) +
  theme(panel.spacing.x = unit(1, "lines"))


ggplot(
  scores  %>% mutate(
    Reference_genome_in_phylum = cut(
      Reference_genome_in_phylum,
      include.lowest = TRUE,
      breaks = c(0, 10, 100, 1000, 10000)
    )
  ),
  aes(x = Distance_to_closest, y = F1, color = Method, linetype = Method, shape = Method)
) +
  facet_grid(rows = vars(Reference_genome_in_phylum), cols = vars(Taxonomic_rank)) +
  geom_point(alpha = 0.1) +
  stat_smooth(se = F, span = 0.7, method = "glm", method.args = list(family = binomial), size = 1) +
  labs(shape = "Approach", colour = "Approach", linetype = "Approach", x = "Distance to the closest", y = "F1") +
  scale_colour_brewer(palette = "Dark2") +
  theme_cowplot(font_size = 17) +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12, angle = 90),) +
  theme( legend.position = "left", legend.justification = "center", legend.direction = "vertical") +
  theme(panel.spacing.x = unit(1, "lines"))



ggplot(
  scores %>%
    group_by(Taxon, Taxonomic_rank, Method) %>%
    summarise(F1 = mean(F1), Distance_to_closest = mean(Distance_to_closest)),
  aes(x = Distance_to_closest, y = F1, group = Method, color = Method)) +
  facet_wrap(vars(Taxonomic_rank), nrow = 2) +
  geom_point(position = "jitter", alpha = 0.25, size = 1) +
  stat_smooth(se = F, span = 0.7, method = "glm", method.args = list(family = binomial), size = 1) +
  scale_colour_brewer(palette = "Dark2") +
  theme_cowplot(font_size = 17) +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_text(size = 12), axis.text.x = element_text(size = 12), aspect.ratio = 1.25) +
  labs(colour = "Approach", x = "Mean distances to closest of the phylum", y = "F1") +
  theme(panel.spacing.x = unit(1, "lines")) +
  theme(legend.position = "bottom", legend.justification = "center", legend.direction = "vertical")

ggplot(
  scores %>%
    mutate(
      dist_bins = cut(
        Distance_to_closest,
        include.lowest = TRUE,
        breaks = c(0, 0.001, 0.01, 0.025, 0.05, 0.1, 0.2, 0.5, 1)
      )
    ) %>%
    group_by(Taxonomic_rank, Method, dist_bins) %>%
    summarise(F1 = mean(F1))
  ) +
  aes(color=grepl("CONS",Method), y = F1, x = reorder(sub(" k-mers", "", sub("k.* ~ ", "", sub(".0.00.*", "", Method))), F1), shape = Method) +
  facet_wrap(vars(Taxonomic_rank), scales = "free", nrow = 1) +
  geom_point(size = 2.5, alpha = 0.85) +
  labs(shape = "Tool & Memory", colour = "Tool & Memory", x = "Distance to the closest", y = "F1") +
  geom_line(aes(group=dist_bins),color="grey60")+
  scale_colour_brewer(palette = "Set2") +
  theme_cowplot(font_size = 16) +
  theme(legend.text = element_text(size = 12.75), legend.title = element_text(size = 14)) +
  theme(
    axis.text.y = element_text(size = 12.25),
    axis.text.x = element_text(size = 12, angle = 90, vjust = 0.5, hjust = 1)
  ) +
  theme(
    legend.position = "bottom",
    legend.justification = "center",
    legend.direction = "vertical"
  )









pr = function(x,k,p,h=0)
  vapply(x,
         function(di) sum(choose(k,0:p)*di^(0:p)*(1-di)^(k-0:p)*(2*(1-h/k)^(0:p)-((1-h/k)^(0:p))^2)),
         c(1))

pr(0.15,30,4)

(pr(0.1,30,0)*vote(30,0)+
    pr(0.1,30,1)*vote(30,1)+
    pr(0.1,30,2)*vote(30,2)+
    pr(0.1,30,3)*vote(30,3)+
    pr(0.1,30,4)*vote(30,4)+
    pr(0.1,30,5)*vote(30,5))*120/10


ggplot()+stat_function(fun=function(x) pr(x,30,14,5))+
  stat_function(fun=function(x) pr(x,32,15,5),color="red")+
  coord_cartesian(xlim=c(0,.635))

vote = function(k,d) (1-d/k)^k

vote(30,4)*2

ggplot()+stat_function(fun=function(x) vote(30,x*5))+
  stat_function(fun=function(x) vote(32,x*5),color="red")+
  coord_cartesian(xlim=c(0,1))

require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)

taxon_dist <- read_tsv("../data/closest_taxon_wrank.txt")

ggplot(taxon_dist) + aes(x=rank, y=dist, color=rank) +
  geom_boxplot() +
  geom_point()

changed = function(k,d) 1-(1-d)^k
atleasttwosame =  function(k,d,n=n) 1 - (dbinom(n-1,n,changed(k,d)) + dbinom(n,n,changed(k,d)))
sharedkmers = function(k,d,n=n,L=L) (n*L-n*L*atleasttwosame(k,d,n=n))/(n*L)

pr = function(x,k,p,h=0)
  vapply(x,
         function(di) sum(choose(k,0:p)*di^(0:p)*(1-di)^(k-0:p)*(2*(1-h/k)^(0:p)-((1-h/k)^(0:p))^2)),
         c(1))

pr(0.15,30,4)

(pr(0.1,30,0)*vote(30,0)+
    pr(0.1,30,1)*vote(30,1)+
    pr(0.1,30,2)*vote(30,2)+
    pr(0.1,30,3)*vote(30,3)+
    pr(0.1,30,4)*vote(30,4)+
    pr(0.1,30,5)*vote(30,5))*120/10


ggplot()+stat_function(fun=function(x) pr(x,30,14,5))+
  stat_function(fun=function(x) pr(x,32,15,5),color="red")+
  coord_cartesian(xlim=c(0,.635))

vote = function(k,d) (1-d/k)^k

vote(30,4)*2

ggplot()+stat_function(fun=function(x) vote(30,x*5))+
  stat_function(fun=function(x) vote(32,x*5),color="red")+
  coord_cartesian(xlim=c(0,1))
