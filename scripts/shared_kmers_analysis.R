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
p1
p1 <- ggplot(taxon_dist %>% filter(!rank %in% c("root","kingdom")), aes(color=rank, x=d)) +
  stat_ecdf(linewidth = 1.5, n = 2000) +
  theme_cowplot(font_size = 23) +
  scale_y_continuous(labels = percent, name = "ECDF") +
  scale_color_discrete(
    type = palette.colors(palette = "Paired", n = 7)[-1], 
    name="Rank",
    labels = c("same phylum, diff. class","same class, diff. order", "same order, diff. family","same family, diff. genus","same genus, diff. species", "same species")) +
  geom_vline(xintercept = c(5,10,15,20,25,33)/100, color="grey50", linetype=3)+
  scale_x_continuous(name="Pairwise distance") +
  labs(x = "Pairwise distance")
p1

ggplot(taxon_dist %>% filter(!rank %in% c("root","kingdom")) %>% filter(rank == "species"), aes(color=rank, x=d)) +
  stat_bin() +
  theme_cowplot(font_size = 16) +
  scale_y_continuous(labels = percent, name = "ECDF") +
  scale_color_discrete(
    type = palette.colors(palette = "Paired", n = 7)[-1], 
    name="Rank",
    labels = c("same phylum, diff. class","same class, diff. order", "same order, diff. family","same family, diff. genus","same genus, diff. species", "same species")) +
  theme(
    axis.text.y = element_text(size = 14.5),
    axis.text.x = element_text(size = 14.5),
  ) +
  # geom_vline(xintercept = c(5,10,15,20,25,33)/100, color="grey50", linetype=3)+
  scale_x_continuous(name="Pairwise distance") +
  labs(x = "Pairwise distance")

p2 <- ggplot(taxon_dist %>% filter(!rank %in% c("root","kingdom")), aes(x=rank, y=JI)) +
  geom_boxplot() +
  theme_cowplot(font_size = 16) +
  scale_color_brewer(palette = "Paired")+
  scale_y_continuous(labels = percent, name="Shared 30-mers", trans = "log2")+
  stat_summary(color="blue")+
  scale_x_discrete(name="",labels = c("Phylum\nClass","Class\nOrder","Order\nFamily","Family\nGenus","Genus\nSpecies")) +
  theme(
    axis.text.y = element_text(size = 14.5),
    axis.text.x = element_text(size = 14.5),
    axis.title.x = element_blank()
    )
p2 <- ggplot(taxon_dist %>% filter(!rank %in% c("root", "kingdom")), aes(color=rank, x=JI)) +
  stat_ecdf(linewidth = 1.5, n = 2000) +
  theme_cowplot(font_size = 23) +
  scale_y_continuous(labels = percent, name = "ECDF") +
  scale_color_discrete(
    type = palette.colors(palette = "Paired", n = 7)[-1], 
    name="",
    labels = c("same phylum, diff. class","same class, diff. order", "same order, diff. family","same family, diff. genus","same genus, diff. species", "same species")) +
  scale_x_continuous(name="Shared 30-mers between pair", trans = "log10", breaks = c(5*10^-5, 5*10^-3, 5*10^-2, 5*10^-1), labels = c("0.0005%", "0.05%", "5%", "50%")) +
  theme(legend.position = "bottom") +  guides(color = guide_legend(nrow = 3)) + 
  geom_vline(xintercept = c(5, 60)/100, color="grey50", linetype=3) +
  geom_hline(yintercept = c(14, 50)/100, color="grey50", linetype=3)
p2

L=10^6
n=32
k=30
prob_change = function(k,d) 1-(1-d)^k
shared_kmers = function(k,d,n=n,L=L)  L*n - L*n*(prob_change(k,d)+((1-prob_change(k,d))*prob_change(k,d)^(n-1)))
shared_kmer_prop = function(k,d,n=n,L=L) (1-d)^k*(1-(1-(1-d)^k)^(n-1))
p3 <- ggplot(aes(x=x), data=data.frame(x=1:1000))+
  stat_function(aes(color="5%"), linewidth = 1.2,fun=function(x) shared_kmer_prop(30,0.025,x+1,10^6))+
  stat_function(aes(color="10%"), linewidth = 1.2,fun=function(x) shared_kmer_prop(30,0.05,x+1,10^6))+
  stat_function(aes(color="15%"), linewidth = 1.2,fun=function(x) shared_kmer_prop(30,0.075,x+1,10^6))+
  stat_function(aes(color="20%"), linewidth = 1.2,fun=function(x) shared_kmer_prop(30,0.1,x+1,10^6))+
  stat_function(aes(color="25%"), linewidth = 1.2,fun=function(x) shared_kmer_prop(30,0.125,x+1,10^6))+
  stat_function(aes(color="33%"), linewidth = 1.2,fun=function(x) shared_kmer_prop(30,1/6,x+1,10^6))+
  theme_cowplot(font_size = 23) +
  labs(y="Expected shared 30-mers", x="Number of reference genomes", color="Within group\ndiversity") +
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(labels = percent,trans = "log10") +
  theme(
    legend.position=c(0.64,0.35), 
  ) + scale_color_manual(
    values = c( "#fee090", "#fdae61", "#f46d43", "#d73027", "#B2DF8A", "#4575b4"),
    breaks=c("5%", "10%", "15%", "20%", "25%", "33%")
    )
p3

p <- plot_grid(p1 + theme(legend.position = "none"), p2 + theme(legend.position = "none"), rel_widths = c(1, 1), labels = c('B', 'C'), label_size = 25)
legend_p <- cowplot::get_plot_component(p1  + theme_cowplot(font_size = 23) + theme(legend.box.margin = margin(-10, 0, 0, 0)) + 
                                          theme(
                                            legend.position = "bottom",
                                            legend.justification = "center",
                                            legend.direction = "horizontal",
                                            legend.box = "vertical",
                                          ) + guides(color = guide_legend(nrow = 2, byrow = TRUE)), 'guide-box-bottom', return_all = TRUE)
p <- plot_grid(p, legend_p, nrow = 2, rel_heights = c(4, 1))
p
c1 <- plot_grid(p, p3, ncol = 2, labels = c(' ', 'D'), rel_widths = c(1.85, 1), label_size = 25)
c1
taxa_counts <- read_tsv("../data/ref_taxa_counts.txt")
taxa_counts$Rank <- factor(
  taxa_counts$Rank,
  levels = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
)

p4 <- ggplot(taxa_counts %>% filter(Rank != "kingdom"),
             aes(x=Rank, y=Count)) +
  geom_violin(aes(color=Rank), linewidth=0.8, scale = "width", color="gray20", trim = TRUE, draw_quantiles = TRUE) +
  geom_point(aes(x=Rank, y=Count, color=Rank), alpha = 0.75, size = 1, position = "jitter") +
  stat_summary(color="black", fun = median) +
  scale_y_continuous(trans = "log2") +
  theme_cowplot(font_size = 23) +
  theme(strip.background = element_rect(fill = "gray")) +
  scale_colour_manual(values = palette.colors(palette = "Paired", n = 7)[-1]) +
  theme(
    panel.spacing.x = unit(1, "lines"),
    axis.title.x = element_blank() 
  ) +
  labs(y = "Number of\ngenomes in taxon", x = "")
p4
cf <- plot_grid(p4+ theme(legend.position = "none"), c1, ncol = 1, labels = c('A', ''), rel_heights = c(2.5, 5), label_size = 25)
cf
ggsave2("../figures/shared_kmers_analysis.pdf", width=18, height = 9)

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
  "discriminative k-mers"
scores$Method[scores$Method == "negative child taxon count ranking ~ k=32 w=35 h=12 b=16 l=2 mer-count const. (0.00)"] <-
  "children discrim."
scores$Method[scores$Method == "positive species count ranking ~ k=32 w=35 h=12 b=16 l=2 mer-count const. (0.00)"] <-
  "shared k-mers"
scores$Method[scores$Method == "positive child taxon count ranking ~ k=32 w=35 h=12 b=16 l=2 mer-count const. (0.00)"] <-
  "children common"
scores$Method[scores$Method == "weighted sum of species counts w.r.t. k-mer coverages ~ k=32 w=35 h=12 b=16 l=2 mer-count const. (0.00)"] <-
  "taxon coverage"
scores_p1 <- scores %>%
  filter(Method %in% c("random", "taxon coverage")) %>%
  # mutate(Distance_to_closest = cut(Distance_to_closest, include.lowest = TRUE, breaks = c(0, 0.01, 0.5, 0.1, 0.25))) %>%
  group_by(Method, Taxonomic_rank) %>%
  summarise(Precision = mean(Precision), Recall = mean(Recall), F1 = mean(F1))
scores_p1$expt = "Species-based ranking (R)"
p1 <- ggplot(
    scores %>% filter(Method %in% c("random", "taxon coverage")) %>% filter(Method %in% c("random", "taxon coverage")) %>% filter(!is.na(Taxonomic_rank)) %>% group_by(Method, Taxonomic_rank) %>%  summarise(Precision = mean(Precision), Recall = mean(Recall), F1 = mean(F1))
) +
  aes(x = Taxonomic_rank, y = F1, color = Method, shape = Method) +
  # geom_boxplot() +
  # geom_point(size = 5,alpha = 1) +
  geom_point(size = 5, alpha = 1, position = position_dodge2(width = 0.5, preserve = "single", padding = 0.2), )+
  # geom_smooth() + 
  # facet_wrap(c("expt")) +
  # facet_wrap(c("Taxonomic_rank")) +
  labs(shape = "Approach", colour = "Approach", x = "Rank", y = "F1") +
  # geom_line(aes(group = Distance_to_closest), color = "grey20") +
  scale_color_manual(values = c("#377eb8", "#ff7f00")) + 
  theme_minimal_grid(font_size = 23) + scale_shape_manual(values = c(15, 19, 17)) +
  theme(
    # axis.text.y = element_text(),
    # axis.text.x = element_text(angle = 20, vjust = 0.75),
  ) + 
  # scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5))
p1

ggplot(scores %>% filter(Method %in% c("random", "shared k-mers", "discriminative k-mers"))  %>% filter(!is.na(Reference_genome_in_phylum))) +
  aes(y = F1, x = Taxonomic_rank, color = Method, shape = Method) +
  # geom_point(position = position_dodge2(width = 0.5), size = 0.75, alpha = 0.25) +
  stat_summary(aes(group = Method, color = Method), geom="line", linewidth = 2, alpha = 0.5, position=position_dodge2(width = 0.5)) +
  stat_summary(size = 1, position=position_dodge2(width = 0.5)) + 
  scale_color_manual(values = c("#e7298a", "#1b9e77", "#e6ab02")) + 
  scale_shape_manual(values = c(17, 16, 15)) + 
  # scale_x_continuous(trans="log10") +
  labs(x="Rank", color="Approach", shape="Approach") +
  theme_cowplot(font_size = 19) +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_text(), axis.text.x = element_text())
