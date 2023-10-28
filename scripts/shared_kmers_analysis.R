require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)
require(scales)
require(reshape2)

taxon_dist <- read_csv("../data/dist_wrt_lastcommonrank.csv")

taxon_dist$rank <- factor(
  taxon_dist$rank,
  levels = c("root", "kingdom", "phylum", "class", "order", "family", "genus", "species")
)
# taxon_dist <- melt(taxon_dist, id = c("g1", "g2", "rank"), measure = c("JI", "d")) 

p1 <- ggplot(taxon_dist %>% filter(rank !="species" & !rank %in% c("root","kingdom")), aes(x=rank, y=d)) +
  geom_boxplot() +
  theme_cowplot(font_size = 18) +
  scale_color_brewer(palette = "Paired")+
  scale_y_continuous(labels = percent, name="Pairwise distance")+
  stat_summary(color="blue")+
  scale_x_discrete(name="",labels = c("Phylum\nClass","Class\nOrder","Order\nFamily","Family\nGenus","Genus\nSpecies")) +
  theme(
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_blank()
  )

p2 <- ggplot(taxon_dist %>% filter(rank !="species" & !rank %in% c("root","kingdom")), aes(x=rank, y=JI)) +
  geom_boxplot() +
  theme_cowplot(font_size = 18) +
  scale_color_brewer(palette = "Paired")+
  scale_y_continuous(labels = percent, name="Shared 30-mers (portion)", trans = "log2")+
  stat_summary(color="blue")+
  scale_x_discrete(name="",labels = c("Phylum\nClass","Class\nOrder","Order\nFamily","Family\nGenus","Genus\nSpecies")) +
  theme(
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_blank()
    )

L=10^6
n=32
k=30
prob_change = function(k,d) 1-(1-d)^k
shared_kmers = function(k,d,n=n,L=L)  L*n - L*n*(prob_change(k,d)+((1-prob_change(k,d))*prob_change(k,d)^(n-1)))
shared_kmer_prop = function(k,d,n=n,L=L) shared_kmers(k,d,n=n,L=L)/(n*L)

p3 <- ggplot(aes(x=x),data=data.frame(x=1:1000))+
  stat_function(aes(color="5%"),fun=function(x) shared_kmer_prop(30,0.025,x+1,10^6))+
  stat_function(aes(color="10%"),fun=function(x) shared_kmer_prop(30,0.05,x+1,10^6))+
  stat_function(aes(color="15%"),fun=function(x) shared_kmer_prop(30,0.075,x+1,10^6))+
  stat_function(aes(color="20%"),fun=function(x) shared_kmer_prop(30,0.1,x+1,10^6))+
  stat_function(aes(color="25%"),fun=function(x) shared_kmer_prop(30,0.125,x+1,10^6))+
  stat_function(aes(color="33%"),fun=function(x) shared_kmer_prop(30,1/6,x+1,10^6))+
  theme_cowplot(font_size = 18) +
  labs(y="Expected shared 30-mers (portion)", x="Number of reference genomes", color="Within group\ndiversity") +
  scale_x_continuous(trans = "log10")+
  scale_y_continuous(labels = percent,trans = "log10") +
  theme(
    axis.text.y = element_text(size = 15),
    axis.text.x = element_text(size = 15),
    legend.position=c(0.7,0.25)
  ) + scale_color_discrete(breaks=c("5%", "10%", "15%", "20%", "25%", "33%"), type = palette.colors(n=6, "Paired"))

plot_grid(p1, p2, p3, ncol = 3, labels = c('A', 'B', 'C', 'D'), rel_widths = c(1, 1, 1.25), label_size = 19)
ggsave2("../figures/shared_kmers_analysis.pdf", width=18, height = 6)