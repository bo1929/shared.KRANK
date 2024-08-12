require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)
require(tidyr)
require(tidyverse)

mm <- read_tsv("../results/resultsCAMI2-marine.tsv")
sm <- read_tsv("../results/resultsCAMI2-strain_madness.tsv")
mm$Dataset <- "Marine dataset"
sm$Dataset <- "Strain-madness dataset"
pm <- sm # rbind(sm, mm)
pm <- pm %>% filter(rank != "strain")
pm <- pm %>% group_by(sample, metric) %>% mutate(p = value)
pm <- pm %>% filter(tool != "Gold standard")
pm  <- pm %>% filter(
  tool %in% c(
    "KRANK v0.3.2",
    "Bracken 2.2",
    "CCMetagen 1.1.3",
    "Centrifuge 1.0.4 beta",
    "CONSULT-II v0.4.0",
    "DUDes 0.08",
    "FOCUS 1.5",
    # "Metalign 0.6.2",
    "MetaPhlAn 2.9.22",
    # "MetaPhyler 1.25",
    "mOTUs 2.5.1_2",
    # "MetaPalette 1.0.0",
    "TIPP 4.3.10"
    # "NBC++"
  )
)
pm$tool[pm$tool == "mOTUs 2.5.1_2"] <- "mOTUs 2.5.1"
pm$rank <- factor(
  pm$rank,
  levels = c("superkingdom", "phylum", "class", "order", "family", "genus", "species")
)
pm$metric <- as.factor(pm$metric)
pm$tool <- factor(pm$tool, levels = c(
    "KRANK v0.3.2",
    "CONSULT-II v0.4.0",
    "Bracken 2.2",
    "CCMetagen 1.1.3",
    "Centrifuge 1.0.4 beta",
    "DUDes 0.08",
    "FOCUS 1.5",
    # "Metalign 0.6.2",
    "MetaPhlAn 2.9.22",
    # "MetaPhyler 1.25",
    "mOTUs 2.5.1",
    "TIPP 4.3.10"
    # "NBC++",
    # "MetaPalette 1.0.0"
  )
)

pm <- rbind(pm, pm %>% filter(metric == "Completeness") %>% mutate(p = 1 - p, metric = "1 - Completeness"))
pm <- rbind(pm, pm %>% filter(metric == "Purity") %>% mutate(p = 1 - p, metric = "1 - Purity"))
pm <- rbind(pm, pm %>% filter(metric == "Purity-f") %>% mutate(p = 1 - p, metric = "1 - Purity-f"))

pm <- pm %>%
  group_by(rank, metric, Dataset, sample) %>%
  mutate(r = dense_rank(p))

pm <- pm %>%
  group_by(tool, metric, rank, Dataset) %>%
  summarise(r = sum(r), p = mean(p))

plot_rankm <- function(xyz){
  pm %>% 
    filter(metric == xyz) %>%
    mutate(toolx = tidytext::reorder_within(tool, p, within = Dataset)) %>% 
    ggplot() + aes(x=p, y=reorder(toolx, p), color=tool, label=tool) + 
    geom_point() + geom_text(hjust=-0.1, vjust=-0.1) + 
    facet_wrap(facets = vars(Dataset), scales = "free") +
    labs(x = xyz, y = "Tool rank") +
    tidytext::scale_y_reordered() +
    scale_colour_manual(values=c("black", "darkgray", RColorBrewer::brewer.pal(12,"Paired"))) +
    theme_cowplot() + 
    theme(legend.position = "none", axis.text.y = element_blank())
}

plot_rankm("Unweighted UniFrac (CAMI)")
plot_rankm("Weighted UniFrac (CAMI)")

plot_ranksva <- function(gmetric){
  pm %>% 
    filter(rank != "superkingdom" & metric == gmetric) %>%
    group_by(tool, Dataset, rank == "species") %>% summarize(p =  mean(p, na.rm = T))  %>%
    mutate(toolx = tidytext::reorder_within(tool, p, within =`rank == "species"`, sep = "_-_")) %>%
    mutate(toolx = tidytext::reorder_within(toolx, p, within = Dataset)) %>%
    ggplot() + aes(x=p, y=reorder(toolx, p), color=tool, label=tool) + 
    geom_point() + geom_text(hjust=-0.1, vjust=-0.1) + 
    labs(x = gmetric, y = "Tool rank") +
    scale_colour_manual(values=c("black", "darkgray", RColorBrewer::brewer.pal(12,"Paired"))) +
    theme_cowplot() + 
    tidytext::scale_y_reordered() +
    # facet_wrap(`rank == "species"` ~ Dataset, scales = "free", labeller = labeller(`rank == "species"` = function(x){print(x);x[c(3,4)]="Higher ranks";x[c(1,2)]="Species";x})) +
    theme(legend.position = "none", axis.text.y = element_blank()) +
    theme(panel.spacing = unit(.05, "lines"), panel.border = element_rect(color = "darkgray", fill = NA, size = 1), strip.background = element_rect(color = "darkgray", size = 1))
}

plot_ranksva("L1 norm error")
plot_ranksva("1 - Completeness")
plot_ranksva("1 - Purity")

plot_rankstog <- function(gmetric){
  pm %>% 
    filter(rank != "superkingdom" & metric == gmetric) %>%
    group_by(tool, Dataset, rank) %>% summarize(p =  mean(p, na.rm = T))  %>%
    mutate(toolx = tidytext::reorder_within(tool, p, within = list(rank, Dataset))) %>%
    ggplot() + aes(x=p, y=reorder(toolx, p), color=tool, label=tool) + 
    geom_point() + geom_text(hjust=-0.1, vjust=-0.1) + 
    tidytext::scale_y_reordered() +
    labs(x = gmetric, y = "Tool rank") +
    facet_wrap(Dataset ~ rank, scales = "free", nrow = 2) +
    scale_colour_manual(values=c("black", "darkgray", RColorBrewer::brewer.pal(12,"Paired"))) +
    theme_cowplot() + 
    theme(legend.position = "none", axis.text.y = element_blank()) +
    theme(panel.spacing = unit(.05, "lines"), panel.border = element_rect(color = "darkgray", fill = NA, size = 1), strip.background = element_rect(color = "darkgray", size = 1))
}

plot_rankstog("L1 norm error")
plot_rankstog("1 - Completeness")
plot_rankstog("1 - Purity")

xm <- pm %>% # filter(rank != "species" | is.na(rank)) %>% 
  group_by(tool, metric, Dataset) %>%
  summarise(r = sum(r), p = mean(p))

xm$metric[xm$metric == "1 - Purity-f"] <- "1 - Purity (1% filtered)"

xm %>% filter(metric %in% c("1 - Completeness", "1 - Purity-f", "L1 norm error", "Weighted UniFrac (CAMI)")) %>% 
  mutate(toolx = tidytext::reorder_within(tool, p, within = list(metric, Dataset))) %>%
  ggplot() + aes(x=r, y=reorder(toolx, r), color=tool, label=tool) + 
  geom_point() + geom_text(hjust=-0.1, vjust=-0.1) + 
  facet_wrap(facets = vars(Dataset), scales = "free") +
  # labs(x = xyz, y = "Tool rank") +
  tidytext::scale_y_reordered() +
  facet_wrap(Dataset ~ metric, scale = "free") +
  theme_cowplot() + 
  scale_colour_manual(values=c("black", "darkgray", RColorBrewer::brewer.pal(12,"Paired"))) +
  theme(panel.spacing = unit(.05, "lines"), panel.border = element_rect(color = "darkgray", fill = NA, size = 1), strip.background = element_rect(color = "darkgray", size = 1)) + 
  theme(legend.position = "none", axis.text.y = element_blank())

xm %>% filter(metric %in% c("1 - Completeness", "1 - Purity", "L1 norm error", "Weighted UniFrac (CAMI)")) %>% 
  mutate(toolx = tidytext::reorder_within(tool, p, within = list( Dataset))) %>% 
  group_by(Dataset,metric) %>%
  mutate(r = dense_rank(r)) %>%
  filter(grepl("Strain",Dataset)) %>%
  mutate(toolx = sub("___.*","",toolx)) %>%
  group_by(metric,toolx) %>%
  mutate(r=mean(r)) %>%
  group_by(toolx) %>%
  mutate(rt=mean(r)) %>%
  ggplot(aes(fill=r,x=metric,y=reorder(toolx,r),label=r)) + 
  geom_tile(aes(x=" mean",fill=rt))+
  geom_text(aes(x=" mean",label=round(rt,1)))+
  geom_tile()+
  #tidytext::scale_y_reordered()+
  #facet_wrap(~Dataset,scales = "free")+
  geom_text()+
  scale_fill_viridis_b(direction = -1 ,breaks = (0:10)+0.000000000000001,name="")+
  xlab("")+
  ylab("")+
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle=90, hjust=0.95, vjust=0.2))
# ggsave2("../../appeal-files/cami2_summary-strain_madness.pdf", width = 9, height = 8)

ggplot() + aes(x=r, y=reorder(toolx, r), color=tool, label=tool) + 
  geom_point() + geom_text(hjust=-0.1, vjust=-0.1) + 
  facet_wrap(facets = vars(Dataset), scales = "free") +
  # labs(x = xyz, y = "Tool rank") +
  tidytext::scale_y_reordered() +
  facet_wrap(Dataset ~ metric, scale = "free") +
  theme_cowplot() + 
  scale_colour_manual(values=c("black", "darkgray", RColorBrewer::brewer.pal(12,"Paired"))) +
  theme(panel.spacing = unit(.05, "lines"), panel.border = element_rect(color = "darkgray", fill = NA, size = 1), strip.background = element_rect(color = "darkgray", size = 1)) + 
  theme(legend.position = "none", axis.text.y = element_blank())

xm %>% filter(metric %in% c("1 - Completeness", "L1 norm error", "Weighted UniFrac (CAMI)")) %>% 
  group_by(tool, Dataset) %>%
  summarise(r = sum(r), p = mean(p)) %>% 
  mutate(toolx = tidytext::reorder_within(tool, p, within = list(Dataset))) %>%
  ggplot() + aes(x=r, y=reorder(toolx, r), color=tool, label=tool) + 
  geom_point() + geom_text(hjust=-0.1, vjust=-0.1) + 
  facet_wrap(facets = vars(Dataset), scales = "free") +
  # labs(x = xyz, y = "Tool rank") +
  tidytext::scale_y_reordered() +
  facet_wrap(vars(Dataset), scale = "free") +
  theme_cowplot() + 
  scale_colour_manual(values=c("black", "darkgray", RColorBrewer::brewer.pal(12,"Paired"))) +
  theme(panel.spacing = unit(.05, "lines"), panel.border = element_rect(color = "darkgray", fill = NA, size = 1), strip.background = element_rect(color = "darkgray", size = 1)) + 
  theme(legend.position = "none", axis.text.y = element_blank())