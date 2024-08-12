require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)
library(tidyr)
library(stringr)

dropLeadingZero <- function(l) {
  str_replace(l, '.00', '')
}
drop_version <- function(l) {
  str_replace(l, " .*", "")
}

mm <- read_tsv("../results/resultsCAMI2-marine.tsv")
sm <- read_tsv("../results/resultsCAMI2-strain_madness.tsv")
mm$Dataset <- "Marine dataset"
sm$Dataset <- "Strain-madness dataset"
pm <- rbind(sm, mm)
pm <- pm %>% filter(rank != "strain")
pm <- pm %>% filter(tool != "Gold standard")
pm <- pm %>%
  group_by(rank, tool, metric, Dataset) %>%
  summarise(lower = min(value),
            upper = max(value),
            p = mean(value))
pm  <- pm %>% filter(
  tool %in% c(
    "KRANK v0.3.2",
    "Bracken 2.2",
    "CCMetagen 1.1.3",
    "Centrifuge 1.0.4 beta",
    "CONSULT-II v0.4.0",
    "DUDes 0.08",
    "FOCUS 1.5",
    "Metalign 0.6.2",
    "MetaPhlAn 2.9.22",
    "MetaPhyler 1.25",
    "mOTUs 2.5.1_2",
    "MetaPalette 1.0.0",
    "TIPP 4.3.10",
    "NBC++"
  )
)
pm$tool[pm$tool == "mOTUs 2.5.1_2"] <- "mOTUs 2.5.1"
pm$rank <- factor(
  pm$rank,
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
pm$metric <- as.factor(pm$metric)
pm$tool <- factor(
  pm$tool,
  levels = c(
    "KRANK v0.3.2",
    "CONSULT-II v0.4.0",
    "Bracken 2.2",
    "CCMetagen 1.1.3",
    "Centrifuge 1.0.4 beta",
    "DUDes 0.08",
    "FOCUS 1.5",
    "Metalign 0.6.2",
    "MetaPhlAn 2.9.22",
    "MetaPhyler 1.25",
    "mOTUs 2.5.1",
    "TIPP 4.3.10",
    "NBC++",
    "MetaPalette 1.0.0"
  )
)

bc_plot <- ggplot(
  pm %>% filter(metric == "Bray-Curtis distance" &
                  rank != "superkingdom"),
  mapping = aes(
    x = p,
    y = reorder(tool, p),
    color = tool,
    shape = grepl("CONS", tool)
  )
) +
  geom_pointrange(
    size = 0.25,
    linewidth = 0.5,
    mapping = aes(xmin = lower, xmax = upper)
  ) +
  facet_wrap(vars(rank), scales = "free_x") +
  xlim(0, NA) +
  theme_minimal_vgrid(font_size = 17) +
  labs(color = "Tool", y = "", x = "Bray-Curtis dissimilarity to true profile") +
  scale_colour_manual(values = c("black", "gray", RColorBrewer::brewer.pal(12, "Paired"))) +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(
          vjust = 0.5,
          hjust = 1,
          size = 9
        )) +
  theme(aspect.ratio = 0.4,
        panel.spacing.x = unit(1.5, "lines")) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_shape_manual(values = c(1, 16), guide = "none")
bc_plot
se_plot <- ggplot(
  pm %>% filter(metric == "Shannon equitability" &
                  rank != "superkingdom"),
  mapping = aes(
    x = p,
    y = tool,
    color = tool,
    shape = grepl("CONS", tool)
  )
) +
  geom_pointrange(
    size = 0.25,
    linewidth = 0.5,
    mapping = aes(xmin = lower, xmax = upper)
  ) +  facet_wrap(vars(rank), scales = "free_x") +
  xlim(0, NA) +
  theme_minimal_vgrid(font_size = 17) +
  labs(color = "Tool", y = "", x = "Shannon's equitability (absolute difference to gold standard)") +
  scale_colour_manual(values = c("black", "gray", RColorBrewer::brewer.pal(12, "Paired"))) +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(
          vjust = 0.5,
          hjust = 1,
          size = 9
        )) +
  theme(aspect.ratio = 0.4,
        panel.spacing.x = unit(1.5, "lines")) +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_shape_manual(values = c(1, 16), guide = "none")
se_plot

prow <- plot_grid(
  bc_plot + theme(legend.position = "none"),
  se_plot + theme(legend.position = "none"),
  ncol = 2,
  rel_heights = c(1, 1),
  vjust = 5,
  hjust = 5
)
legend <- get_legend(# create some space to the left of the legend
  se_plot + theme(legend.box.margin = margin(0, 0, 0, 12)))
plot_grid(prow, legend, rel_widths = c(3, .45))

tm <- pivot_wider(
  pm,
  id_cols = c(rank, tool, Dataset),
  names_from = metric,
  values_from = p
)
p1 <- ggplot(
  tm %>% filter(rank != "na" & rank != "superkingdom"),
  aes(
    y = `Completeness`,
    x = `Purity`,
    color = tool,
    shape = grepl("CONS|KRANK", tool)
  )
) +
  labs(
    x = "Purity",
    y = "Completeness",
    shape = "",
    color = "Tool"
  ) +
  geom_point(size = 3, alpha = 0.9) +
  scale_shape_discrete(guide = "none") +
  theme_cowplot(font_size = 20) +
  scale_colour_manual(
    values = c(
      "#1F78B4",
      "#33A02C",
      "#E31A1C",
      "salmon",
      "darkgray",
      "#F8CEE3",
      "#B2DF8A",
      "#FB9A99",
      "#FDBF6F",
      "#FF7F00",
      "#CAB2D6",
      "#6A3D9A",
      "#FFFF99",
      "#B15928"
    )
  ) +
  scale_x_continuous(labels = dropLeadingZero) +
  scale_y_continuous(labels = dropLeadingZero) +
  facet_grid(rows = vars(Dataset),
             cols = vars(rank),
             scales = "free_x") +
  theme(
    panel.spacing = unit(.05, "lines"),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      size = 1
    ),
    strip.background = element_rect(color = "black", size = 1),
    legend.position = "bottom"
  )
p1
ggsave2(
  "../figures/profiling-cami2_combined-tool_comparison-completeness_purity.pdf",
  width = 12,
  height = 7
)

tm$`Unweighted UniFrac (CAMI)` <- rep((tm %>% filter(is.na(rank)))$`Unweighted UniFrac (CAMI)`, 8)
tm$`Weighted UniFrac (CAMI)` <- rep((tm %>% filter(is.na(rank)))$`Weighted UniFrac (CAMI)`, 8)
equal_breaks <- function(n = 3, s = 0.1, ...) {
  function(x) {
    d <- s * diff(range(x)) / (1 + 2 * s)
    seq = seq(min(x) + d, max(x) - d, length = n)
    round(seq, -floor(log10(abs(seq[2] - seq[1]))))
  }
}

tm$Dataset <- factor(tm$Dataset,
                     levels = c("Strain-madness dataset", "Marine dataset"))

p2 <- tm %>% group_by(tool, rank) %>% summarize(
  `Weighted UniFrac (CAMI)` = mean(`Weighted UniFrac (CAMI)`),
  `L1 norm error` = mean(`L1 norm error`, na.rm = T)
)  %>%
  ggplot(aes(
    x = `L1 norm error`,
    y = reorder(tool, -`Weighted UniFrac (CAMI)`),
    color = tool,
    fill = tool,
    shape = grepl("CONS|KRANK", tool)
  )) +
  geom_point(size = 3, alpha = 0.9) +
  labs(
    x = "Purity",
    y = "Completeness",
    shape = "",
    color = "Tool"
  ) +
  scale_shape_discrete(guide = "none") +
  facet_wrap(vars(rank)) +
  scale_color_manual(
    values = c(
      "#1F78B4",
      "#33A02C",
      "#E31A1C",
      "salmon",
      "darkgray",
      "#F8CEE3",
      "#B2DF8A",
      "#FB9A99",
      "#FDBF6F",
      "#FF7F00",
      "#CAB2D6",
      "#6A3D9A",
      "#FFFF99",
      "#B15928"
    )
  ) +
  labs(
    title = "Strain-madness dataset of CAMI-II",
    x = "L1 norm error",
    y = "",
    fill = "",
    shape = "",
    color = ""
  ) +
  ylim(8, NA) +
  xlim(0, 2) +
  scale_shape_discrete(guide = "none") +
  theme_cowplot(font_size = 20) +
  scale_fill_manual(
    values = c(
      "#1F78B4",
      "#33A02C",
      "#E31A1C",
      "salmon",
      "darkgray",
      "#F8CEE3",
      "#B2DF8A",
      "#FB9A99",
      "#FDBF6F",
      "#FF7F00",
      "#CAB2D6",
      "#6A3D9A",
      "#FFFF99",
      "#B15928"
    )
  ) +
  scale_x_continuous(
    labels = dropLeadingZero,
    limits = c(0, 1.0),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0)
  ) +
  scale_y_discrete(labels = drop_version)  +
  theme(
    legend.text = element_text(size = 14),
    panel.spacing = unit(1.0, "lines"),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      size = 1
    ),
    strip.background = element_rect(color = "black", size = 1),
    legend.position = "none",
    legend.box.margin = margin(0, 0, 0, -40)
  ) +
  guides(fill = guide_legend(nrow = 4), shape = "none")
p2
ggsave2("../figures/profiling-cami2_combined-tool_comparison-l1_unifrac-avg_all.pdf", width = 8, height = 9)
tm$isSpecies <- tm$rank == "species"
tm$isSpecies[tm$isSpecies] = "species"
tm$isSpecies[tm$isSpecies == FALSE] = "higher ranks"

p2 <- tm %>% filter(!is.na(isSpecies)) %>% filter(Dataset == "Strain-madness dataset") %>%
  group_by(tool, isSpecies, Dataset) %>% summarize(
    `Weighted UniFrac (CAMI)` = mean(`Weighted UniFrac (CAMI)`, na.rm = T),
    `L1 norm error` = mean(`L1 norm error`, na.rm = T)
  )  %>%
  ggplot(aes(
    x = `L1 norm error`,
    y = `Weighted UniFrac (CAMI)`,
    color = tool,
    shape = grepl("CONS|KRANK", tool)
  )) +
  geom_point(size = 4, alpha = 0.9) +
  facet_wrap(~`rank == "species"`, labeller = function(x) {x[1,1]="Higher ranks";x[2,1]="Species";x})+
  facet_grid(isSpecies ~ Dataset) +
  facet_wrap(facets = vars(isSpecies), scale = "free_x") +
  scale_colour_manual(
    values = c(
      "#1F78B4",
      "#33A02C",
      "#E31A1C",
      "salmon",
      "darkgray",
      "#F8CEE3",
      "#B2DF8A",
      "#FB9A99",
      "#FDBF6F",
      "#FF7F00",
      "#CAB2D6",
      "#6A3D9A",
      "#FFFF99",
      "#B15928"
    )
  ) +
  labs(title="Strain-madness dataset of CAMI-II", x="L1 norm error", y="weighted UniFrac error", shape="", color="Tool") +
  ylim(8, NA) +
  xlim(1, 1.95) +
  scale_shape_discrete(guide = "none") +
  theme_cowplot(font_size = 15) +
  scale_colour_manual(
    values = c(
      "#1F78B4",
      "#33A02C",
      "#E31A1C",
      "salmon",
      "darkgray",
      "#F8CEE3",
      "#B2DF8A",
      "#FB9A99",
      "#FDBF6F",
      "#FF7F00",
      "#CAB2D6",
      "#6A3D9A",
      "#FFFF99",
      "#B15928"
    ),
    labels = drop_version
  ) +
  scale_x_continuous(breaks = equal_breaks(n = 3, s = 0.25)) +
  scale_y_continuous(labels = dropLeadingZero) +
  theme(
    panel.spacing = unit(.01, "lines"),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      size = 1
    ),
    strip.background = element_rect(color = "black", size = 1),
    legend.position = "right"
  )
p2
ggsave2(
  "../figures/profiling-cami2_combined-tool_comparison-l1_unifrac.pdf",
  width = 12,
  height = 7
)

p3 <- tm %>%  filter(!is.na(rank) & rank != "superkingdom") %>%
  # group_by(tool,rank=="species") %>% summarize( `Weighted UniFrac (CAMI)`=mean(`Weighted UniFrac (CAMI)`),`L1 norm error`= mean(`L1 norm error`,na.rm = T))  %>%
  ggplot(aes(
    x = 2 - `L1 norm error`,
    y = 16 - `Weighted UniFrac (CAMI)`,
    color = tool,
    shape = grepl("CONS|KRANK", tool)
  )) +
  geom_point(size = 5, alpha = 0.9) +
  scale_shape_discrete(guide = "none") +
  theme_cowplot(font_size = 20) +
  scale_colour_manual(
    values = c(
      "#1F78B4",
      "#33A02C",
      "#E31A1C",
      "salmon",
      "darkgray",
      "#F8CEE3",
      "#B2DF8A",
      "#FB9A99",
      "#FDBF6F",
      "#FF7F00",
      "#CAB2D6",
      "#6A3D9A",
      "#FFFF99",
      "#B15928"
    )
  ) +
  scale_x_continuous(labels = dropLeadingZero) +
  scale_y_continuous(labels = dropLeadingZero) +
  facet_grid(rows = vars(Dataset),
             cols = vars(rank),
             scales = "free_x") +
  theme(
    legend.title = element_blank(),
    panel.spacing = unit(.05, "lines"),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      size = 1
    ),
    strip.background = element_rect(color = "black", size = 1),
    legend.position = "bottom"
  ) + labs(
    x = "2 - L1 norm error",
    y = "16 - weighted UniFrac error",
    shape = "",
    color = ""
  ) +
  ylim(8, NA) +
  xlim(1, NA) +
  scale_x_continuous(breaks = equal_breaks(n = 3, s = 0.25)) +
  # scale_x_continuous(labels = dropLeadingZero) +
  # scale_y_continuous(labels = dropLeadingZero)
  scale_shape_discrete(guide = "none") +
  guides(colour = guide_legend(nrow = 3), shape = "none")
p3
ggsave2(
  "../figures/profiling-cami2_combined-tool_comparison-l1_unifrac.pdf",
  width = 13,
  height = 8
)
plot_grid(
  rel_heights = c(8, 1),
  plot_grid(
    p3 + theme(legend.position = "none"),
    p2 + theme(legend.position = "none"),
    ncol = 2
  ),
  get_legend(
    p2 + theme(
      legend.box.margin = margin(0, 0, 0, 0),
      legend.position = "bottom",
      legend.justification = "center"
    ) +  guides(colour = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))
  ),
  nrow = 2
)
ggsave2(
  "../figures/profiling-cami2_combined-tool_comparison.pdf",
  width = 16,
  height = 4.5
)
