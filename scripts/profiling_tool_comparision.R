require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)

pm <- read_tsv("../results/all_tools-profiling_evaluation-CAMI1_hc.tsv")
pm <- pm %>% filter(metric == "Bray-Curtis distance" | metric == "Shannon equitability")
pm <- pm %>% filter(rank != "kingdom")  %>% filter(rank != "strain")
pm <- pm %>%
  group_by(sample, metric) %>%
  mutate(dvalue = abs(value - value[tool == "Gold standard"]))
pm <- pm %>% filter(tool != "Gold standard")
pm$rank <- factor(
  pm$rank,
  levels = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
)
pm <- pm %>%
  group_by(rank, tool, metric) %>%
  summarise(lower = min(dvalue), upper = max(dvalue), p = mean(dvalue))

bc_plot <- ggplot(pm %>% filter(metric == "Bray-Curtis distance"), mapping = aes(x = p, y = tool, color = tool)) +
  geom_pointrange(size = 0.7, linewidth = 1, mapping = aes(xmin = lower, xmax = upper)) +
  facet_wrap(vars(rank), scales = "free_x") +
  xlim(0, NA) +
  theme_minimal_vgrid(font_size = 17) +
  scale_colour_manual(values=c( "#e31a1c", "#ff7f00", "#33a02c", "#1f78b4")) +
  labs(color = "Tool", y = "", x = "Bray-Curtis dissimilarity to true profile") +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(vjust = 0.5, hjust = 1)) +
  theme(panel.spacing.x = unit(1.5, "lines"), panel.spacing.y = unit(1.5, "lines"))

se_plot <- ggplot(pm %>% filter(metric == "Shannon equitability"), mapping = aes(x = p, y = tool, color = tool)) +
  geom_pointrange(size =  0.7, linewidth = 1, mapping = aes(xmin = lower, xmax = upper)) +
  facet_wrap(vars(rank), scales = "free_x") +
  xlim(0, NA) +
  theme_minimal_vgrid(font_size = 17) +
  scale_colour_manual(values=c( "#e31a1c", "#ff7f00", "#33a02c", "#1f78b4")) +
  labs(color = "Tool", y = "", x = "Shannon's equitability (absolute difference to gold standard)") +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(vjust = 0.5, hjust = 1)) +
  theme(panel.spacing.x = unit(1.5, "lines"), panel.spacing.y = unit(1.5, "lines"))

prow <- plot_grid(
  bc_plot + theme(legend.position = "none"),
  se_plot + theme(legend.position = "none"),
  ncol = 1, vjust = 2,
  rel_heights = c(1, 1), labels = c("A", "B"), label_size = 19
)
legend <- get_legend(
  # create some space to the left of the legend
  se_plot + theme(legend.position = "bottom", legend.justification = "center", legend.margin = margin(-30,0,0,0))
)
plot_grid(prow, legend, nrow=2, rel_heights = c(8, 1))

ggsave2("../figures/profiling_tool_comparison.pdf", width = 13, height = 6.5)
