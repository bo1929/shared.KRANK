require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)
require(reshape2)
require(tidyr)

pma <- read_tsv("../results/all_tools-profiling_evaluation-CAMI1_hc.tsv")
pma <- pma %>% filter(metric == "Bray-Curtis distance" | metric == "Shannon equitability")
pma <- pma %>% filter(rank != "kingdom")  %>% filter(rank != "strain") %>% filter(rank != "superkingdom")
pma <- pma %>%
  group_by(sample, metric) %>%
  mutate(dvalue = abs(value - value[tool == "Gold standard"]))
pma <- pma %>% filter(tool != "Gold standard")
pma$rank <- factor(
  pma$rank,
  levels = c("kingdom", "phylum", "class", "order", "family", "genus", "species")
)
pma <- pma %>%
  group_by(rank, tool, metric) %>%
  summarise(lower = min(dvalue), upper = max(dvalue), p = mean(dvalue))
pm <- pma %>% filter(!tool %in% c("KRANK-hs (Eq. 7, 51.2Gb)", "KRANK-hs (Eq. 6, 51.2Gb)", "CONSULT-II (Eq. 6, 140.7Gb)"))
pm$tool[pm$tool == "Bracken (46.5Gb)"] <- "Bracken v2.8 (46.5GB)"
pm$tool[pm$tool == "CLARK (149.6Gb)"] <- "CLARK v1.2.6.1 (149.6GB)"
# pm$tool[pm$tool == "KRANK-hs (Eq. 7, 51.2Gb)"] <- "KRANK-hs (51.2GB)" # v0.4.0
pm$tool[pm$tool == "KRANK-hs (v0.5.0, 51.2Gb)"] <- "KRANK-hs v0.3.2 (51.2GB)" # v0.5.0
pm$tool[pm$tool == "CONSULT-II (Eq. 7, 140.7Gb)"] <- "CONSULT-II v0.4.0 (140.7GB)"

bc_plot <- ggplot(pm %>% filter(metric == "Bray-Curtis distance"), mapping = aes(x = p, y = tool, color = tool)) +
  geom_pointrange(size = 1.2, linewidth = 1, mapping = aes(xmin = lower, xmax = upper)) +
  facet_wrap(vars(rank), scales = "free_x") +
  xlim(0, NA) +
  theme_minimal_vgrid(font_size = 23) +
  scale_colour_manual(values=c("#e31a1c", "#ff7f00" , "#33a02c", # "#ff1a1c",  
                               "#1f78b4" #  "#4f78b4"
                               )) + labs(color = "Tool", y = "", x = "Bray-Curtis dissimilarity to true profile") +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(vjust = 0.5, hjust = 1)) +
  theme(panel.spacing.x = unit(1.5, "lines"), panel.spacing.y = unit(1.5, "lines"))

se_plot <- ggplot(pm %>% filter(metric == "Shannon equitability"), mapping = aes(x = p, y = tool, color = tool)) +
  geom_pointrange(size =  1.2, linewidth = 1, mapping = aes(xmin = lower, xmax = upper)) +
  facet_wrap(vars(rank), scales = "free_x") +
  xlim(0, NA) +
  theme_minimal_vgrid(font_size = 23) +
  scale_colour_manual(values=c("#e31a1c", "#ff7f00" , "#33a02c", # "#ff1a1c",  
                               "#1f78b4" #  "#4f78b4"
  )) +  labs(color = "", y = "", x = "Shannon's equitability (absolute difference to gold standard)") +
  theme(strip.background = element_rect(fill = "gray")) +
  theme(axis.text.y = element_blank(), axis.text.x = element_text(vjust = 0.5, hjust = 1)) +
  theme(panel.spacing.x = unit(1.35, "lines"), panel.spacing.y = unit(1.35, "lines"))

prow <- plot_grid(
  bc_plot + theme(legend.position = "none"),
  se_plot + theme(legend.position = "bottom", legend.box.margin = margin(0,0,0,-45)),
  ncol = 1, vjust = 2,
  rel_heights = c(1, 1.2), labels = c("A", "B"), label_size = 25
)
prow
ggsave2("../figures/profiling_tool_comparison.pdf", width = 15, height = 9)

pm <- pma %>% filter(tool !=  "KRANK-hs (v0.5.0, 51.2Gb)") # %>% filter(tool %in% c("KRANK-hs (old)", "CONSULT-II (new)", "KRANK-hs (new)", "CONSULT-II (old)"))
pm$tool[pm$tool == "KRANK-hs (Eq. 6, 51.2Gb)"] <- "KRANK-hs v0.3.2 + CONSULT-II v0.4.0 (with v0.3.0 - 51.2Gb)"
pm$tool[pm$tool == "CONSULT-II (Eq. 6, 140.7Gb)"] <- "CONSULT-II v0.4.0 (with v0.3.0 - 140.7Gb)"
pm$tool[pm$tool == "KRANK-hs (Eq. 7, 51.2Gb)"] <- "KRANK-hs v0.3.2 + CONSULT-II v0.4.0 (with v0.4.0 - 51.2Gb)"
pm$tool[pm$tool == "CONSULT-II (Eq. 7, 140.7Gb)"] <- "CONSULT-II v0.4.0 (with v0.4.0 - 140.7Gb)"
pm$tool[pm$tool == "Bracken (46.5Gb)"] <- "Bracken v2.8 (46.5Gb)"
pm$tool[pm$tool == "CLARK (149.6Gb)"] <- "CLARK v1.2.6.1 (149.6Gb)"
pm %>% mutate(old = grepl("with v0.4.0", tool),tool = sub("[(].*[)]","", tool))  %>% 
  select(rank, tool, metric, p, old) %>%
  pivot_wider(values_from = p,names_from = c(old,metric) ) %>%
ggplot(
       mapping = aes(x = `FALSE_Bray-Curtis distance`, xend = `TRUE_Bray-Curtis distance`, 
                     y =`FALSE_Shannon equitability` , yend = `TRUE_Shannon equitability`, color=tool)) +
  geom_segment(size = 1, linewidth = 0.75, arrow = arrow(length = unit(0.35, "cm")), show.legend = FALSE) +
  geom_point(size=2) +
  scale_colour_manual(values=c("#e31a1c", "#ff7f00" , "#33a02c", "#1f78b4" ))+ xlim(0, NA) + ylim(0, NA) +
  labs(title = "CONSULT-II profiling v0.3.0 versus v0.4.0 with genome size correction", color = "Tool", x = "Bray-Curtis dissimilarity to true profile", y = "Shannon's equitability (absolute diff. to gold standard)") +
  facet_wrap(vars(rank), scales = "free") + theme_minimal_grid(font_size = 17) +   theme(panel.spacing.x = unit(1.5, "lines"), panel.spacing.y = unit(1.5, "lines"), legend.position = "bottom")
ggsave2("../figures/improvement_profiling-genome_size_correction.pdf", width = 14, height = 7)

pm <- pma %>% filter(!tool %in% c("KRANK-hs (Eq. 6, 51.2Gb)", "CONSULT-II (Eq. 6, 140.7Gb)"))
pm$tool[pm$tool == "KRANK-hs (Eq. 7, 51.2Gb)"] <- "KRANK-hs v0.3.2 + CONSULT-II v0.5.0 (v0.4.0 - 51.2Gb)"
pm$tool[pm$tool == "Bracken (46.5Gb)"] <- "Bracken v2.8 (46.5Gb)"
pm$tool[pm$tool == "CLARK (149.6Gb)"] <- "CLARK v1.2.6.1 (149.6Gb)"
pm$tool[pm$tool == "KRANK-hs (v0.5.0, 51.2Gb)"] <- "KRANK-hs v0.3.2 + CONSULT-II v0.5.0 (v0.5.0 - 51.2Gb)"
pm$tool[pm$tool == "CONSULT-II (Eq. 7, 140.7Gb)"] <- "CONSULT-II v0.4.0 (v0.4.0 - 140.7Gb)"
pm %>% mutate(old = grepl("v0.5.0 -", tool),tool = sub("[(].*[)]","", tool))  %>% 
  select(rank, tool, metric, p, old) %>%
  pivot_wider(values_from = p, names_from = c(old,metric) ) %>%
  ggplot(
    mapping = aes(x = `FALSE_Bray-Curtis distance`, xend = `TRUE_Bray-Curtis distance`, 
                  y =`FALSE_Shannon equitability` , yend = `TRUE_Shannon equitability`, color=tool)) +
  geom_segment(size = 1, linewidth = 0.75, arrow = arrow(length = unit(0.35, "cm")), show.legend = FALSE) +
  geom_point(size=2) +
  scale_colour_manual(values=c("#e31a1c", "#ff7f00" , "#33a02c", "#1f78b4" ))+ xlim(0, NA) + ylim(0, NA) +
  labs(title = "KRANK uses v0.5.0 of CONSULT-II profiling algorithm instead of v0.4.0", color = "Tool", x = "Bray-Curtis dissimilarity to true profile", y = "Shannon's equitability (absolute diff. to gold standard)") +
  facet_wrap(vars(rank), scales = "free") + theme_minimal_grid(font_size = 17) +   theme(panel.spacing.x = unit(1.5, "lines"), panel.spacing.y = unit(1.5, "lines"), legend.position = "bottom")
ggsave2("../figures/improvement_profiling-new_method.pdf", width = 14, height = 7)

pm <- pma %>% filter(tool %in% c("CONSULT-II (Eq. 7, 140.7Gb)", "CLARK (149.6Gb)", "Bracken (46.5Gb)", "KRANK-hs (Eq. 6, 51.2Gb)", "KRANK-hs (v0.5.0, 51.2Gb)"))
pm$tool[pm$tool == "KRANK-hs (Eq. 6, 51.2Gb)"] <- "KRANK-hs + CONSULT-II v0.5.0 (initial version)"
pm$tool[pm$tool == "KRANK-hs (v0.5.0, 51.2Gb)"] <- "KRANK-hs + CONSULT-II v0.5.0 (current version)"
pm$tool[pm$tool == "Bracken (46.5Gb)"] <- "Bracken v2.8 (46.5Gb)"
pm$tool[pm$tool == "CLARK (149.6Gb)"] <- "CLARK v1.2.6.1 (149.6Gb)"
pm$tool[pm$tool == "CONSULT-II (Eq. 7, 140.7Gb)"] <- "CONSULT-II v0.4.0 (v0.4.0 - 140.7Gb)"
pm %>% mutate(old = grepl("current", tool),tool = sub("[(].*[)]","", tool))  %>% 
  select(rank, tool, metric, p, old) %>%
  pivot_wider(values_from = p, names_from = c(old,metric) ) %>%
  ggplot(
    mapping = aes(x = `FALSE_Bray-Curtis distance`, xend = `TRUE_Bray-Curtis distance`, 
                  y =`FALSE_Shannon equitability` , yend = `TRUE_Shannon equitability`, color=tool)) +
  geom_segment(size = 1, linewidth = 0.75, arrow = arrow(length = unit(0.35, "cm")), show.legend = FALSE) +
  geom_point(size=2) +
  scale_colour_manual(values=c("#e31a1c", "#ff7f00" , "#33a02c", "#1f78b4" ))+ xlim(0, NA) + ylim(0, NA) +
  labs(title = "KRANK using v0.5.0 of CONSULT-II profiling algorithm instead of v0.3.0", color = "Tool", x = "Bray-Curtis dissimilarity to true profile", y = "Shannon's equitability (absolute diff. to gold standard)") +
  facet_wrap(vars(rank), scales = "free") + theme_minimal_grid(font_size = 17) +   theme(panel.spacing.x = unit(1.5, "lines"), panel.spacing.y = unit(1.5, "lines"), legend.position = "bottom")
