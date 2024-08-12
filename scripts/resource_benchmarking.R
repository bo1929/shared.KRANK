require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)
require(reshape2)
require(tidyr)

rtdf <- read_tsv("../results/running_times-query.tsv")
rtdf$NumQuerySR = rtdf$NumQueries * 66667
rtdf %>% ggplot() + aes(x=NumQuerySR, y=RunningTime, colour=Tool, shape=Tool) +
  geom_line() + geom_point(size=3) +
  scale_color_manual(values=c("#ff7f00", "#33a02c",  "#e31a1c" , "#1f78b4", "#a6cee3")) +
  scale_shape_manual(values=c(18, 17,15,16,16)) +
  theme_minimal_hgrid(font_size = 20) +
  labs(x = "Number of short reads (million)", y = "Running-time (seconds)") +
  scale_x_continuous(trans="log2", breaks = sort(unique(rtdf$NumQuerySR)), labels = c("4.2", "8.4", "16.8", "33.6", "67.1")) +
  scale_y_continuous(trans="log2", n.breaks = 6)

ggsave("../figures/running_time-query.pdf", width=7, height = 6)
