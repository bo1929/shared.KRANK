require(ggplot2)
require(readr)
require(cowplot)
require(dplyr)
library(reshape2)

matches <- read_csv("../results/histogramMatchDistances-10kSpecies.txt")
matches <- melt(matches, id = c("Method"))

ggplot(matches, aes(color=Method, x=value)) +
  facet_wrap(c("variable")) +
  stat_ecdf() +
  scale_y_continuous(name = "ECDF") +
  scale_x_continuous(trans = "log10") +
  scale_colour_brewer(palette = "Dark2") +
  theme_cowplot(font_size = 17)

ggplot(matches, aes(color=Method, y=value, x=Method)) +
  facet_wrap(c("variable"), scales = "free") +
  geom_boxplot() +
  geom_point(position = "jitter", alpha = 0.5) +
  scale_y_continuous(trans = "log10") +
  scale_colour_brewer(palette = "Dark2") +
  theme_cowplot(font_size = 17) +
  theme(axis.text.x = element_blank())
