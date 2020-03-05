library(tidyverse)
library(gridExtra)
usage_df <- read_csv("../../2020-02-28_sim_hogwash_viz_continuous/data/hogwash_resource_usage.csv", 
                     col_types = c("numeric", "numeric", "numeric", "numeric", "character", "character"))

memory_plot <- usage_df %>% 
  ggplot(aes(x = test, y = `Memory (GB)`)) + 
  geom_boxplot(aes(fill = test)) + 
  geom_jitter() + 
  theme_bw() + 
  xlab("") + 
  theme(legend.position = "none") + 
  ylim(0, max(usage_df$`Memory (GB)`)) + 
  theme(text = element_text(size = 15), 
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) + 
  scale_fill_manual(values = c("dodgerblue", "orchid", "orange"))

time_plot <- usage_df %>% 
  ggplot(aes(x = test, y = `Time (Hours)`)) + 
  geom_boxplot(aes(fill = test)) + 
  geom_jitter() + 
  theme_bw() + 
  xlab("") + 
  ylim(0, max(usage_df$`Time (Hours)`)) + 
  theme(legend.position = "none") + 
  theme(text = element_text(size = 15), 
        axis.text.x = element_text(color = "black"),
        axis.text.y = element_text(color = "black")) + 
  scale_fill_manual(values = c("dodgerblue", "orchid", "orange"))

pdf("../img/hogwash_resource_usage.pdf", width = 8, height = 4)
grid.arrange(memory_plot, time_plot, ncol = 2)
dev.off()
