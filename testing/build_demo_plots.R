devtools::load_all()

library(tasic2016data)
options(stringsAsFactors = F)

anno <- tasic_2016_anno
data <- tasic_2016_rpkm
data_df <- cbind(sample_name = colnames(data),
                 as.data.frame(t(data[c("Pvalb","Sst","Rorb"),])))


sample_bar_plot(data_df, 
                 anno, 
                 genes = c("Pvalb","Sst","Rorb"), 
                 grouping = "primary_type", 
                 log_scale = FALSE,
                 font_size = 5)
ggsave("man/figures/sample_bar_plot.png", width = 4, height = 2.5)

sample_fire_plot(data_df, 
                  anno, 
                  genes = c("Pvalb","Sst","Rorb"), 
                  grouping = "primary_type", 
                  log_scale = TRUE,
                  top_value = "lowest",
                 font_size = 5)
ggsave("man/figures/sample_fire_plot.png", width = 4, height = 2.5)

sample_heatmap_plot(data_df, 
                 anno, 
                 genes = c("Pvalb","Sst","Rorb"), 
                 grouping = "primary_type", 
                 log_scale = TRUE,
                 font_size = 5)
ggsave("man/figures/sample_heatmap_plot.png", width = 4, height = 2.5)

group_dot_plot(data_df, 
                    anno, 
                    genes = c("Pvalb","Sst","Rorb"), 
                    grouping = "primary_type", 
                    log_scale = TRUE,
                    font_size = 5,
               rotate_counts = TRUE)
ggsave("man/figures/group_dot_plot.png", width = 4, height = 2.5)

group_violin_plot(data_df, 
               anno, 
               genes = c("Pvalb","Sst","Rorb"), 
               grouping = "primary_type", 
               log_scale = FALSE,
               font_size = 5,
               rotate_counts = TRUE)
ggsave("man/figures/group_violin_plot.png", width = 4, height = 2.5)

group_quasirandom_plot(data_df, 
                  anno, 
                  genes = c("Pvalb","Sst","Rorb"), 
                  grouping = "primary_type", 
                  log_scale = FALSE,
                  font_size = 5,
                  rotate_counts = TRUE)
ggsave("man/figures/group_quasirandom_plot.png", width = 4, height = 2.5)

group_box_plot(data_df, 
                  anno, 
                  genes = c("Pvalb","Sst","Rorb"), 
                  grouping = "primary_type", 
                  log_scale = FALSE,
                  font_size = 5,
               rotate_counts = TRUE)
ggsave("man/figures/group_box_plot.png", width = 4, height = 2.5)

group_heatmap_plot(data_df, 
               anno, 
               genes = c("Pvalb","Sst","Rorb"), 
               grouping = "primary_type", 
               stat = "tmean",
               log_scale = TRUE,
               font_size = 5,
               rotate_counts = TRUE)
ggsave("man/figures/group_heatmap_plot.png", width = 4, height = 2.5)
