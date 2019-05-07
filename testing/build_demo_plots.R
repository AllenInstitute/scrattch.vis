library(tasic2016data)
library(scrattch.vis)
library(scrattch.io)
options(stringsAsFactors = F)

anno <- tasic_2016_anno
anno <- anno[anno$primary_type_id > 0,]
anno <- anno %>%
  mutate(core_intermediate = ifelse(core_intermediate == "intermediate","int",core_intermediate)) %>%
  annotate_cat(core_intermediate)
data <- tasic_2016_rpkm
data_df <- cbind(sample_name = colnames(data),
                 as.data.frame(t(data[c("Pvalb","Sst","Rorb"),])))

sample_bar_plot(data_df, 
                anno, 
                genes = c("Pvalb","Sst","Rorb"), 
                grouping = "primary_type", 
                group_order = c(1,6,2,9,4),
                log_scale = FALSE,
                font_size = 5,
                label_type = "angle")
ggsave("man/figures/sample_bar_plot.png", width = 4, height = 2.5, type = "cairo")

sample_fire_plot(data_df, 
                 anno, 
                 genes = c("Pvalb","Sst","Rorb"), 
                 grouping = "primary_type", 
                 log_scale = TRUE,
                 top_value = "lowest",
                 font_size = 5)
ggsave("man/figures/sample_fire_plot.png", width = 4, height = 2.5, type = "cairo")

sample_heatmap_plot(data_df, 
                    anno, 
                    genes = c("Pvalb","Sst","Rorb"), 
                    grouping = "primary_type", 
                    log_scale = TRUE,
                    font_size = 5)
ggsave("man/figures/sample_heatmap_plot.png", width = 4, height = 2.5, type = "cairo")

group_dot_plot(data_df, 
               anno, 
               genes = c("Pvalb","Sst","Rorb"), 
               grouping = "primary_type", 
               log_scale = TRUE,
               font_size = 5,
               max_size = 5,
               rotate_counts = TRUE)
ggsave("man/figures/group_dot_plot.png", width = 4, height = 2.5, type = "cairo")

group_split_dot_plot(data_df, 
               anno, 
               genes = c("Pvalb","Sst","Rorb"), 
               grouping = "primary_type", 
               split_by = "core_intermediate",
               log_scale = TRUE,
               font_size = 5,
               max_size = 5,
               show_counts = FALSE)
ggsave("man/figures/group_split_dot_plot.png", width = 4, height = 2.5, type = "cairo")

group_violin_plot(data_df, 
                  anno, 
                  genes = c("Pvalb","Sst","Rorb"), 
                  grouping = "primary_type", 
                  log_scale = FALSE,
                  font_size = 5,
                  rotate_counts = TRUE)
ggsave("man/figures/group_violin_plot.png", width = 4, height = 2.5, type = "cairo")

group_quasirandom_plot(data_df, 
                       anno, 
                       genes = c("Pvalb","Sst","Rorb"), 
                       grouping = "primary_type", 
                       log_scale = FALSE,
                       font_size = 5,
                       rotate_counts = TRUE)
ggsave("man/figures/group_quasirandom_plot.png", width = 4, height = 2.5, type = "cairo")

group_box_plot(data_df, 
               anno, 
               genes = c("Pvalb","Sst","Rorb"), 
               grouping = "primary_type", 
               log_scale = FALSE,
               font_size = 5,
               rotate_counts = TRUE)
ggsave("man/figures/group_box_plot.png", width = 4, height = 2.5, type = "cairo")

group_heatmap_plot(data_df, 
                   anno, 
                   genes = c("Pvalb","Sst","Rorb"), 
                   grouping = "primary_type", 
                   stat = "tmean",
                   log_scale = TRUE,
                   font_size = 5,
                   rotate_counts = TRUE)
ggsave("man/figures/group_heatmap_plot.png", width = 4, height = 2.5, type = "cairo")
