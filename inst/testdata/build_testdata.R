library(tasic2016data)
library(scrattch.vis)

all_genes <- rownames(tasic_2016_counts)

saveRDS(all_genes, file = "all_genes.RData")

genes <- c("Snap25","Pvalb","Rorb",
           "A930038C07Rik","1110059E24Rik","Nkx2-1")

saveRDS(genes, file = "test_genes.RData")

test_matrix <- tasic_2016_counts[genes,]

saveRDS(test_matrix, "test_data_matrix.RData")

test_df <- mat_to_data_df(test_matrix, cols_are = "sample_name")

saveRDS(test_df, "test_data_df.RData")

test_anno <- tasic_2016_anno
names(test_anno)[names(test_anno) == "primary_type"] <- "primary_type_label"

saveRDS(test_anno, "test_anno_incomplete.RData")

set.seed(42)
test_base <- runif(1)
test_sci <- test_base * c(1, 10, 100, 1000, 10000, 100000, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001)
test_sci <- c(0, 1, 10, 100, 1000, 10000, test_sci)
test_sci <- c(test_sci, -1 * test_sci)

test_sci_3_ggplot <- c("0.00%*%10^0","1.00%*%10^0","1.00%*%10^1","1.00%*%10^2","1.00%*%10^3","1.00%*%10^4",
                       "9.15%*%10^-1","9.15%*%10^0","9.15%*%10^1","9.15%*%10^2","9.15%*%10^3","9.15%*%10^4",
                       "9.15%*%10^-2","9.15%*%10^-3","9.15%*%10^-4","9.15%*%10^-5","9.15%*%10^-6","9.15%*%10^-7",
                       "0.00%*%10^0","-1.00%*%10^0","-1.00%*%10^1","-1.00%*%10^2","-1.00%*%10^3","-1.00%*%10^4",
                       "-9.15%*%10^-1","-9.15%*%10^0","-9.15%*%10^1","-9.15%*%10^2","-9.15%*%10^3","-9.15%*%10^4",
                       "-9.15%*%10^-2","-9.15%*%10^-3","-9.15%*%10^-4","-9.15%*%10^-5","-9.15%*%10^-6","-9.15%*%10^-7"
)

test_sci_3_dt <- c("0.00\u271510<sup>0</sup>","1.00\u271510<sup>0</sup>","1.00\u271510<sup>1</sup>","1.00\u271510<sup>2</sup>","1.00\u271510<sup>3</sup>","1.00\u271510<sup>4</sup>",
                   "9.15\u271510<sup>-1</sup>","9.15\u271510<sup>0</sup>","9.15\u271510<sup>1</sup>","9.15\u271510<sup>2</sup>","9.15\u271510<sup>3</sup>","9.15\u271510<sup>4</sup>",
                   "9.15\u271510<sup>-2</sup>","9.15\u271510<sup>-3</sup>","9.15\u271510<sup>-4</sup>","9.15\u271510<sup>-5</sup>","9.15\u271510<sup>-6</sup>","9.15\u271510<sup>-7</sup>",
                   "0.00\u271510<sup>0</sup>","-1.00\u271510<sup>0</sup>","-1.00\u271510<sup>1</sup>","-1.00\u271510<sup>2</sup>","-1.00\u271510<sup>3</sup>","-1.00\u271510<sup>4</sup>",
                   "-9.15\u271510<sup>-1</sup>","-9.15\u271510<sup>0</sup>","-9.15\u271510<sup>1</sup>","-9.15\u271510<sup>2</sup>","-9.15\u271510<sup>3</sup>","-9.15\u271510<sup>4</sup>",
                   "-9.15\u271510<sup>-2</sup>","-9.15\u271510<sup>-3</sup>","-9.15\u271510<sup>-4</sup>","-9.15\u271510<sup>-5</sup>","-9.15\u271510<sup>-6</sup>","-9.15\u271510<sup>-7</sup>"
)

test_sci_3_text <- c("0.00E0","1.00E0","1.00E1","1.00E2","1.00E3","1.00E4",
                     "9.15E-1","9.15E0","9.15E1","9.15E2","9.15E3","9.15E4",
                     "9.15E-2","9.15E-3","9.15E-4","9.15E-5","9.15E-6","9.15E-7",
                     "0.00E0","-1.00E0","-1.00E1","-1.00E2","-1.00E3","-1.00E4",
                     "-9.15E-1","-9.15E0","-9.15E1","-9.15E2","-9.15E3","-9.15E4",
                     "-9.15E-2","-9.15E-3","-9.15E-4","-9.15E-5","-9.15E-6","-9.15E-7"
)

test_sci_label <- list(input = test_sci,
                       output_3_ggplot2 = test_sci_3_ggplot,
                       output_3_dt = test_sci_3_dt,
                       output_3_text = test_sci_3_text)
saveRDS(test_sci_label, file = "test_sci_label_data.RData")


sigline <- sigline() 
saveRDS(sigline, file = "helper_sigline.RData")


plot <- heatmap_legend_plot()
saveRDS(plot, file = "helper_heatmap_legend_plot.RData")


riverplot <- build_river_plot(tasic_2016_anno, c("primary_type", "secondary_type"))
saveRDS(riverplot, file = "helper_riverplot_tasicdata.RData")
