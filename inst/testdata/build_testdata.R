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

