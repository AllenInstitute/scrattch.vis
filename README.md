# scrattch.vis

Single-cell RNA-seq visualization tools.  

This package is under heavy development and optimization.  

<img src="https://upload.wikimedia.org/wikipedia/commons/c/c0/Animal-Cat-Black-cat-watching-fish-bowl.jpg" alt="Seeing into the fishbowl" width="200px"/>

## Installation

All of the dependencies for `scrattch.vis` are available from CRAN. You can install scrattch.vis from Github using:

```
install.packages("devtools")
devtools::install_github("AllenInstitute/scrattch.vis")
```

## Plot Types

Once installed, `scrattch.vis` provides a variety of functions for visualizing scRNA-seq data based on sample annotations. Each of these returns a `ggplot2` plot object:

Load packages and data:
```
library(tasic2016data)
library(scrattch.io)
library(scrattch.vis)
options(stringsAsFactors = F)

anno <- tasic_2016_anno
anno <- anno[anno$cluster_id > 0,]

anno$core_intermediate <- sub("ermediate","",anno$core_intermediate)
anno <- anno %>%
  annotate_cat(core_intermediate)

data <- tasic_2016_rpkm
data_df <- cbind(sample_name = colnames(data),
                 as.data.frame(t(data[c("Pvalb","Sst","Rorb"),])))
```

### Sample-centric plots:

sample_bar_plot()
```
sample_bar_plot(data_df, 
                anno, 
                genes = c("Pvalb","Sst","Rorb"), 
                grouping = "primary_type", 
                group_order = c(1,6,2,9,4),
                log_scale = FALSE,
                font_size = 5,
                label_type = "angle")
```
![](man/figures/sample_bar_plot.png?raw=true)  

sample_heatmap_plot()
```
sample_heatmap_plot(data_df, 
                    anno, 
                    genes = c("Pvalb","Sst","Rorb"), 
                    grouping = "primary_type", 
                    log_scale = TRUE,
                    font_size = 5)
```
![](man/figures/sample_heatmap_plot.png?raw=true)  

sample_fire_plot()
```
sample_fire_plot(data_df, 
                 anno, 
                 genes = c("Pvalb","Sst","Rorb"), 
                 grouping = "primary_type", 
                 log_scale = TRUE,
                 top_value = "lowest",
                 font_size = 5)
```
![](man/figures/sample_fire_plot.png?raw=true)

### Group-centric plots:

group_violin_plot()
```
group_violin_plot(data_df, 
                  anno, 
                  genes = c("Pvalb","Sst","Rorb"), 
                  grouping = "primary_type", 
                  log_scale = FALSE,
                  font_size = 5,
                  rotate_counts = TRUE)
```

![](man/figures/group_violin_plot.png?raw=true)  

group_quasirandom_plot()
```
group_quasirandom_plot(data_df, 
                       anno, 
                       genes = c("Pvalb","Sst","Rorb"), 
                       grouping = "primary_type", 
                       log_scale = FALSE,
                       font_size = 5,
                       rotate_counts = TRUE)
```
![](man/figures/group_quasirandom_plot.png?raw=true)  

group_dot_plot()
```
group_dot_plot(data_df, 
               anno, 
               genes = c("Pvalb","Sst","Rorb"), 
               grouping = "primary_type", 
               log_scale = TRUE,
               font_size = 5,
               max_size = 5,
               rotate_counts = TRUE)
```
![](man/figures/group_dot_plot.png?raw=true)  

group_split_dot_plot()
```
group_dot_plot(data_df, 
               anno, 
               genes = c("Pvalb","Sst","Rorb"), 
               grouping = "primary_type",
               split_by = "core_intermediate"
               log_scale = TRUE,
               font_size = 5,
               max_size = 5,
               show_counts = FALSE)
```
![](man/figures/group_split_dot_plot.png?raw=true)  

group_box_plot()
```
group_box_plot(data_df, 
               anno, 
               genes = c("Pvalb","Sst","Rorb"), 
               grouping = "primary_type", 
               log_scale = FALSE,
               font_size = 5,
               rotate_counts = TRUE)
```
![](man/figures/group_box_plot.png?raw=true)  

group_heatmap_plot()
```
group_heatmap_plot(data_df, 
                   anno, 
                   genes = c("Pvalb","Sst","Rorb"), 
                   grouping = "primary_type", 
                   stat = "tmean",
                   log_scale = TRUE,
                   font_size = 5,
                   rotate_counts = TRUE)
```
![](man/figures/group_heatmap_plot.png?raw=true)  

## The `scrattch` suite

`scrattch.vis` is one component of the [scrattch](https://github.com/AllenInstitute/scrattch/) suite of packages for Single Cell RNA-seq Analysis for Transcriptomic Type CHaracterization from the Allen Institute.

## License

The license for this package is available on Github at: https://github.com/AllenInstitute/scrattch.vis/blob/master/LICENSE

## Level of Support

We are planning on occasional updating this tool with no fixed schedule. Community involvement is encouraged through both issues and pull requests.

## Contribution Agreement

If you contribute code to this repository through pull requests or other mechanisms, you are subject to the Allen Institute Contribution Agreement, which is available in full at: https://github.com/AllenInstitute/scrattch.vis/blob/master/CONTRIBUTION
