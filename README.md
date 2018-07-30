# scrattch.vis

Single-cell RNA-seq visualization tools.  

## Installation

All of the dependencies for `scrattch.vis` are available from CRAN. You can install scrattch.vis from Github using:

```
install.packages("devtools")
devtools::install_github("AllenInstitute/scrattch.vis")
```

## Plot Types

Once installed, `scrattch.vis` provides a variety of functions for visualizing scRNA-seq data based on sample annotations. Each of these returns a `ggplot2` plot object:

### Sample-centric plots:

sample_bar_plot()
![](man/figures/sample_bar_plot.png?raw=true)  

sample_heatmap_plot()
![](man/figures/sample_heatmap_plot.png?raw=true)  

sample_fire_plot()
![]man/figures/sample_fire_plot.png?raw=true)

### Group-centric plots:

group_violing_plot()
![](man/figures/group_violin_plot.png?raw=true)  

group_quasirandom_plot()
![](man/figures/group_quasirandom_plot.png?raw=true)  

group_dot_plot()
![](man/figures/group_dot_plot.png?raw=true)  

group_box_plot()
![](man/figures/group_box_plot.png?raw=true)  

group_heatmap_plot()
![](man/figures/group_heatmap_plot.png?raw=true)  

## The `scrattch` suite

`scrattch.vis` is one component of the [scrattch](https://github.com/AllenInstitute/scrattch/) suite of packages for Single Cell RNA-seq Analysis for Transcriptomic Type CHaracterization from the Allen Institute.

## License

The license for this package is available on Github at: https://github.com/AllenInstitute/scrattch.vis/blob/master/LICENSE

## Level of Support

We are planning on occasional updating this tool with no fixed schedule. Community involvement is encouraged through both issues and pull requests.

## Contribution Agreement

If you contribute code to this repository through pull requests or other mechanisms, you are subject to the Allen Institute Contribution Agreement, which is available in full at: https://github.com/AllenInstitute/scrattch.vis/blob/master/CONTRIBUTION