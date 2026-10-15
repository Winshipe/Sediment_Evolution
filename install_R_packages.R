# Installs every R package used by the notebooks in analysis/.
# Run once:  Rscript install_R_packages.R

cran_pkgs <- c("tidyverse", "here", "rmarkdown", "gt", "colorspace", "cowplot", "ggpubr",
               "gganimate", "gifski", "hexbin", "extrafont", "magick", "rsvg", "remotes", "BiocManager")
missing <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(missing) > 0) install.packages(missing)

bioc_pkgs <- c("treeio", "ggtree")
missing <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(missing) > 0) BiocManager::install(missing, update = FALSE)

if (!requireNamespace("ampvis2", quietly = TRUE)) remotes::install_github("kasperskytte/ampvis2", upgrade = "never")
