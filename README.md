# FUNCellA <img src='man/logo/logo.jpg' align="right" height="140" />
Here is created a package for single-sample pathway enrichment with GMM cell clustering

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GSVA", "AUCell"))

install.packages("devtools")
devtools::install_github("ZAEDPolSl/dpGMM")
devtools::install_github("ZAEDPolSl/SinglePathGMM")

