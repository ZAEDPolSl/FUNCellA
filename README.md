# FUNCellA <img src='man/logo/FUNCellaA.png' align="right" height="150" />
The presented package provides a wrapper solution for single-sample pathway enrichment algorithms. Additional functionality includes sample-level thresholding.

Implemented single-sample enrichment algorithms:
1) AUCell scoring (Aibar et al. 2017) - scRNA-Seq only
2) BINA (Zyla et al. 2025?) - scRNA-Seq only
3) CERNO AUC (Zyla et al. 2019)
4) JASMINE (Noureen et al. 2022) - scRNA-Seq only
5) Mean
6) ssGSEA (Barbie et al. 2009, Hänzelmann et al. 2013)
7) Z-score (Lee et al. 2008)

Implemented thresholding solutions:
1) AUCell package thresholding (Aibar et al. 2017)
2) k-means
3) GMM thresholding with Top 1 and k-means adjustment (Zyla et al. 2025?)

## Installation
As the FUNCellA uses other packages from Biocunductor please install at first GSVA and AUCell pacakges.
``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("GSVA", "AUCell"))
```
Next install package from git by devtools.

``` r
install.packages("devtools")
devtools::install_github("ZAEDPolSl/dpGMM")
devtools::install_github("ZAEDPolSl/FUNCellA")
```
## Example of run
### Data
At first let's perpare a data. You need a matrix or data.frame of your molecular biology data e.g. scRNA-Seq counts. You can use your own data or check the example from Seurat Pacakge.
It is very important that in rownames you will have features names (genes/transripts) of the same name as in the pathways you would like to chcek.
``` r
library(FUNCellA)
# data - A numeric matrix or data.frame with genes/features as rows and samples as columns.
# Row names (gene identifiers) must be provided.
rownames(data)
```
Next you need pathways to analysis. You can use example in package or create your own.
``` r
data(pathways)  # examplary pathways in the package
paths <- list(
  path1 = c("MCAM", "THY1", "KLF4", "PDGFRB", "SOX2", "VCAM1", "ITGB1"),
  path2 = c("PECAM1", "CD34", "KIT", "NT5E", "CD44")
)
```
### Single-sample enrichment
Now let's run transformation from genes/feature level to pahway level using gene2path function. By default it will use CERNO approach and apply some filtration; for details see help.
``` r
path_level<-gene2path(data, pathways)
?gene2path
```
### Thresholding
Finally, you can group samples by their pathways activity (each pathway separately). For AUCell package thresholding and simple K-means you use results from previous step.
```r
# extracte only 3 pathwasy to speed up
TA<-thr_AUCell(path_level[1:3,])
TKM<-thr_KM(path_level[1:3,])
```
For GMM solution at first you need to perform decomposition of signal and then perform thresholding search.
```r
# extracte only 3 pathwasy to speed up
res_gmm<-GMMdecomp(path_level[1:3,], K=10, multiply = T)
TGMM<-thr_GMM(res_gmm)
```

## REFERENCES
Aibar, S., Bravo González-Blas, C., Moerman, T., Huynh-Thu, V.A., Imrichová, H., Hulselmans, G., Rambow, F., Marine, J.C., Geurts, P., Aerts, J., van den Oord, J., Kalender Atak, Z., Wouters, J., & Aerts, S (2017). SCENIC: Single-cell regulatory network inference and clustering. *Nature Methods*, *14*, 1083–1086.\
Barbie, D.A., Tamayo, P., Boehm, J.S., et al. (2009). Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1. *Nature*, *462*(7273), 108–112.\
Hänzelmann, S., Castelo, R., & Guinney, J. (2013). GSVA: gene set variation analysis for microarray and RNA-seq data. *BMC Bioinformatics*, *14*, 7.\
Lee, E., Chuang, H.Y., Kim, J.W., Ideker, T., & Lee, D. (2008). Inferring pathway activity toward precise disease classification. *PLoS Computational Biology*, *4*(11), e1000217.\
Noureen, N., Ye, Z., Chen, Y., Wang, X., & Zheng, S. (2022). Signature-scoring methods developed for bulk samples are not adequate for cancer single-cell RNA sequencing data. *Elife*, *11*, e71994.\
Zyla, J., Marczyk, M., Domaszewska, T., Kaufmann, S. H., Polanska, J., & Weiner III, J. (2019). Gene set enrichment for reproducible science: comparison of CERNO and eight other algorithms. *Bioinformatics*, *35*(24), 5146–5154. 

