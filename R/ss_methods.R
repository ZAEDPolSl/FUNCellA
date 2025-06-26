#' Wrapper to ssGSEA pathway enrichment
#'
#' The function ryns Single sample Gene Set Enrichment Analysis (ssGSEA) a non-parametric method
#'   that computes an enrichment score for each gene set and sample by calculating
#'   the normalized difference between the empirical cumulative distribution functions (CDFs)
#'   of gene expression ranks inside and outside the gene set (\code{\link{pathssGSEA}}).
#'
#'   The implementation used here leverages the \pkg{GSVA} (Hänzelmann et al. (2013)) package and follows the
#'   final normalization step described by Barbie et al. (2009), where pathway scores
#'   are normalized by dividing by the range of calculated values.
#'
#' @param X A numeric matrix or data.frame with genes/features as rows and samples as columns.
#'   Row names (gene identifiers) must be provided.
#' @param pathway A named list of pathways, where each element is a character vector of gene identifiers.
#'
#' @return A data.frame with pathways as rows and samples as columns containing pathway activity scores.
#'
#' @import GSVA
#'
#' @source \url{https://github.com/rcastelo/GSVA}
#'
#' @references
#' Barbie, D.A., Tamayo, P., Boehm, J.S., et al. (2009).
#' Systematic RNA interference reveals that oncogenic KRAS-driven cancers require TBK1.
#' *Nature*, *462*(7273), 108–112.
#' \doi{10.1038/nature08460}
#'
#' Hänzelmann, S., Castelo, R., & Guinney, J. (2013).
#' GSVA: gene set variation analysis for microarray and RNA-seq data.
#' *BMC Bioinformatics*, *14*, 7.
#' \doi{10.1186/1471-2105-14-7}
#'
#' @import cli
#' @export
#'
pathssGSEA<-function(X,pathway){
  df_enrich <- ssgseaParam(as.matrix(X), pathway)
  df_enrich <- as.data.frame(gsva(df_enrich))
  rownames(df_enrich)<-names(pathway)
  cli_alert_info("ssGSEA enrichment calculated")
return(df_enrich)
}

#' Wrapper to AUCell pathway enrichment
#'
#'  Computes pathway activity using the AUCell method, which estimates gene set enrichment per sample by calculating the Area Under the Curve (AUC)
#'    of the recovery curve of gene set members across a ranked gene expression profile.
#'    Conceptually, the AUC reflects the proportion of pathway genes found within the top-ranked
#'    genes of each sample. The number of top genes considered is controlled is set to 5% of the total number of genes.
#'    For details see Aibar et al. (2017).
#'
#' @param X A numeric matrix or data.frame with genes/features as rows and samples as columns.
#'   Row names (gene identifiers) must be provided.
#' @param pathway A named list of pathways, where each element is a character vector of gene identifiers.
#'
#' @return A data.frame with pathways as rows and samples as columns containing pathway activity scores.
#'
#' @import AUCell
#' @import cli
#'
#' @source \url{https://www.bioconductor.org/packages/release/bioc/html/AUCell.html}
#' @references
#' Aibar, S., Bravo González-Blas, C., Moerman, T., Huynh-Thu, V.A., Imrichová, H., Hulselmans, G.,
#' Rambow, F., Marine, J.C., Geurts, P., Aerts, J., van den Oord, J., Kalender Atak, Z., Wouters, J., & Aerts, S. (2017).
#' SCENIC: Single-cell regulatory network inference and clustering.
#' *Nature Methods*, *14*, 1083–1086.
#' \doi{10.1038/nmeth.4463}
#'
#' @import AUCell
#' @import cli
#'
#' @export
#'
pathAUCell<-function(X,pathway){
  cells_rankings <- AUCell_buildRankings(as.matrix(X))
  AUCs <- AUCell_calcAUC(pathway, cells_rankings, aucMaxRank=0.25*nrow(X), nCores=1)
  df_enrich<-as.data.frame(AUCs@assays@data@listData$AUC)
  rownames(df_enrich)<-names(pathway)
  cli_alert_info("AUCell enrichment calculated")
  return(df_enrich)
}


#' Wrapper to Mean Pathway Enrichment
#'
#' Calculates the mean expression of pathway genes per sample, providing a simple aggregate score.
#' This function reflects the behavior of \code{AddModuleScore} from the \pkg{Seurat} package.
#'
#' @param X A numeric matrix or data.frame with genes/features as rows and samples as columns.
#'   Row names (gene identifiers) must be provided.
#' @param pathway A named list of pathways, where each element is a character vector of gene identifiers.
#'
#' @return A data.frame with pathways as rows and samples as columns containing pathway activity scores.
#'
#'
#' @import cli
#' @export
#'
pathMean <- function(X,pathway){
  idprog <- cli_progress_bar("Calculating mean scores", total=length(pathway))
  df_enrich <- as.data.frame(do.call(rbind, lapply(pathway, function(path) {
    cli_progress_update(id = idprog)
    colMeans(extract_pathway(X, path))
  })))
  cli_progress_done(idprog)
  df_enrich <- as.data.frame(df_enrich)
  rownames(df_enrich) <- names(pathway)
  cli_alert_success("Mean scores calculated")

  return(df_enrich)
}

#' Wrapper to JASMINE Pathway Enrichment Scoring Function
#'
#' This function computes pathway enrichment scores using the JASMINE method to the pathway list. The JASMINE is designed for single-cell data.
#' It integrates a dropout-based gene ranking component with an effect size adjustment using either the odds ratio
#' (default) or a likelihood-based statistic. The scores are normalized via min-max scaling, and the final score is the
#' average of the ranking and the effect size components.
#'
#' @param X A numeric matrix or data.frame with genes/features as rows and samples as columns.
#'   Row names (gene identifiers) must be provided.
#' @param pathway A named list of pathways, where each element is a character vector of gene identifiers.
#' @param type A character string indicating the type of effect size adjustment to use.
#'   Valid options are \code{"oddsratio"} (default) and \code{"likelihood"}.
#'
#' @return A data.frame with pathways as rows and samples as columns containing pathway activity scores.
#'
#' @details
#' The JASMINE score is computed by:
#' \itemize{
#'   \item Ranking expressed genes (non-zero) using a dropout-aware strategy.
#'   \item Calculating an effect size (odds ratio or likelihood) to quantify overrepresentation.
#'   \item Scaling both components via min-max normalization.
#'   \item Averaging the two normalized components for the final score.
#' }
#' If the effect size cannot be computed (e.g., all values are \code{NA}), the score is based solely on the ranking component.
#'
#' @seealso \code{\link{JASMINE}}
#'
#' @source \url{https://github.com/NNoureen/JASMINE}
#' @references
#' Noureen, N., Ye, Z., Chen, Y., Wang, X., & Zheng, S. (2022).
#' Signature-scoring methods developed for bulk samples are not adequate for cancer single-cell RNA sequencing data.
#' *Elife*, *11*, e71994.
#' \doi{10.7554/eLife.71994}
#'
#' @import cli
#' @export
#'
pathJASMINE<-function(X,pathway,type="oddsratio"){

  idprog <- cli_progress_bar("Calculating JASMINE scores", total=length(pathway))
  df_enrich <- as.data.frame(do.call(rbind, lapply(pathway, function(path) {
    result <- as.vector(JASMINE(X, path, type))
    cli_progress_update(id = idprog)
    return(result)
  })))

  cli_progress_done(idprog)
  df_enrich <- as.data.frame(df_enrich)
  rownames(df_enrich) <- names(pathway)
  colnames(df_enrich) <- colnames(X)
  cli_alert_success('JASMINE scores calculated')
  return(df_enrich)
}

#' Wrapper to CERNO pathway enrichment
#'
#' Implements the CERNO pathway enrichment test, which calculates pathway activity
#'   using a non-parametric approach based on gene expression ranks. For each gene set, the method evaluates whether
#'   the genes in the set are enriched toward the top (or bottom) of a ranked gene list for each sample.
#'   Specifically, it computes an Area Under the Curve (AUC)-like statistic derived from the Mann–Whitney U test.
#'
#'  The AUC score is computed as:
#'   \deqn{
#'   \text{AUC}_G = \frac{n_G (n - n_G) + \frac{n_G(n_G + 1)}{2} - R_G}{n_G (n - n_G)}
#'   }
#'   where \( n \) is the total number of genes, \( n_G \) is the number of genes in the pathway,
#'   and \( R_G \) is the sum of their ranks in a given sample.
#'
#'   Scores closer to 1 indicate that genes in the pathway tend to have higher ranks (greater activity),
#'   while values closer to 0 suggest lower relative activity.
#'
#'   This approach enables robust and distribution-free inference of pathway activity, particularly useful
#'   in transcriptomic analyses. Its extensive implementation was presented in Zyla et al. (2019).
#'
#' @param X A numeric matrix or data.frame with genes/features as rows and samples as columns.
#'   Row names (gene identifiers) must be provided.
#' @param pathway A named list of pathways, where each element is a character vector of gene identifiers.
#'
#' @return A data.frame with pathways as rows and samples as columns containing pathway activity scores.
#'
#' @references
#' Zyla, J., Marczyk, M., Domaszewska, T., Kaufmann, S. H., Polanska, J., & Weiner III, J. (2019).
#' Gene set enrichment for reproducible science: comparison of CERNO and eight other algorithms.
#' *Bioinformatics*, *35*(24), 5146–5154.
#' \doi{10.1093/bioinformatics/btz447}
#'
#' @import cli
#' @export
#'
pathCERNO<- function(X, pathway) {
  cli_alert_info("Ranking calculation")
  X_ranked<-Rank_data(X)
  cli_alert_success('Ranks calculated')

  idprog <- cli_progress_bar("Calculating CERNO scores", total=length(pathway))
  df_enrich <- do.call(rbind, lapply(pathway, function(path) {
    df_path <- extract_pathway(X_ranked,path)
    row_AUC <- Calc_AUC(X_ranked,df_path)
    cli_progress_update(id = idprog)
    return(row_AUC)
  }))
  cli_progress_done(idprog)
  df_enrich <- as.data.frame(df_enrich)
  rownames(df_enrich) <- names(pathway)
  cli_alert_success('CERNO scores calculated')
  return(df_enrich)
}


#' Wrapper to Z-score pathway enrichment
#'
#' The Z-score method standardizes gene expression values across samples and
#'   calculates a combined enrichment score for each gene set using the Stouffer integration of Z scores.
#'
#'   Specifically, for a gene set \eqn{G = \{1, \dots, k\}} with standardized
#'   expression values \eqn{z_1, \dots, z_k} for a given sample, the combined
#'   Z-score \eqn{Z_{\G}} is computed as:
#'   \deqn{
#'     Z_{\G} = \frac{\sum_{i=1}^{k} z_i}{\sqrt{k}}
#'   }
#'
#'   This approach is based on Lee et al. (2008) and aggregates standardized gene
#'   expression values to reflect pathway activity while accounting for gene set size. The aggregation itself is known as Stouffer integration method.
#'
#' @param X A numeric matrix or data.frame with genes/features as rows and samples as columns.
#'   Row names (gene identifiers) must be provided.
#' @param pathway A named list of pathways, where each element is a character vector of gene identifiers.
#'
#' @return A data.frame with pathways as rows and samples as columns containing pathway activity scores.
#'
#' @references
#' Lee, E., Chuang, H.Y., Kim, J.W., Ideker, T., & Lee, D. (2008).
#' Inferring pathway activity toward precise disease classification.
#' *PLoS Computational Biology*, *4*(11), e1000217.
#' \doi{10.1371/journal.pcbi.1000217}
#'
#' @import cli
#' @export
#'
pathZScore<- function(X, pathway) {
  cli_alert_info("Data normalization")
  X_normed <- scale_zscore(X)
  cli_alert_success('Data normalized')

  idprog <- cli_progress_bar("Calculating Z-score enrichment", total=length(pathway))
  df_enrich <- do.call(rbind, lapply(pathway, function(path) {
    df_path <- extract_pathway(X_normed,path)
    row_score <- colSums(df_path)/sqrt(nrow(df_path))
    cli_progress_update(id = idprog)
    return(row_score)
  }))

  cli_progress_done(idprog)
  df_enrich <- as.data.frame(df_enrich)
  rownames(df_enrich) <- names(pathway)
  cli_alert_success('Z-score enrichment calculated')
  return(df_enrich)
}

#' Wrapper to BINA pathway enrichment
#'
#'  Binary scoring method based on the proportion of pathway genes expressed (non-zero) in each sample, and run
#'    logit-transformation to generate scores.
#'
#'
#' @param X A numeric matrix or data.frame with genes/features as rows and samples as columns.
#'   Row names (gene identifiers) must be provided.
#' @param pathway A named list of pathways, where each element is a character vector of gene identifiers.
#'
#' @return A data.frame with pathways as rows and samples as columns containing pathway activity scores.
#'
#'
#' @import cli
#' @export
#'
pathBINA <- function(X, pathway) {
  idprog <- cli_progress_bar("Calculating BINA scores", total = length(pathway))
  df_enrich <- do.call(rbind, lapply(pathway, function(path) {
    df <- extract_pathway(X, path)
    row <- colSums((df != 0) / nrow(df))
    row_logit <- log((row + 0.1) / (1 - row + 0.1))
    cli_progress_update(id = idprog)
    return(row_logit)
  }))
  cli_progress_done(idprog)
  df_enrich <- as.data.frame(df_enrich)
  rownames(df_enrich) <- names(pathway)
  colnames(df_enrich) <- colnames(X)
  cli_alert_success("BINA scores calculated")
  return(df_enrich)
}
