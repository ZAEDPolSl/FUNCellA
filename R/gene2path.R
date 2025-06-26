#' Transform Gene-level Data to Pathway-level Scores Using Single-Sample Methods
#'
#' This function reduces gene-level data to pathway-level activity scores using a variety of
#' single-sample pathway enrichment methods. It supports filtering pathways based on coverage
#' and size, and allows adjustment of specific scoring methods.
#'
#' @param X A numeric matrix or data.frame with genes/features as rows and samples as columns.
#'   Row names (gene identifiers) must be provided.
#' @param pathway A named list of pathways, where each element is a character vector of gene identifiers.
#' @param method A character string specifying the single-sample enrichment method to use.
#'   Possible values are: \code{"ssGSEA"}, \code{"Mean"}, \code{"BINA"}, \code{"CERNO"}, \code{"ZSc"}, \code{"JASMINE"}, and \code{"AUCell"}.
#'   Default is \code{"ssGSEA"}.
#' @param filt_cov Numeric value between 0 and 1 indicating the minimum fraction of pathway genes
#'   that must be present in \code{X} to include the pathway. Default is 0 (no coverage filtering).
#' @param filt_min Integer specifying the minimum size (number of genes) a pathway must have after mapping to \code{X}.
#'   Default is 15.
#' @param filt_max Integer specifying the maximum size (number of genes) a pathway can have after mapping to \code{X}.
#'   Default is 500.
#' @param type Character string indicating the type of adjustment used in the JASMINE method.
#'   Possible values are \code{"oddsratio"} (default) or \code{"likelihood"}. This parameter is only used when \code{method = "JASMINE"}.
#'
#' @return A data.frame with pathways as rows and samples as columns containing pathway activity scores.
#' @details
#' This function filters and scores pathways using one of several single-sample enrichment methods.
#' Each method uses a different scoring strategy. See below for method descriptions.
#'
#' Method summaries:
#' \describe{
#' \item{ssGSEA}{
#'   Single sample Gene Set Enrichment Analysis (ssGSEA) is a non-parametric method
#'   that computes an enrichment score for each gene set and sample by calculating
#'   the normalized difference between the empirical cumulative distribution functions (CDFs)
#'   of gene expression ranks inside and outside the gene set (\code{\link{pathssGSEA}}).
#'
#'   The implementation used here leverages the \pkg{GSVA} (Hänzelmann et al. (2013)) package and follows the
#'   final normalization step described by Barbie et al. (2009), where pathway scores
#'   are normalized by dividing by the range of calculated values.
#' }
#'
#' \item{Mean}{
#'    Calculates the mean expression of pathway genes per sample, providing a simple aggregate score
#' (\code{\link{pathMean}}). Function reflects \code{addModuleScore} from Seurat package.
#' }
#'
#' \item{BINA}{
#'    Binary scoring method based on the proportion of pathway genes expressed (non-zero) in each sample, and run
#'    logit-transformation to generate scores (\code{\link{pathBINA}}).
#'    }
#'
#' \item{CERNO}{
#'   Implements the CERNO pathway enrichment test (\code{\link{pathCERNO}}), which calculates pathway activity
#'   using a non-parametric approach based on gene expression ranks. For each gene set, the method evaluates whether
#'   the genes in the set are enriched toward the top (or bottom) of a ranked gene list for each sample.
#'   Specifically, it computes an Area Under the Curve (AUC)-like statistic derived from the Mann–Whitney U test.
#'
#'   The AUC score is computed as:
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
#'   }
#'
#' \item{JASMINE}{
#'   Applies the JASMINE scoring method (\code{\link{pathJASMINE}}), designed specifically for single-cell expression data.
#'   It integrates two main components:
#'   \enumerate{
#'     \item Dropout-based gene ranking using \code{\link{Rank_dropout}}, which estimates expression strength based on the rank of non-zero values.
#'     \item Effect size estimation via either odds ratio (\code{\link{Calc_OR}}) or likelihood (\code{\link{Calc_Likelihood}}), capturing differential dropout patterns.
#'   }
#'   Both components are normalized using min-max scaling, and their average is used as the final score.
#'   By default, the method uses the odds ratio for effect size estimation.
#' }
#'
#'  \item{AUCell}{
#'  Computes pathway activity using the AUCell method (\code{\link{pathAUCell}}),
#'    which estimates gene set enrichment per sample by calculating the Area Under the Curve (AUC)
#'    of the recovery curve of gene set members across a ranked gene expression profile.
#'    Conceptually, the AUC reflects the proportion of pathway genes found within the top-ranked
#'    genes of each sample. The number of top genes considered is controlled is set to 5% of the total number of genes.
#'    For details see Aibar et al. (2017).
#'    }
#'
#'  \item{ZSc}{
#'   The Z-score method standardizes gene expression values across samples and
#'   calculates a combined enrichment score for each gene set using the Stouffer integration of Z scores.
#'
#'   Specifically, for a gene set \eqn{G = \{1, \dots, k\}} with standardized
#'   expression values \eqn{z_1, \dots, z_k} for a given sample, the combined
#'   Z-score \eqn{Z_G} is computed as:
#'   \deqn{
#'     Z_G = \frac{\sum_{i=1}^{k} z_i}{\sqrt{k}}
#'   }
#'
#'   This approach is based on Lee et al. (2008) and aggregates standardized gene
#'   expression values to reflect pathway activity while accounting for gene set size. The aggregation itself is known as Stouffer integration method.
#' }
#' }
#' @import cli
#'
#' @examples
#' \dontrun{
#'  # X: gene expression matrix
#'  data(pathways)
#'  pathway_scores <- gene2path(X, pathway, method = "ssGSEA", filt_cov = 0.65)
#' }
#'
#' @references
#' Aibar, S., Bravo González-Blas, C., Moerman, T., Huynh-Thu, V.A., Imrichová, H., Hulselmans, G.,
#' Rambow, F., Marine, J.C., Geurts, P., Aerts, J., van den Oord, J., Kalender Atak, Z., Wouters, J., & Aerts, S. (2017).
#' SCENIC: Single-cell regulatory network inference and clustering.
#' *Nature Methods*, *14*, 1083–1086.
#' \doi{10.1038/nmeth.4463}
#'
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
#' Lee, E., Chuang, H.Y., Kim, J.W., Ideker, T., & Lee, D. (2008).
#' Inferring pathway activity toward precise disease classification.
#' *PLoS Computational Biology*, *4*(11), e1000217.
#' \doi{10.1371/journal.pcbi.1000217}
#'
#' Noureen, N., Ye, Z., Chen, Y., Wang, X., & Zheng, S. (2022).
#' Signature-scoring methods developed for bulk samples are not adequate for cancer single-cell RNA sequencing data.
#' *Elife*, *11*, e71994.
#' \doi{10.7554/eLife.71994}
#'
#' Zyla, J., Marczyk, M., Domaszewska, T., Kaufmann, S. H., Polanska, J., & Weiner III, J. (2019).
#' Gene set enrichment for reproducible science: comparison of CERNO and eight other algorithms.
#' *Bioinformatics*, *35*(24), 5146–5154.
#' \doi{10.1093/bioinformatics/btz447}
#'
#' @export
gene2path<- function(X, pathway, method = "ssGSEA",filt_cov=0,filt_min=15,filt_max=500,type="oddsratio"){
  # ---- Parameter validation ----
  if (missing(X) || is.null(X)) {
    cli_abort(c("x" = "No data provided (X)."))
  }

  if (is.null(rownames(X))) {
    cli_abort(c("x" = "The input assay object doesn't have row names."))
  }

  if (missing(pathway) || is.null(pathway)) {
    cli_abort(c("x" = "No pathway list provided."))
  }

  if (!is.numeric(filt_cov) || filt_cov < 0 || filt_cov > 1) {
    cli_abort(c("x" = "filt_cov must be a numeric value between 0 and 1."))
  }

  valid_methods <- c('ssGSEA', 'Mean', 'BINA', "CERNO", "ZSc", 'JASMINE',"AUCell")
  method <- match.arg(method, choices = valid_methods)

  # ---- Filtering steps ----
  # Filter pathways by size
  pathway<-filter_minmax(pathway,filt_min,filt_max)

  # Coverage-based filtering
  if (filt_cov != 0) {
    pathway <- filter_cover(X, pathway, filt_cov)
  } else {
    invisible(cli_alert_success("No filtration due to coverage"))
  }


  # ---- Gene-2-pathway transformation ----
  switch(method,
      "ssGSEA"= df_path<-pathssGSEA(X,pathway),
      "Mean" =  df_path<-pathMean(X,pathway),
      "JASMINE" = df_path<-pathJASMINE(X,pathway,type),
      "CERNO" = df_path<-pathCERNO(X,pathway),
      "ZSc" = df_path<-pathZScore(X,pathway),
      "BINA" = df_path<-pathBINA(X,pathway),
      "AUCell" = df_path<-pathAUCell(X,pathway)
    )


return(df_path)
}
