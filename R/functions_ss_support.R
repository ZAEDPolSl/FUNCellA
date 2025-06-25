#' Rank Genes Excluding Zero Counts
#'
#' This function ranks non-zero gene expression values within a sample and calculates
#' the average normalized rank for a given set of genes (pathway).
#' It excludes zero counts (dropouts) from ranking to focus on expressed genes.
#'
#' @param X A numeric matrix or data.frame with genes/features as rows and samples as columns.
#' @param g_vec Character vector of gene names representing the gene set/pathway.
#'
#' @return A numeric value representing the average normalized rank of genes in \code{g_vec}
#'   relative to all non-zero genes in \code{X}.
#'
#' @source \url{https://github.com/NNoureen/JASMINE}
#'
#' @references
#' Noureen, N., Ye, Z., Chen, Y., Wang, X., & Zheng, S. (2022).
#' Signature-scoring methods developed for bulk samples are not adequate for cancer single-cell RNA sequencing data.
#' *Elife*, *11*, e71994.
#' \doi{10.7554/eLife.71994}
#'
Rank_dropout<- function(X,g_vec){

  subdata = X[X!=0]
  rank_subdata=rank(subdata)
  rank_path <- rank_subdata[match(g_vec, names(rank_subdata), nomatch = 0)]
  CumSum = ifelse(length(rank_path),mean(rank_path,na.rm = TRUE),0 )
  final_rank = CumSum/length(subdata)
  return(final_rank)
}

#' Rank Genes Within Each Sample
#'
#' This function ranks gene expression values across genes within each sample (column-wise),
#' using average ranks for ties.
#'
#' @param X A numeric matrix or data.frame with genes/features as rows and samples as columns.
#'
#' @return A data.frame of the same dimensions as \code{X}, containing ranked values
#'   for each gene within each sample (higher expression values receive lower ranks).
#'
#' @importFrom matrixStats colRanks
#'
#'
Rank_data <- function(X) {
  X<--X
  rank_matrix<- colRanks(as.matrix(X), ties.method="average", preserveShape=TRUE, decreasing=TRUE)
  return(as.data.frame(rank_matrix, row.names = rownames(X)))
}

#' Function to extract matrix of genes in pathway
#'
#' This function extracts and returns the expression matrix corresponding to genes in a specified pathway.
#' Genes with zero or undefined variance across samples are removed to ensure informative input
#' for downstream analysis.
#'
#' @param X A numeric matrix or data.frame with genes/features as rows and samples as columns.
#'   Row names (gene identifiers) must be provided.
#' @param g_vec Character vector of gene names representing the gene set/pathway.
#'
#' @return A data.frame containing only the genes in \code{g_vec} that are present in \code{X},
#'   excluding genes with zero or missing variance.
#'
#' @importFrom matrixStats rowVars
#'
#'
extract_pathway <- function(X,g_vec){
  g_vec<-unlist(g_vec)
  poz<-match(g_vec,rownames(X))
  poz<-poz[!is.na(poz)]

  data_path <- X[poz,]
  tmp <- rowVars(as.matrix(data_path))

  data_path <- data_path[!(is.na(tmp)==T | tmp == 0),]
  return(data_path)
}


#' Calculate AUC for CERNO pathway enrichment
#'
#' This function calculates the Area Under the Curve (AUC) for pathway gene ranks
#' using the CERNO test, which is based on the Mann–Whitney U-statistic. It compares
#' the distribution of gene ranks in a pathway against all other genes in the ranked expression matrix.
#'
#'
#' @param X A numeric matrix or data.frame of ranked expression values, where rows are genes/features
#'   and columns are samples. Typically obtained via \code{\link{Rank_data}}.
#' @param df A data.frame containing a subset of \code{X} with ranked expression values for genes
#'   in a specific pathway (i.e., a result of \code{\link{extract_pathway}} applied to ranked data).
#'
#' @return A numeric vector of AUC values, one per sample, representing CERNO enrichment scores.
#'
#' @details
#' The AUC score is computed as:
#'   \deqn{
#'   \text{AUC}_G = \frac{n_G (n - n_G) + \frac{n_G(n_G + 1)}{2} - R_G}{n_G (n - n_G)}
#'   }
#'   where n is the total number of genes, n_G is the number of genes in the pathway,
#'   and R_G is the sum of their ranks in a given sample.
#'
#'   Scores closer to 1 indicate that genes in the pathway tend to have higher ranks (greater activity),
#'   while values closer to 0 suggest lower relative activity.
#'
#' @seealso \code{\link{pathCERNO}}, \code{\link{Rank_data}}, \code{\link{extract_pathway}}
#'
#' @references
#' Zyla, J., Marczyk, M., Domaszewska, T., Kaufmann, S. H., Polanska, J., & Weiner III, J. (2019).
#' Gene set enrichment for reproducible science: comparison of CERNO and eight other algorithms.
#' *Bioinformatics*, *35*(24), 5146–5154.
#' \doi{10.1093/bioinformatics/btz447}
#'
Calc_AUC <- function(X,df) {
  row_AUC <- apply(as.matrix(df), 2, function(col) {
    cell <- as.numeric(col)
    AUC <- (nrow(df) * (nrow(X) - nrow(df)) + (nrow(df) * (nrow(df) + 1) / 2) - sum(cell)) / (nrow(df) * (nrow(X) - nrow(df)))
    return(AUC)
  })
  return(row_AUC)
}

#' Calculate Odds Ratios for JASMINE Pathway Enrichment
#'
#' This function computes the odds ratio (OR) for each sample, assessing the enrichment
#' of a pathway by comparing the presence (non-zero expression) of genes in the pathway
#' versus those not in the pathway. It is used as an effect size component in the
#' JASMINE scoring method for single-cell data.
#'
#' @param X A numeric matrix or data.frame with genes/features as rows and samples as columns.
#'   Row names (gene identifiers) must be provided.
#' @param g_vec Character vector of gene names representing the gene set/pathway.
#'
#' @return A numeric vector of odds ratios, one per sample. Higher values indicate stronger
#'   enrichment of the pathway's genes in the expressed gene set of each sample.
#'
#' @details
#' The odds ratio is computed based on the contingency table of expressed vs. non-expressed
#' genes for both pathway and non-pathway genes. A pseudocount is used to prevent division
#' by zero in cases with no expression.
#'
#' @source \url{https://github.com/NNoureen/JASMINE}
#'
#' @references
#' Noureen, N., Ye, Z., Chen, Y., Wang, X., & Zheng, S. (2022).
#' Signature-scoring methods developed for bulk samples are not adequate for cancer single-cell RNA sequencing data.
#' *Elife*, *11*, e71994.
#' \doi{10.7554/eLife.71994}
#'
Calc_OR <- function(X,g_vec){
  sig_gene_indices <- which(rownames(X) %in% g_vec)
  non_sig_gene_indices <- setdiff(seq_len(nrow(X)), sig_gene_indices)

  SigGenesExp <- colSums(X[sig_gene_indices, , drop = FALSE] != 0)
  NSigGenesExp <- colSums(X[non_sig_gene_indices, , drop = FALSE] != 0)

  SigGenesNE <- pmax(length(sig_gene_indices) - SigGenesExp, 1)
  NSigGenesNE <- pmax(nrow(X) - (NSigGenesExp + SigGenesExp) - SigGenesNE, 0)

  OR <- (SigGenesExp * NSigGenesNE) / (SigGenesNE * NSigGenesExp)
  return(OR)
}

#' Function to calculate Likelihood score in JASMINE pathway enrichment
#'
#' This function computes the likelihood ratio statistic for pathway enrichment,
#' comparing the expression of genes in the provided gene set  against non-pathway genes
#' across all samples. The output reflects the relative likelihood that genes in the pathway
#' are enriched in the non-zero expression profile of each sample.
#'
#' @param X A numeric matrix or data.frame with genes/features as rows and samples as columns.
#'   Row names (gene identifiers) must be provided.
#' @param g_vec Character vector of gene names representing the gene set/pathway.
#'
#' @return A numeric vector of likelihood ratios, one per sample. Higher values indicate stronger
#'   enrichment of the pathway's genes in the expressed gene set of each sample.
#'
#' @source \url{https://github.com/NNoureen/JASMINE}
#'
#' @references
#' Noureen, N., Ye, Z., Chen, Y., Wang, X., & Zheng, S. (2022).
#' Signature-scoring methods developed for bulk samples are not adequate for cancer single-cell RNA sequencing data.
#' *Elife*, *11*, e71994.
#' \doi{10.7554/eLife.71994}
#'
Calc_Likelihood <- function(X, g_vec) {
  sig_gene_indices <- which(rownames(X) %in% g_vec)
  non_sig_gene_indices <- setdiff(seq_len(nrow(X)), sig_gene_indices)

  SigGenesExp <- colSums(X[sig_gene_indices, , drop = FALSE] != 0)
  NSigGenesExp <- colSums(X[non_sig_gene_indices, , drop = FALSE] != 0)

  SigGenesNE <- pmax(length(sig_gene_indices) - SigGenesExp, 1)
  NSigGenesNE <- pmax(nrow(X) - (NSigGenesExp + SigGenesExp) - SigGenesNE, 0)

  LR1 <- SigGenesExp * (NSigGenesExp + NSigGenesNE)
  LR2 <- NSigGenesExp * (SigGenesExp + SigGenesNE)

  LR <- LR1 / LR2

  return(LR)
}

#' Function for Min-Max Normalization
#'
#' This function applies min-max normalization to a numeric vector, scaling the values to the range [0, 1].
#' It is commonly used to standardize feature values before further analysis or scoring.
#'
#' @param x A numeric vector of data to be normalized.
#'
#' @return A numeric vector of the same length as \code{x}.
#'
#' @examples
#' scale_minmax(c(2, 4, 6, 8))
#' # Returns: c(0, 0.333, 0.667, 1)
#'
scale_minmax <- function(x)
{
  x_scale = (x - min(x))/(max(x)- min(x))
  return(x_scale)
}

#' Function for Z-Score Normalization
#'
#' This function applies z-score normalization (standardization) across each row of a numeric matrix or to a numeric vector.
#' For matrices, it standardizes each row to have a mean of 0 and standard deviation of 1.
#' For vectors, it standardizes the values globally.
#'
#' @param X A numeric matrix or vector. For matrices, rows are standardized independently.
#'
#' @return A numeric matrix or vector of the same shape as \code{X}, containing z-score normalized values.
#'
#' @importFrom stats sd
#'
#' @examples
#' mat <- matrix(1:9, nrow = 3, byrow = TRUE)
#' scale_zscore(mat)
#'
#' vec <- c(1, 2, 3, 4, 5)
#' scale_zscore(vec)
#'
scale_zscore <- function(X) {
  if (is.vector(X)) {
    x_mean <- mean(X)
    x_sd <- sd(X)
    return((X - x_mean) / x_sd)
  } else {
    row_means <- rowMeans(X)
    row_sds <- apply(X, 1, sd)
    x_scale <- (X - row_means) / row_sds
    return(x_scale)
  }
}

#' JASMINE Pathway Enrichment Scoring Function
#'
#' This function computes pathway enrichment scores using the JASMINE method, which is designed for single-cell data.
#' JASMINE integrates a dropout-based gene ranking component with an effect size adjustment using either the odds ratio
#' (default) or a likelihood-based statistic. The scores are normalized via min-max scaling, and the final score is the
#' average of the ranking and the effect size components.
#'
#' @param X A numeric matrix or data.frame with genes/features as rows and samples as columns.
#'   Row names (gene identifiers) must be provided.
#' @param g_vec Character vector of gene names representing the gene set/pathway.
#' @param type A character string indicating the type of effect size adjustment to use.
#'   Valid options are \code{"oddsratio"} (default) and \code{"likelihood"}.
#'
#' @return A numeric vector of normalized JASMINE scores, one per sample/cell. Higher scores indicate greater pathway activity.
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

#' @source \url{https://github.com/NNoureen/JASMINE}
#'
#' @references
#' Noureen, N., Ye, Z., Chen, Y., Wang, X., & Zheng, S. (2022).
#' Signature-scoring methods developed for bulk samples are not adequate for cancer single-cell RNA sequencing data.
#' *Elife*, *11*, e71994.
#' \doi{10.7554/eLife.71994}
#'
JASMINE <- function(X,g_vec,type="oddsratio")
{
  g_vec<-unlist(g_vec)
  idx <- which(rownames(X) %in% g_vec)
  if(length(idx)> 1){
    RM = apply(X,2,function(x) Rank_dropout(x,g_vec))
    RM = scale_minmax(RM)

    # effect size calculation
    if(type == "oddsratio"){
      ES = Calc_OR(X,g_vec)
    }else if(type == "likelihood"){
      ES = Calc_Likelihood(X,g_vec)
    }
    ES = scale_minmax(ES)

    if (sum(is.na(ES))==length(ES)){
      FinalScores = RM
    } else{
      FinalScores = (RM + ES)/2
    }
    return(FinalScores)
  }
}
