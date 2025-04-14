#' Function to rank data without zero counts
#'
#' The function...
#'
#' @param X matrix of data.
#' @param g_vec vector of gene names from pathway.
#'
#' @return A vector of ranked genes.
#'
#' @source \url{https://github.com/NNoureen/JASMINE}
#'
Rank_dropout<- function(X,g_vec){

  subdata = X[X!=0]
  rank_subdata=rank(subdata)
  rank_path <- rank_subdata[match(g_vec, names(rank_subdata), nomatch = 0)]
  CumSum = ifelse(length(rank_path),mean(rank_path,na.rm = TRUE),0 )
  final_rank = CumSum/length(subdata)
  return(final_rank)
}

#' Function to rank data within sample
#'
#' The function...
#'
#' @param X matrix or data.frame of data.
#'
#' @return A data.frame with ranked genes rows and samples in columns.
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
#' The function...
#'
#' @param X matrix of data.
#' @param g_vec vector of gene names from pathway.
#'
#' @return A data.frame with pathway genes in rows and samples in columns.
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


#' Function to calculate AUC for CERNO pathway enrichment
#'
#' The function...
#'
#' @param X matrix of data.
#' @param df data.frame of ranked genes in pathway.
#'
#' @return A vector of CERNO AUC calculation.
#'
#'
Calc_AUC <- function(X,df) {
  row_AUC <- apply(as.matrix(df), 2, function(col) {
    cell <- as.numeric(col)
    AUC <- (nrow(df) * (nrow(X) - nrow(df)) + (nrow(df) * (nrow(df) + 1) / 2) - sum(cell)) / (nrow(df) * (nrow(X) - nrow(df)))
    return(AUC)
  })
  return(row_AUC)
}

#' Function to calculated Odds Ratio in JASMINE pathway enrichment
#'
#' The function...
#'
#' @param X matrix of data.
#' @param g_vec vector of gene names from pathway.
#'
#' @return A vector of odds ratios.
#'
#' @source \url{https://github.com/NNoureen/JASMINE}
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

#' Function to calculated Likelihood in JASMINE pathway enrichment
#'
#' The function...
#'
#' @param X matrix of data.
#' @param g_vec vector of gene names from pathway.
#'
#' @return A vector of likelihoods.
#'
#' @source \url{https://github.com/NNoureen/JASMINE}
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

#' Function to min-max normalization
#'
#' The function...
#'
#' @param x vector of of data.
#'
#' @return A vector of normalized data.
#'
#'
scale_minmax <- function(x)
{
  x_scale = (x - min(x))/(max(x)- min(x))
  return(x_scale)
}

#' Function to z-score normalization
#'
#' The function...
#'
#' @param X matrix or vector of of data.
#'
#' @return A matrix or vector of normalized data.
#'
#' @importFrom stats sd
#'
scale_zscore <- function(X)
{
  row_means <- rowMeans(X)
  row_sds <- apply(X, 1, sd)
  x_scale <- (X - row_means) / row_sds
  return(x_scale)
}

#' Function to calculate JASMINE score pathway enrichment
#'
#' The function...
#'
#' @param X matrix of data.
#' @param g_vec vector of gene names from pathway.
#' @param type type of adjustment of JASMINE score. By default 'oddsratio", another possible input is "likelihood". Parameter only valid for JASMINE method.
#'
#' @return A vector of normalized JASMINE score.
#'
#' @source \url{https://github.com/NNoureen/JASMINE}
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
