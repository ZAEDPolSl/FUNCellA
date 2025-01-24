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
  rank_path = rank_subdata[which(names(rank_subdata) %in% g_vec)]
  CumSum = ifelse(length(rank_path),mean(rank_path,na.rm = TRUE),0 )
  final_rank = CumSum/length(subdata)
  return(final_rank)
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

    GE <- X[sig_gene_indices, , drop = FALSE]
    NGE <- X[non_sig_gene_indices, , drop = FALSE]

    SigGenesExp <- colSums(GE != 0)
    NSigGenesExp <- colSums(NGE != 0)

    SigGenesNE <- pmax(nrow(GE) - SigGenesExp, 1)  # Replace 0 with 1
    NSigGenesExp <- pmax(NSigGenesExp, 1)         # Replace 0 with 1
    NSigGenesNE <- pmax(nrow(X) - (NSigGenesExp + SigGenesExp) - SigGenesNE, 0)

    # Calculate the Odds Ratio (Enrichment)
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

  GE <- X[sig_gene_indices, , drop = FALSE]
  NGE <- X[non_sig_gene_indices, , drop = FALSE]

  SigGenesExp <- colSums(GE != 0)
  NSigGenesExp <- colSums(NGE != 0)

  SigGenesNE <- pmax(nrow(GE) - SigGenesExp, 1)  # Replace 0 with 1
  NSigGenesExp <- pmax(NSigGenesExp, 1)         # Replace 0 with 1

  # Calculate the number of not expressed non-signature genes per cell
  NSigGenesNE <- pmax(nrow(X) - (NSigGenesExp + SigGenesExp) - SigGenesNE, 0)

  # Calculate Likelihood Ratios
  LR1 <- SigGenesExp * (NSigGenesExp + NSigGenesNE)
  LR2 <- NSigGenesExp * (SigGenesExp + SigGenesNE)
  LR <- LR1 / LR2

  return(LR)
}

#' Function to normalize JASMINE score
#'
#' The function...
#'
#' @param JAS_Scores vector of of data.
#'
#' @return A vector of normalized JASMINE score.
#'
#' @source \url{https://github.com/NNoureen/JASMINE}
#'
Norm_JASMINE <- function(data)
{
  data_norm = (data - min(data))/(max(data)- min(data))
  return(data_norm)
}

#' Function to calculate JASMINE score pathway enrichment
#'
#' The function...
#'
#' @param X matrix of data.
#' @param g_vec vector of gene names from pathway.
#' @param type type of adjustemnt of JASMINE score. By default 'oddsratio", another possible input is "likelihood". Parameter only valid for JASMINE method.
#'
#' @return A vector of normalized JASMINE score.
#'
#' @source \url{https://github.com/NNoureen/JASMINE}
#'
JASMINE <- function(X,g_vec,type="oddsratio")
{
  g_vec<-unlist(g_vec)
  idx = match(g_vec,rownames(X))
  idx = idx[!is.na(idx)]
  if(length(idx)> 1){
    RM = apply(X,2,function(x) Rank_dropout(x,g_vec))
    RM = Norm_JASMINE(RM)

    if(type == "oddsratio"){
      OR = Calc_OR(X,g_vec)
      JAS_Scores = (RM + OR)/2
    }else if(type == "likelihood"){

      LR = Calc_Likelihood(X,g_vec)
      LR = Norm_JASMINE(LR)
      JAS_Scores = (RM + LR)/2
    }
    FinalScores = data.frame(names(RM),JAS_Scores)
    return(FinalScores)
  }
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
  data_path <- X[rownames(X) %in% g_vec,]
  tmp <- rowVars(as.matrix(data_path))

  data_path <- data_path[!(is.na(tmp)==T | tmp == 0),]
  return(data_path)
}

#' Function to rank data
#'
#' The function...
#'
#' @param X matrix of data.
#'
#' @return A data.frame with ranked genes rows and samples in columns.
#'
#' @importFrom dplyr desc
#'
#'
Rank_data <- function(X){
  rank_df <- as.data.frame(lapply(X, function(col) {
    rank(desc(as.numeric(col)), ties.method = "average")
  }))
  rownames(rank_df) <- rownames(X)
  colnames(rank_df) <- colnames(X)
  return(rank_df)
}

#' Function to calculate AUC for CERNO pathway enrichment
#'
#' The function...
#'
#' @param X matrix of data.
#' @param df data.frame of ranked genes in pathway.
#'
#' @return a vector of CERNO AUC calculation.
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
