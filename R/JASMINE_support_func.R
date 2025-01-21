#' Function to rank data on JASMINE pathway enrichment
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
RankCalculation <- function(X,g_vec){

  subdata = X[X!=0]                                                                      ### Removing Dropouts from single cell
  DataRanksUpdated=rank(subdata)                                                         ### Calculating ranks of each signature gene per cell
  DataRanksSigGenes = DataRanksUpdated[which(names(DataRanksUpdated) %in% g_vec)]        ### Shortling rank vector for signature genes
  CumSum = ifelse(length(DataRanksSigGenes),mean(DataRanksSigGenes,na.rm = TRUE),0 )     ### Calculating Mean of ranks for signature genes
  FinalRawRank = CumSum/length(subdata)                                                  ### Normalizing Means by total coverage
  return(FinalRawRank)
}

#' Function to calculated Odds Ratio in JASMINE pathway enrichment
#'
#' The function...
#'
#' @param X matrix of data.
#' @param g_vec vector of gene names from pathway.
#'
#' @return A vector of odds ratio.
#'
#' @source \url{https://github.com/NNoureen/JASMINE}
#'
ORCalculation <- function(X,g_vec){
  GE = X[which(rownames(X) %in% g_vec),]                                          ### Subsetting data for signature genes
  NGE = X[-which(rownames(X) %in% g_vec),]                                        ### Subsetting data for non-signature genes
  SigGenesExp = apply(GE,2,function(x) length(x[x!=0]))                                 ### Calculating Number of expressed Signature Genes per cell
  NSigGenesExp = apply(NGE,2,function(x) length(x[x!=0]))                               ### Calculating Number of expressed Non-Signature Genes per cell
  SigGenesNE = nrow(GE) - SigGenesExp                                                   ### Calculating Number of Not expressed Signature Genes per cell
  SigGenesNE = replace(SigGenesNE,SigGenesNE==0,1)									                    ### Replacing Zero's with 1
  NSigGenesExp = replace(NSigGenesExp,NSigGenesExp==0,1)                                ### Replacing Zero's with 1
  NSigGenesNE = nrow(X) - (NSigGenesExp + SigGenesExp)                               ### Calculating Number of Not expressed Non-Signature Genes per cell
  NSigGenesNE = NSigGenesNE - SigGenesNE
  OR = (SigGenesExp * NSigGenesNE) / (SigGenesNE * NSigGenesExp)                        ### Calculating Enrichment (Odds Ratio)
  return(OR)
}

#' Function to calculated Likelihood in JASMINE pathway enrichment
#'
#' The function...
#'
#' @param X matrix of data.
#' @param g_vec vector of gene names from pathway.
#'
#' @return A vector of likelihood.
#'
#' @source \url{https://github.com/NNoureen/JASMINE}
#'
LikelihoodCalculation <- function(X,g_vec){
  GE = X[which(rownames(X) %in% g_vec),]
  NGE = X[-which(rownames(X) %in% g_vec),]
  SigGenesExp = apply(GE,2,function(x) length(x[x!=0]))
  NSigGenesExp = apply(NGE,2,function(x) length(x[x!=0]))
  SigGenesNE = nrow(GE) - SigGenesExp
  SigGenesNE = replace(SigGenesNE,SigGenesNE==0,1)
  NSigGenesExp = replace(NSigGenesExp,NSigGenesExp==0,1)
  NSigGenesNE = nrow(X) - (NSigGenesExp + SigGenesExp)
  NSigGenesNE = NSigGenesNE - SigGenesNE
  LR1 = SigGenesExp*(NSigGenesExp + NSigGenesNE)
  LR2 = NSigGenesExp * (SigGenesExp + SigGenesNE)
  LR = LR1/LR2
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
NormalizationJAS <- function(JAS_Scores)
{
  JAS_Scores = (JAS_Scores - min(JAS_Scores))/(max(JAS_Scores)- min(JAS_Scores))
  return(JAS_Scores)
}

#' Function to calculate JASMINE score pathway enrichment
#'
#' The function...
#'
#' @param X matrix of data.
#' @param g_vec vector of gene names from pathway.
#'
#' @return A vector of normalized JASMINE score.
#'
#' @source \url{https://github.com/NNoureen/JASMINE}
#'
JASMINE <- function(X,g_vec,type="oddsratio")
{
  idx = match(g_vec,rownames(X))
  idx = idx[!is.na(idx)]
  if(length(idx)> 1){
    RM = apply(X,2,function(x) RankCalculation(x,g_vec))                              ### Mean RankCalculation for single cell data matrix
    RM = NormalizationJAS(RM)                                                            ### Normalizing Mean Ranks

    if(type == "oddsratio"){
      OR = ORCalculation(X,g_vec)			                                             ### Signature Enrichment Calculation for single cell data matrix (OR)
      OR = NormalizationJAS(OR)															 ### Normalizing Enrichment Scores (OR)
      JAS_Scores = (RM + OR)/2
    }else if(type == "likelihood"){

      LR = LikelihoodCalculation(X,g_vec)			                                     ### Signature Enrichment Calculation for single cell data matrix  (LR)
      LR = NormalizationJAS(LR)															 ### Normalizing Enrichment Scores (LR)
      JAS_Scores = (RM + LR)/2
    }
    FinalScores = data.frame(names(RM),JAS_Scores)                                       ### JASMINE scores
    colnames(FinalScores)[1]='SampleID'
    return(FinalScores)
  }
}
