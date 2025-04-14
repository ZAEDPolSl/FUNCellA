#' Function to perform ssGSEA pathway enrichment
#'
#' The function...
#'
#' @param X matrix or data.frame of data (rows: genes/features, columns: samples).
#' @param pathway list of pathways to analysis.
#'
#' @return A data.frame with pathways in rows and samples in columns.
#'
#' @import GSVA
#'
#' @source \url{https://github.com/rcastelo/GSVA}
#' @references Barbie, D.A. et al. Systematic RNA interference reveals that
#' oncogenic KRAS-driven cancers require TBK1.
#' Nature, 462(5):108-112, 2009.
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

#' Function to perform AUCell pathway enrichment
#'
#' The function...
#'
#' @param X matrix or data.frame of data (rows: genes/features, columns: samples).
#' @param pathway list of pathways to analysis.
#'
#' @return A data.frame with pathways in rows and samples in columns.
#'
#' @import AUCell
#' @import cli
#'
#' @references to fill
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


#' Function to perform Mean pathway enrichment
#'
#' The function...
#'
#' @param X matrix or data.frame of data (rows: genes/features, columns: samples).
#' @param pathway list of pathways to analysis.
#'
#' @return A data.frame with pathways in rows and samples in columns.
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

#' Function to perform JASMINE pathway enrichment
#'
#' The function...
#'
#' @param X matrix or data.frame of data (rows: genes/features, columns: samples).
#' @param pathway list of pathways to analysis.
#' @param type type of adjustment of JASMINE score. By default 'oddsratio", another possible input is "likelihood". Parameter only valid for JASMINE method.
#'
#' @return A data.frame with pathways in rows and samples in columns.
#'
#' @source \url{https://github.com/NNoureen/JASMINE}
#' @references Noureen, N. et al. Signature-scoring methods developed for bulk samples are
#' not adequate for cancer single-cell RNA sequencing data. Elife, 2022, 11: e71994.
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

#' Function to perform CERNO pathway enrichment
#'
#' The function...
#'
#' @param X matrix or data.frame of data (rows: genes/features, columns: samples).
#' @param pathway list of pathways to analysis.
#'
#' @return A data.frame with pathways in rows and samples in columns.
#'
#' @references Zyla, J. et al. Gene set enrichment for reproducible science: comparison of CERNO and
#' eight other algorithms. Bioinformatics, 2019, 35.24: 5146-5154.
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


#' Function to perform Z-score integration pathway enrichment
#'
#' The function...
#'
#' @param X matrix or data.frame of data (rows: genes/features, columns: samples).
#' @param pathway list of pathways to analysis.
#'
#' @return A data.frame with pathways in rows and samples in columns.
#'
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

#' Function to perform BINA pathway enrichment
#'
#' The function...
#'
#' @param X matrix or data.frame of data (rows: genes/features, columns: samples).
#' @param pathway list of pathways to analysis.
#'
#' @return A data.frame with pathways in rows and samples in columns.
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
