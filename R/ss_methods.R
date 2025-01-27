#' Function to perform ssGSEA pathway enrichment
#'
#' The function...
#'
#' @param X matrix of data.
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
#' @export
#'
pathssGSEA<-function(X,pathway){
  df_path <- ssgseaParam(as.matrix(X), pathway)
  df_path <- as.data.frame(gsva(df_path))
  rownames(df_path)<-names(pathway)
return(df_path)
}

#' Function to perform Mean pathway enrichment
#'
#' The function...
#'
#' @param X matrix of data.
#' @param pathway list of pathways to analysis.
#'
#' @return A data.frame with pathways in rows and samples in columns.
#'
#'
#' @export
#'
pathMean <- function(X,pathway){
  df_path <- as.data.frame(do.call(rbind, lapply(pathway, function(path) colMeans(extract_pathway(X,path)))))
  rownames(df_path) <- names(pathway)
  return(df_path)
}


#' Function to perform JASMINE pathway enrichment
#'
#' The function...
#'
#' @param X matrix of data.
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
  df_path <- as.data.frame(do.call(rbind, lapply(pathway, function(path) {
    result <- as.vector(JASMINE(X, path, type))
    cli_progress_update(id = idprog)
    return(result)
  })))

  cli_progress_done(idprog)
  rownames(df_path) <- names(pathway)
  colnames(df_path) <- colnames(X)
  cli_alert_success('JASMINE scores calculated')
  return(df_path)
}

#' Function to perform CERNO pathway enrichment
#'
#' The function...
#'
#' @param X matrix of data.
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

  idprog <- cli_progress_bar("Calculating CERNO scores", total=length(pathway))
  dfCerno <- do.call(rbind, lapply(pathway, function(path) {
    df_path <- extract_pathway(X_ranked,path)
    row_AUC <- Calc_AUC(X_ranked,df_path)
    cli_progress_update(id = idprog)
    return(row_AUC)
  }))
  cli_progress_done(idprog)
  dfCerno<-as.data.frame(dfCerno)
  rownames(dfCerno) <- names(pathway)
  cli_alert_success('CERNO scores calculated')
  return(dfCerno)
}
