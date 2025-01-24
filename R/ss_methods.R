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
#' @export
Path_ssGSEA<-function(X,pathway){
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
Path_Mean <- function(X,pathway){
  df_path <- as.data.frame(do.call(rbind, lapply(pathway, function(path) colMeans(Path_df(path, X)))))
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
#'
#'
#' @import cli
#'
#' @export
Path_JASMINE<-function(X,pathway,type="oddsratio"){


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
#' @import cli
#'
#' @export
Path_Cerno <- function(X, pathway) {
  cli_alert_info("Ranking calculation")
  X_ranked<-Rank_data(X)
  cli_alert_success('Ranks calculated')

  idprog <- cli_progress_bar("Calculating CERNO scores", total=length(pathway))
  dfCerno <- do.call(rbind, lapply(pathway, function(path) {
    df_path <- Path_extract(X_ranked,path)
    row_AUC <- rowAUC(X_ranked,df_path)
    cli_progress_update(id = idprog)
    return(row_AUC)
  }))
  cli_progress_done(idprog)
  dfCerno<-as.data.frame(dfCerno)
  rownames(dfCerno) <- names(pathway)
  cli_alert_success('CERNO scores calculated')
  return(dfCerno)
}
