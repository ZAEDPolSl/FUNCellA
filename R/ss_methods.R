#' Function to perform ssGSEA pathway enrichment
#'
#' The function...
#'
#' @param X matrix of data.
#' @param pathway list of pathways to analysis
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
#' @param pathway list of pathways to analysis
#'
#' @return A data.frame with pathways in rows and samples in columns.
#'
#'
#' @export
Path_Mean <- function(pathway, X){
  df_path <- as.data.frame(do.call(rbind, lapply(pathway, function(path) colMeans(Path_df(path, X)))))
  rownames(df_path) <- names(pathway)
  return(df_path)
}


#' Function to perform JASMINE pathway enrichment
#'
#' The function...
#'
#' @param X matrix of data.
#' @param pathway list of pathways to analysis
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
