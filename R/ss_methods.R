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
  df_path <- ssgseaParam(as.matrix(data), pathway)
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

  df_path <- matrix(NA,length(pathway),ncol(data))
  idprog <- cli_progress_bar("Calculating JASMINE scores", total=length(pathway))
  for (i in 1:length(pathway)){
    df_path[i,]<-as.vector(JASMINE(data,pathway[[i]],type)$JAS_Scores)
    cli_progress_update(id=idprog)
  }
  cli_progress_done(idprog)

  df_path<-as.data.frame(df_path)
  rownames(df_path) <- names(pathway)
  colnames(df_path) <- colnames(data)
  cli_alert_success('JASMINE scores calculated')
  return(df_path)
}
