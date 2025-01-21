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
ssGSEA<-function(X,pathway){
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
calMean <- function(pathways, cells){
  dfMean <- as.data.frame(do.call(rbind, lapply(pathways, function(path) colMeans(Path_df(path, cells)))))
  return(dfMean)
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
#' @export
JAS_path<-function(X,pathway,type="oddsratio"){

  df_path <- matrix(NA,length(pathway),ncol(data))
  for (i in 1:length(pathway)){
    df_path[i,]<-as.vector(JASMINE(data,pathway[[i]],method =type)$JAS_Scores)
  }

  df_path<-as.data.frame(df_path)
  rownames(df_path) <- names(pathway)
  colnames(df_path) <- colnames(data)
  return(df_path)
}
