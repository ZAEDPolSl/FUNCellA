#' Function to performe ssGSEA pathway enrichment
#'
#' The function...
#'
#' @param X matrix of data.
#' @param pathway list of pathways to analsysis
#'
#' @returns A data.frame with pathways in rows and samples in columns.
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
