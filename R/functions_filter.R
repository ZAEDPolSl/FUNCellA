#' Filter Pathways by Gene Coverage in Expression Data
#'
#' This function filters out pathways from a list based on the fraction of their genes that
#' are present (mapped) in the expression matrix \code{X}. Only pathways meeting the minimum coverage
#' threshold \code{filt_cov} are retained.
#'
#' @param X A numeric matrix or data.frame with genes/features as rows and samples as columns.
#'   Row names (gene identifiers) must be provided.
#' @param pathway A named list of pathways, where each element is a character vector of gene identifiers.
#' @param filt_cov Numeric value between 0 and 1 indicating the minimum fraction of pathway genes
#'   that must be present in \code{X} to include the pathway. Default is 0 (no coverage filtering).
#'
#' @return A list of filtered pathways.
#'
#' @import cli
#'
filter_cover<- function(X,pathway,filt_cov){
  cli_alert_info("PATHWAY COVERAGE FILTRATION")
  name <- c()
  X<-as.matrix(X)
  for (i in seq_along(pathway)){

    poz<-match(pathway[[i]],rownames(X))
    poz<-poz[!is.na(poz)]
    df <- X[poz,]
     if (is.matrix(df) && nrow(df) != 0){
          cof<-round(nrow(df)/length(pathway[[i]]), digits = 2)
          if(cof < filt_cov){
            name <- append(name, names(pathway[i]))
          }
     } else if(is.vector(df)){
       cof<-round(length(df)/length(pathway[[i]]), digits = 2)
       if(cof < filt_cov){
         name <- append(name, names(pathway[i]))
       }
     }else {
       name <- c(name, names(pathway[i]))
    }
  }
  pathway[name] = NULL

  cli_alert_success(paste0('In total removed: ',length(name)," pathways"))
  return(pathway)
}

#' Filter Pathways by Size
#'
#' This function removes pathways from a list based on the number of genes they contain. It is typically used
#' prior to enrichment analysis to exclude pathways that are too small or too large, which may bias the results.
#'
#' @param pathway A named list of pathways, where each element is a character vector of gene identifiers.
#' @param filt_min Integer specifying the minimum size (number of genes) a pathway must have after mapping to \code{X}.
#'   Default is 15.
#' @param filt_max Integer specifying the maximum size (number of genes) a pathway can have after mapping to \code{X}.
#'   Default is 500.
#'
#' @return A list of filtered pathways.
#'
#' @import cli
#'
filter_minmax<- function(pathway,filt_min,filt_max){
  cli_alert_info("PATHWAY SIZE FILTRATION")
  name <- c()
  for (i in 1:length(pathway)){
    siz<-length(pathway[[i]])
    if (siz<filt_min | siz>filt_max){
      name <- append(name, names(pathway[i]))}

  }
  cli_alert_success(paste0('In total removed: ',length(name)," pathways"))
  pathway[name] = NULL
  return(pathway)
}
