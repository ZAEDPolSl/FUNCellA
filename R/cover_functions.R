#' Function to check coverage of gene-pathway mapping.
#'
#' The function...
#'
#' @param X matrix of data.
#' @param pathway list of pathways to analysis
#' @param filt_cov numeric vector of length 1. Minimum % of pathway coverage after gene mapping. By default is 0 (no filtration), the 0.65 will indicate that only pathway with coverage of genes larger than 65% will be taken.
#'
#' @return A list of filtered pathways.
#'
#' @import cli
#'
cover_check<- function(X,pathway,filt_cov){
  cli_alert_info("PATHWAY COVERAGE FILTRATION")
  name <- c()
  for (i in 1:length(pathway)){
    df <- X[rownames(X) %in% pathway[[i]],]
    if(nrow(df) != 0){
      cof<-round(nrow(df)/length(pathway[[i]]), digits = 2)
      if(cof < filt_cov){
        #print(paste0("Removed path:",names(pathway[i]) ," with coverage: ",cof*100,"%"))
        name <- append(name, names(pathway[i]))
      }
    }else{
      name <- c(name, names(pathway[i]))
    }
  }
  pathway[name] = NULL

  cli_alert_success(paste0('In total removed: ',length(name)," pathways"))
  return(pathway)
}

#' Function to remove pathways from analysis by their size.
#'
#' The function...
#'
#' @param pathway list of pathways to analysis
#' @param filt_min numeric vector of length 1. Minimum size of the resulting pathway after gene identifier mapping. By default, the minimum size is 15.
#' @param filt_max numeric vector of length 1. Maximum size of the resulting pathway after gene identifier mapping. By default, the minimum size is 500.
#'
#' @return A list of filtered pathways.
#'
#' @import cli
#'
min_max_filter<- function(pathway,filt_min,filt_max){
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
