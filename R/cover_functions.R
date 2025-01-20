#' Function to check coverage of gene-pathway mapping.
#'
#' The function...
#'
#' @param X matrix of data.
#' @param pathway list of pathways to analysis
#' @param filt_cov numeric vector of length 1. Minimum % of pathway coverage after gene mapping. By default is 0 (no filtration), the 0.65 will indicate that only pathway with coverage of genes larger than 65% will be taken.
#'
#' @returns A list of filtered pathways.
cover_check<- function(X,pathway,filt_cov){
  name <- c()
  for (i in 1:length(pathway)){
    df <- X[rownames(X) %in% pathway[[i]],]
    if(nrow(df) != 0){
      cof<-round(nrow(df)/length(pathway[[i]]), digits = 2)
      if(cof < filt_cov){
        print(paste0("Removed path:",names(pathway[i]) ," with coverage: ",cof*100,"%"))
        name <- append(name, names(pathway[i]))
      }
    }else{
      name <- c(name, names(pathway[i]))
    }
  }
  print(paste0('In total removed: ',length(name)," pathways"))
  pathway[name] = NULL
  return(pathway)
}
