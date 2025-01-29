#' Function for selecting thresholds ofor pathway activity
#'
#' Function reduce genes to pathways.
#'
#' @param gmms output from GMMdecomp function.
#'
#' @return Function returns...
#'
#' @import cli
#'
#' @export
GMMthr<-function(gmms){

  idprog <- cli_progress_bar("Calculating thresholds", total=length(gmms))
  ret<-list()
  for(i in 1:length(gmms)){
    thr<-gmms[[i]]$threshold
    params<-gmms[[i]]$model

    # Top 1 threshold
    thrT<-max(thr,na.rm=T)

    #K-measn threshold
    if (length(thr)>0){
      # remove component with NA
      if  (sum(is.na(thr))){
        params<-params[-c(which(is.na(thrs))+1),]
        thr<-thr[!is.na(thr)]
      }
      if (nrow(params)==2){
        thrK<-thr
      } else if ((length(thr)+1)!=nrow(params)){
        thrK<-max(thr,na.rm=T)
      }else{
        thrK<-km_search(params,thr)
      }
    } else{thrK<--Inf}


    ret[[i]]<-list(Kmeans_thr=thrK,Top1_thr=thrT,All_thr=thr)
  }
  cli_progress_done(idprog)
  cli_alert_success('Thresholds calculated')

  names(ret)<-names(gmms)
  return(ret)
}
