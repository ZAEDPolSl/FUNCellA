#' Function for selecting thresholds for pathway activity
#'
#' Function search for thresholds which indicates samples with relative significant pathway activity.
#' The process is performed for GMM decomposition output.
#'
#' @param gmms output from GMMdecomp function or list of outputs from dpGMM package.
#' Each first list level should refer to the decomposition of one pathway.
#' Second level of list should include variable threshold (vector of thresholds) and model (matrix with model parameters).
#'
#' @return List of thresholds for each investigated pathway based on GMM decomposition.
#'
#' @import cli
#'
#' @export
GMMthr<-function(gmms){

  idprog <- cli_progress_bar("Calculating thresholds", total=length(gmms))
  ret <- vector("list", length(gmms))

  for (i in 1:length(gmms)) {
    cli_progress_update(id = idprog)

    thr<-gmms[[i]]$threshold
    params<-gmms[[i]]$model

    # Top 1 threshold
    thrT<-max(thr,na.rm=T)

    # K-means threshold
    if (length(thr)>0){
      # Remove NA components
      idx <- !is.na(thr)
      thr <- thr[idx]

      if (sum(idx ) < nrow(params)) {
        params <- params[idx, , drop = FALSE]
      }

      if (nrow(params)==2){
        thrK<-thr
      } else if ((length(thr)+1)!=nrow(params)){
        thrK<-max(thr,na.rm=T)
      }else{
        thrK<-km_search(params,thr)
      }
    } else{
      thrK<--Inf
    }


    ret[[i]]<-list(Kmeans_thr=thrK,Top1_thr=thrT,All_thr=thr)

  }
  cli_progress_done(idprog)
  cli_alert_success('Thresholds calculated')

  names(ret)<-names(gmms)
  return(ret)
}
