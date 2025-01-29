#' Function for GMM decomposition of pathways
#'
#' Function reduce genes to pathways.
#'
#' @param X matrix of data enrichment scores rows pathways columns samples.
#' @param K maximum number of components for GMM decomposition.
#' @param parallel default FALSE. Lunching parallel computing.
#'
#' @return Function returns...
#'
#' @import dpGMM cli
#'
#' @export
GMMdecomp<-function(X, K=10, parallel=F){

  opt<-dpGMM::GMM_1D_opts
  opt$max_iter<-1000
  opt$KS<-K
  opt$plot<-F
  opt$quick_stop<-F
  opt$SW<-0.05
  opt$sigmas.dev<-0

  ret<-list()
  idprog <- cli_progress_bar("Calculating GMM decompositons", total=nrow(X))
  for(i in 1:nrow(X)){
    tmp<-as.numeric(X[i,])
    ret[[i]]<-runGMM(tmp,opts = opt)
    cli_progress_update(id = idprog)
  }
  cli_progress_done(idprog)
  cli_alert_success('GMM calculated')

  names(ret)<-rownames(X)
  return(ret)
}
