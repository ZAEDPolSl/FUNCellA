#' Function to transform genes to pathway level by single-sample approaches
#'
#' Function reduce genes to pathways.
#'
#' @param X matrix or data.frame of data (rows: genes/features, columns: samples).
#' @param pathway list of pathways to analysis.
#' @param method methods for single-sample pathway enrichment analysis.
#' @param filt_cov numeric vector of length 1. Minimum \% of pathway coverage after gene mapping. By default is 0 (no filtration), the 0.65 will indicate that only pathway with coverage of genes larger than 65\% will be taken.
#' @param filt_min numeric vector of length 1. Minimum size of the resulting pathway after gene identifier mapping. By default, the minimum size is 15.
#' @param filt_max numeric vector of length 1. Maximum size of the resulting pathway after gene identifier mapping. By default, the minimum size is 500.
#' @param type type of adjustment of JASMINE score. By default 'oddsratio", another possible input is "likelihood". Parameter only valid for JASMINE method.
#'
#' @return Function returns...
#'
#' @import cli
#'
#'
#' @export
gene2path<- function(X, pathway, method = "ssGSEA",filt_cov=0,filt_min=15,filt_max=500,type="oddsratio"){
  # ---- Parameter validation ----
  if (missing(X) || is.null(X)) {
    cli_abort(c("x" = "No data provided (X)."))
  }

  if (is.null(rownames(X))) {
    cli_abort(c("x" = "The input assay object doesn't have row names."))
  }

  if (missing(pathway) || is.null(pathway)) {
    cli_abort(c("x" = "No pathway list provided."))
  }

  if (!is.numeric(filt_cov) || filt_cov < 0 || filt_cov > 1) {
    cli_abort(c("x" = "filt_cov must be a numeric value between 0 and 1."))
  }

  valid_methods <- c('ssGSEA', 'Mean', 'BINA', "CERNO", "ZSc", 'JASMINE')
  method <- match.arg(method, choices = valid_methods)

  # ---- Filtering steps ----
  # Coverage-based filtering
  if (filt_cov != 0) {
    pathway <- filter_cover(X, pathway, filt_cov)
  } else {
    invisible(cli_alert_success("No filtration due to coverage"))
  }
  # Filter pathways by size
  pathway<-filter_minmax(pathway,filt_min,filt_max)

  # ---- Gene-2-pathway transformation ----
  switch(method,
      "ssGSEA"= df_path<-pathssGSEA(X,pathway),
      "Mean" =  df_path<-pathMean(X,pathway),
      "JASMINE" = df_path<-pathJASMINE(X,pathway,type),
      "CERNO" = df_path<-pathCERNO(X,pathway),
      "ZSc" = df_path<-pathZScore(X,pathway),
      "BINA" = df_path<-pathBINA(X,pathway)
    )


return(df_path)
}
