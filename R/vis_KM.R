#' Visualize k-means thresholding results
#'
#' This function creates visualizations of cluster assignments obtained from
#' \code{thr_KM()}, highlighting both the binary assignment (0/1 for the cluster
#' with the highest centroid) and the full k-means clustering.
#'
#' @param PR A numeric matrix or data.frame with two columns (first and second dimensions
#' of the reduced space, e.g. PCA, UMAP or t-SNE coordinates).
#' @param df_path A numeric matrix or data.frame of pathway activity scores,
#'   where rows correspond to pathways and columns to samples. Row names should correspond to pathway identifiers. Its is output of function gene2path().
#' @param thr_km The result of \code{thr_KM()}, containing both binary (0/1) assignments
#'   and full cluster labels.
#'
#' @return A list of length equal to the number of pathways.
#' Each element is itself a list containing three ggplot2 objects:
#'   \itemize{
#'     \item \code{binary} – points colored by binary assignment (0 = non-selected, 1 = selected cluster).
#'     \item \code{clusters} – points colored by all k-means clusters.
#'     \item \code{PAS} – scatterplot colored by PAS values.
#'   }
#' @export
#'
#' @import ggplot2
#' @import RColorBrewer
#' @import cli
#'
vis_KM <- function(PR, df_path, thr_km) {
  # --- Argument checks ---
  if (missing(PR) || missing(df_path) || missing(thr_km)) {
    cli_abort("All arguments (PR, df_path, thr_km) must be provided.")
  }
  if (!is.data.frame(PR) && !is.matrix(PR)) {
    cli_abort("{.arg PR} must be a data.frame or matrix.")
  }
  if (!is.matrix(df_path) && !is.data.frame(df_path)) {
    cli_abort("{.arg df_path} must be a data.frame or matrix.")
  }
  if (!is.list(thr_km)) {
    cli_abort("{.arg thr_km} must be a list (output from thr_KM).")
  }

  # --- Prep ---
  PR <- as.data.frame(PR)
  colnames(PR) <- c("V1", "V2")

  all_plots <- vector("list", nrow(df_path))

  cli_alert_info("Generating {nrow(df_path)} KM visualizations...")

  for (i in seq_len(nrow(df_path))) {
    df <- PR
    tit <- rownames(df_path)[i]

    # add assignments from thr_KM
    df$binary   <- factor(thr_km$binary[i,], labels = c("0 = non-selected", "1 = selected"))
    df$clusters <- factor(thr_km$clusters[i,], labels = paste0("Cluster ", sort(unique(thr_km$clusters[i,]))))
    df$PAS <- as.numeric(df_path[i, ])

    # binary assignment plot
    p_bin <- ggplot(df, aes(x = V1, y = V2, color = binary)) +
      geom_point(size = 1.2) + theme_bw(base_size = 12) +
      labs(title = paste0(tit, " (binary assignment)"), x = "Dim. 1", y = "Dim. 2") +
      scale_color_manual(values = c("grey70", "red3")) +
      theme(legend.position = "bottom")

    # full cluster assignment plot
    num_clusters <- length(unique(thr_km$clusters[i, ]))
    cluster_palette <- suppressWarnings(colorRampPalette(brewer.pal(min(num_clusters, 8), "Dark2"))(num_clusters))

    p_clust <- ggplot(df, aes(x = V1, y = V2, color = clusters)) +
      geom_point(size = 1.2) + theme_bw(base_size = 12) +
      labs(title = paste0(tit, " (all clusters)"), x = "Dim. 1", y = "Dim. 2") +
      scale_color_manual(values = cluster_palette) +
      theme(legend.position = "bottom")

    # PAS continuous values plot
    p_pas <- ggplot(df, aes(x = V1, y = V2, color = PAS)) +
      geom_point(size = 1.2) + theme_bw(base_size = 12) +
      labs(title = paste0(tit, " (PAS values)"), x = "Dim. 1", y = "Dim. 2") +
      scale_color_viridis_c() +
      theme(legend.position = "bottom")

    all_plots[[i]] <- list(binary = p_bin, clusters = p_clust, PAS = p_pas)
  }

  names(all_plots)<-rownames(df_path)
  cli_alert_success("KM visualizations generated")
  return(all_plots)
}
