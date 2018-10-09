#' Clustering by Matrix Features
#'
#' This function cluster patients by extracted matrix features.
#'
#' @param FeatureMatrix n*p features matrix with Patient names as rownames
#' @param nClust number of clusters
#' @return clustering labels
#' @export
#'
CopulaCluster <- function(FeatureMatrix, nClust){
  df = FeatureMatrix

  library(factoextra)

  Dist = get_dist(df, method = "pearson")

  hls <- hclust(d=Dist, method= "complete")
  par(cex=0.8, mar=c(5, 8, 4, 1))
  plot(as.dendrogram(hls), horiz=F)

  hls.clust <- cutree(hls,nClust)
  ClusterLabel = hls.clust

  return(ClusterLabel)
}
