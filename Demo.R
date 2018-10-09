library(CopulaClust)

Features = data.frame()
for (k in 1:50) {
  x1 = rnorm(1000)
  x2 = rnorm(1000)
  nGrid = 10
  CMatrix = CopulaMatrix(x1, x2, nGrid)
  Features = rbind(Features, data.frame(Matrix2Features(CMatrix)))
  rownames(Features)[nrow(Features)] = paste0("Patient",k)

}

ClusterLabel = CopulaCluster(Features, nClust = 3)


