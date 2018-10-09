#' Scatter data to Copula
#'
#' This function tranforms orginal joint distribution to Copula and extract discret co-occurent matrix.
#' @param x1 Value array of first modality
#' @param x2 Value array of second modality
#' @return An discrete copula matrix with the size of nGrid*nGrid
#' @export
#'
CopulaMatrix <- function(x1, x2, nGrid){

  library(copula)
  library(squash)

  dx = 1/nGrid

  data = as.matrix(cbind(x1,x2))
  copuladata = pobs(data)

  u <- copuladata[,1]
  v <- copuladata[,2]

  tmp = hist2(x= u, y= v, xbreaks = seq(0,1,dx), ybreaks = seq(0,1,dx), plot = FALSE)
  tmp$z[is.na(tmp$z)] = 0

  Res = t(tmp$z)/sum(tmp$z)
  colnames(Res) <- paste('u', colnames(Res))
  rownames(Res) <- paste('v', rownames(Res))

  return(Res)

}
