#' Matrix Statistics
#'
#' This function computes statistics features of an input Co-occurence Matrix.
#' Features are adapted from matlab toolbox RADIOMICS \url{https://github.com/mvallieres/radiomics/}.
#'
#' @param dataMatrix co-occurent matrix
#' @return a dataframe of co-occurent matrix features
#' @export
#' @references
#' [1] Haralick, R. M., Shanmugam, K., & Dinstein, I. (1973). Textural features for image classification. IEEE Transactions on Systems, Man and Cybernetics, smc 3(6), 610â€“621.
#'
#' [2] Assefa, D., Keller, H., MÃ©nard, C., Laperriere, N., Ferrari, R. J., & Yeung, I. (2010). Robust texture features for response monitoring of glioblastoma multiforme on T1 -weighted and T2 -FLAIR MR images: A preliminary investigation in terms of identification and segmentation. Medical Physics, 37(4), 1722â€“1736.
#'
#' [3] Thibault, G. (2009). Indices de formes et de textures: de la 2D vers la 3D. Application au classement de noyaux de cellules. PhD Thesis, UniversitÃ© AIX-Marseille: p.172.
#'
#' [4] Aerts, H.J.W.L. et al. Decoding tumour phenotype by noninvasive imaging using a quantitative radiomics approach. Nat. Commun. 5:4006 doi: 10.1038/ncomms5006 (2014).

Matrix2Features <- function(dataMatrix){
  Res = list()
  # 1. Energy, Ref.[1]
  Res$Energy = sum(dataMatrix*dataMatrix)

  # 2. Contrast, Ref.[1]
  contrast = 0.0
  nL = dim(dataMatrix)[1]
  for (n in 1:nL-1) {
    temp = 0
    for (i in 1:nL){
      for (j in 1:nL){
        if (abs(i-j)==n)
          temp = temp + dataMatrix[i,j]
      }
    }
    contrast = contrast + n^2*temp;
  }
  Res$contrast = contrast


  # 3. Entropy, Ref.[1]
  tmpMatrix = dataMatrix[dataMatrix>0]
  Res$Entropy = -sum(sum(tmpMatrix*log2(tmpMatrix)));

  # 4. Homogeneity, adapted from Ref.[1]
  temp = 0.0
    for (i in 1:nL){
      for (j in 1:nL){
          temp = temp + dataMatrix[i,j]/(1+abs(i-j))
      }
    }
  Res$Homogeneity = temp


  # 5. Correlation, adapted from Ref. [1] (this definition from MATLAB is preferred from the original one in [1])
  mr = 0
  mc = 0
  for (i in 1:nL){
    for (j in 1:nL){
      mr = mr + i*dataMatrix[i,j]
      mc = mc + j*dataMatrix[i,j]

    }
  }
  Sr = 0
  Sc = 0
  for (i in 1:nL){
    for (j in 1:nL){
      Sr = Sr + (i-mr)^2*dataMatrix[i,j]
      Sc = Sc + (j-mc)^2*dataMatrix[i,j]
    }
  }
  Sr = sqrt(Sr)
  Sc = sqrt(Sc)
  term2 = 0
  for (i in 1:nL){
    for (j in 1:nL){
      term2 = term2 + (i-mr)*(j-mc)*dataMatrix[i,j]
    }
  }
  Corr = term2 / (Sr * Sc);
  Res$Correlation = Corr



  # 6. Variance, Ref.[2]; and 7. SumAverage, Ref.[2]. (adapted from Variance and SumAverage metrics defined by Haralick in Ref. [1])
  # However, in order to compare GLCMs of different sizes, the metrics
  # are divided by the total number of elements in the GLCM (nL*nL). Also,
  # there is probably an error in Assefa's paper [2]: in the variance equation,
  # 'u' should appropriately be replaced by 'ux' and 'uy' as calculated in A1
  # and A2 of the same paper (read 'ui' and 'uj' in our code).
  tempS = 0
  tempV = 0
  for (i in 1:nL){
    for (j in 1:nL){
      tempS = tempS + i*dataMatrix[i,j] + j*dataMatrix[i,j]
      tempV = tempV + (i-mr)^2*dataMatrix[i,j] + (j-mc)^2*dataMatrix[i,j]
    }
  }
  Res$SumAverage = 0.5*sum(tempS)/(nL^2)
  Res$Variance = 0.5*sum(tempV)/(nL^2)

  # 8. Dissimilarity, Ref.[3]
  Dissimilarity = 0
  for (i in 1:nL){
    for (j in 1:nL){
      Dissimilarity = Dissimilarity + abs(i-j)*dataMatrix[i,j]
    }
  }
  Res$Dissimilarity = Dissimilarity


  # 9. AutoCorrelation, Ref.[4]
  AutoCorrelation = 0
  for (i in 1:nL){
    for (j in 1:nL){
      AutoCorrelation = AutoCorrelation + i*j*dataMatrix[i,j]
    }
  }
  Res$AutoCorrelation = AutoCorrelation

  return(Res)

}
