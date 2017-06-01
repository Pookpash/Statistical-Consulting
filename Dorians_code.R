conv2mat <- function(plist, N, var = F){
  plist <- plist[1:(length(plist)-2)]
  vecm <- c(plist$mu1,plist$mu2,plist$mu3,plist$mu,plist$theta)
  matm <- matrix(vecm, ncol = N, nrow = (length(plist)/2),byrow=T)
  vecv <- c(plist$sigma1,plist$sigma2,plist$sigma3,plist$kappa,plist$phi)
  matv <- matrix(vecv, ncol = N, nrow = (length(plist)/2),byrow=T)
  if (var == T){
    result <- matv
  }
  else {
    result <- matm
  }
  return(result)
}
