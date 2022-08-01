
#' @title  Compute TSGS moments
#' 
#' @description 
#' This is a function making use of TSGS function from GSE package to compute
#' the Two-Step Generalized S-Estimate, a robust estimate of location 
#' and scatter for data with cell-wise and case-wise contamination.
#' 
#' @param R xts object of asset returns
#' @param filter the filter to be used in the first step. Available choices are 
#'               "UBF-DDC","UBF","DDC","UF". The default one is "UBF-DDC".
#' @param partial.impute whether partial imputation is used prior to estimation.
#'                       The default is FALSE.
#' @param tol tolerance for the convergence criterion. Default is 1e-4.
#' @param maxiter maximum number of iterations. Default is 150.
#' @param method loss function to use, "bisquare" or "rocke". Default is "bisquare"
#' @param init type of initial estimator. Options include "emve", "qc", "huber","imputed","emve_c"
#'
#' @return estimators of first and second moments
#' 
#' @references Claudio Agostinelli, Andy Leung, "Robust estimation of multivariate 
#'             location and scatter in the presence of cellwise and casewise contamination",
#'             2014.
#' @export
#' @examples
#' 
tsgs.moments <- function(R, filter=control$filter, 
                         partial.impute=control$partial.impute, 
                         tol=control$tol, maxiter=control$maxiter, 
                         method=control$method,
                         init=control$init,
                         control = tsgs.control()){
  
  tsgsRob <- GSE::TSGS(x=R, filter=filter,
                       partial.impute=partial.impute, tol=tol, 
                       maxiter=maxiter, method=method,
                       init=init)
  
  return(list(mu=tsgsRob@mu, sig=tsgsRob@S))
  
}


#' Title
#'
#' @param filter 
#' @param partial.impute 
#' @param tol 
#' @param maxiter 
#' @param method 
#' @param init 
#'
#' @return
#' @export
#'
#' @examples
tsgs.control <- function(filter=c("UBF-DDC","UBF","DDC","UF"),
                         partial.impute=FALSE, tol=1e-4, maxiter=150, 
                         method=c("bisquare","rocke"),
                         init=c("emve","qc","huber","imputed","emve_c")){
  filter <- match.arg(filter)
  method <- match.arg(method)
  init <- match.arg(init)
  
  return(list(filter=filter, partial.impute=partial.impute, 
              tol=tol, maxiter=as.integer(maxiter), 
              method=method,init))
  
  
}