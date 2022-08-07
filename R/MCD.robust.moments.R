#' @title Compute moments using MCD robust method
#' 
#' @description 
#' This is a function that makes use of covMcd function from robustbase package. 
#' For more details of the parameters and function please refer to the manual of 
#' robustbase on cran.
#'
#' @param R xts object of asset returns
#' @param alpha numeric parameter controlling the size of the subsets over 
#'              which the determinant is minimized. Allowed values are between 
#'              0.5 and 1 and the default is 0.5.
#' @param nsamp number of subsets used for initial estimates or "best", "exact", 
#'              or "deterministic". Default is nsamp = 500. For nsamp = "best" 
#'              exhaustive enumeration is done, as long as the number of trials 
#'              does not exceed 100'000, which is the value of nlarge. For "exact", 
#'              exhaustive enumeration will be attempted however many samples are needed.
#'              In this case a warning message may be displayed saying that 
#'              the computation can take a very long time.
#'              For "deterministic", the deterministic MCD is computed; 
#'              as proposed by Hubert et al. (2012) it starts from the h most 
#'              central observations of six (deterministic) estimators.
#' @param nmini,kimi for n >= 2*n0, n0 := nmini, the algorithm splits the data 
#'                    into maximally kmini (by default 5) subsets, of size approximately, 
#'                    but at least nmini. When nmini*kmini < n, the initial search 
#'                    uses only a subsample of size nmini*kmini. The original algorithm 
#'                    had nmini = 300 and kmini = 5 hard coded.
#' @param scalefn function to compute a robust scale estimate or character string 
#'                specifying a rule determining such a function for the deterministic MCD.
#'                The default is "hrv2012". Another option value is "v2014".
#' @param maxcsteps maximal number of concentration steps in the deterministic MCD
#' @param initHsets NULL or a Kxh integer matrix of initial subsets of observations of 
#'                  size h.
#' @param seed initial seed for random generator
#' @param tolSolve numeric tolerance to be used for inversion of the covariance matrix
#' @param wgtFUN a character string or function, specifying how the weights for 
#'               the reweighting step should be computed. Default is "01.originalz".
#' @param use.correction whether to use finite sample correction factors; defaults to TRUE.
#' @param control a list with estimation options - this includes those above provided 
#'                in the function specification, see mcd.control() function for 
#'                the defaults.
#'
#' @return estimators of first and second moments
#' @export
#'
#' @examples
MCD.robust.moment <- function(R, alpha=mcd.control$alpha, nsamp=mcd.control$nsamp, 
                              nmini=mcd.control$nmini, kmini=mcd.control$kmini,
                              scalefn=mcd.control$scalefn, maxcsteps=mcd.control$maxcsteps, 
                              initHsets = NULL, 
                              seed=mcd.control$seed, tolSolve=mcd.control$tolSolve, 
                              wgtFUN=mcd.control$wgtFUN, 
                              use.correction=mcd.control$use.correction,
                              control=mcd.control()){
  
  robustMCD <- robustbase::covMcd(x=R,
         alpha=alpha, nsamp=nsamp, nmini=nmini, kmini=kmini,
         scalefn=scalefn, maxcsteps=maxcsteps,
         initHsets = NULL, seed=seed, tolSolve=tolSolve,
         wgtFUN=wgtFUN, control=control, use.correction=use.correction)
  
  return(list(rbMCD.Mu = robustMCD$center, rbMCD.Sig = robustMCD$cov))
}

#' @title 
#' Control settings for MCD robust moments
#' 
#' @description 
#' Auxiliary function for passing the estimation options as parameters 
#' to the estimation function MCD.robust.moment
#' 
#' @param alpha numeric parameter controlling the size of the subsets over 
#'              which the determinant is minimized. Allowed values are between 
#'              0.5 and 1 and the default is 0.5.
#' @param nsamp number of subsets used for initial estimates or "best", "exact", 
#'              or "deterministic". Default is nsamp = 500. For nsamp = "best" 
#'              exhaustive enumeration is done, as long as the number of trials 
#'              does not exceed 100'000, which is the value of nlarge. For "exact", 
#'              exhaustive enumeration will be attempted however many samples are needed.
#'              In this case a warning message may be displayed saying that 
#'              the computation can take a very long time.
#'              For "deterministic", the deterministic MCD is computed; 
#'              as proposed by Hubert et al. (2012) it starts from the h most 
#'              central observations of six (deterministic) estimators.
#' @param nmini,kimi for n >= 2*n0, n0 := nmini, the algorithm splits the data 
#'                    into maximally kmini (by default 5) subsets, of size approximately, 
#'                    but at least nmini. When nmini*kmini < n, the initial search 
#'                    uses only a subsample of size nmini*kmini. The original algorithm 
#'                    had nmini = 300 and kmini = 5 hard coded.
#' @param scalefn function to compute a robust scale estimate or character string 
#'                specifying a rule determining such a function for the deterministic MCD.
#'                The default is "hrv2012". Another option value is "v2014".
#' @param maxcsteps maximal number of concentration steps in the deterministic MCD
#' @param seed initial seed for random generator
#' @param tolSolve numeric tolerance to be used for inversion of the covariance matrix
#' @param wgtFUN a character string or function, specifying how the weights for 
#'               the reweighting step should be computed. Default is "01.originalz".
#' @param beta a quantile, experimentally used for some of the prespecified wgtFUNs. For our 
#'             MCD method, the default is 0.975.
#' @param use.correction whether to use finite sample correction factors; defaults to TRUE.
#' @return a list of passed parameters
#' @export
#'
#' @examples
mcd.control <-
  function(alpha = 1/2,
           nsamp = 500, nmini = 300, kmini = 5,
           scalefn = "hrv2012", maxcsteps = 200,
           seed = NULL, tolSolve = 1e-14,
           wgtFUN = "01.original", beta,
           use.correction=TRUE
           )
  {
    if(missing(beta) || !is.numeric(beta))
      beta <- 0.975
    list(alpha=alpha, nsamp=nsamp, nmini=as.integer(nmini), kmini=as.integer(kmini),
         seed = as.integer(seed),
         tolSolve=tolSolve, scalefn=scalefn, maxcsteps=as.integer(maxcsteps),
         wgtFUN=wgtFUN, beta=beta,
         use.correction=use.correction)
  }
