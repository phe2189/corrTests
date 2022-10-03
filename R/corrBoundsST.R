#' Rejection Boundaries Of 2+ Correlated Tests for Time-To-Event Endpoints at a Single Time
#' 
#' This function calculates the rejection boundary in p value 
#' (significance level) and z value for each of multiple correlated logrank test.
#' The computation is based on mvtnorm R package.
#' 
#' @param alpha  Overall one-sided type I error for the correlated hypothesis tests. Default 0.025.
#' @param w Weights for alpha allocation for the correlated tests. 
#'          The sum of w must be 1.
#' @param corr Correlation matrix for the correlated tests.          
#' @param eps A vector of efficiency factors for testing all correlated hypotheses.
#' eps is required when method = "Customized Allocation". Must have at least one
#' component as NA. For example, eps = c(1, 1, NA) means allocating all gained
#' efficiency to H3. eps = c(NA, NA, NA) means equally allocating the gained 
#' efficiency in all tests. For all j, eps[j] >= 1.
#'
#' @return A data frame including variables
#'  \itemize{
#'     \item p0 Rejection boundary in P without consideration of correlation
#'     \item z0 Rejection boundary in Z without consideration of correlation 
#'     \item p  Rejection boundary in P with consideration of correlation
#'     \item z  Rejection boundary in Z with consideration of correlation
#'     \item eps  Efficiency factor
#'     \item max.eps   Maximum possible efficiency factor for each test
#'  }
#'  @references 
#'  
#' @examples 
#' #Example 2. Two subgroups and overall population.
#' corr = matrix(1, nrow=3, ncol=3)
#' corr[1, 2] = corr[2, 1] = corrZ(e.strat = list(AandB1=120, AnotB1=100, AandB2=120, BnotA2=80),
#'    e.unstr = list(A1=220, B2=200, minAorB=300),                 
#'    r.strat = list(AandB = 1/2, AnotB=1/2, BnotA=1/2),
#'    r.unstr = list(A=1/2, B=1/2, AandB=1/2), pAandB.unstr = 0.5,
#'    strat = c("Y"), method = c("Observed"))
#'    
#' corr[1, 3] = corr[3, 1] = corrZ(e.strat = list(AandB1=220, AnotB1=0, AandB2=220, BnotA2=280),
#'    e.unstr = list(A1=220, B2=500, minAorB=500),                 
#'    r.strat = list(AandB = 1/2, AnotB=1/2, BnotA=1/2),
#'    r.unstr = list(A=1/2, B=1/2, AandB=1/2), pAandB.unstr = 0.5,
#'    strat = c("Y"), method = c("Observed"))    
#' corr[2, 3] = corr[3, 2] = corrZ(e.strat = list(AandB1=200, AnotB1=0, AandB2=200, BnotA2=300),
#'    e.unstr = list(A1=200, B2=500, minAorB=500),                 
#'    r.strat = list(AandB = 1/2, AnotB=1/2, BnotA=1/2),
#'    r.unstr = list(A=1/2, B=1/2, AandB=1/2), pAandB.unstr = 0.5,
#'    strat = c("Y"), method = c("Observed"))    
#'    
#' #Strategy 1: equal allocation of efficiency to all tests 
#' corrBoundsST(alpha = 0.025, w = c(1/4, 1/4, 1/2), eps = c(NA,NA,NA), corr=corr)
#' 
#' #Strategy 2: equal allocation of efficiency to H1 and H2, no efficiency to H3 
#' corrBoundsST(alpha = 0.025, w = c(1/4, 1/4, 1/2), eps = c(NA,NA,1), corr=corr)
#' 
#' #Strategy 3: maximize power to H3 
#' corrBoundsST(alpha = 0.025, w = c(1/4, 1/4, 1/2), eps = c(1,1,NA), corr=corr)
#' 
#' #Strategy 4: maximize power to H1 
#' corrBoundsST(alpha = 0.025, w = c(1/4, 1/4, 1/2), eps = c(NA,1,1), corr=corr)
#'    
#' #Strategy 5: maximize power to H2 
#' corrBoundsST(alpha = 0.025, w = c(1/4, 1/4, 1/2), eps = c(1,NA,1), corr=corr)
#'    
#' @export
#' 
corrBoundsST = function(alpha = 0.025, w = c(1/4, 1/4, 1/2), 
                        eps = c(NA,NA,NA), corr=diag(3)){


  #Boundaries without consideration of correlations (P and Z)
  p0 = w * alpha
  z0 = qnorm(1 - p0)

  #Dimension
  D = length(w)
  
  #Find the max epsilons, i.e., setting all eps = 1 except for one test
  max.eps = rep(NA, D)
  for (i in 1:D){
    f.eps.max = function(x){
      if (i == 1) {
        u = qnorm(c(1-x*p0[1], 1-p0[2:D]))
      } else if (i > 1 && i < D){
        u = qnorm(c(1-p0[1:(i-1)], 1-x*p0[i], 1-p0[(i+1):D]))
      } else {
        u = qnorm(c(1-p0[1:(D-1)], 1-x*p0[D]))
      }
      I = mvtnorm::pmvnorm(lower = rep(-Inf, D), upper = u, 
                           corr = corr, abseps = 1e-8, maxpts=100000)[1]
      return(1 - I - alpha)
    }
    max.eps[i] = uniroot(f=f.eps.max, interval=c(1, 10), tol=1e-8)$root
  }
  
  #For customized efficiency strategy, solve for the corresponding eps
  f.eps = function(x){
    for (i in 1:D) if(is.na(eps[i])){eps[i] = x}
    
    I = mvtnorm::pmvnorm(lower = rep(-Inf, D), upper = qnorm(1-eps*p0), 
                         corr = corr, abseps = 1e-8, maxpts=100000)[1]
    return(1 - I - alpha)
  }
  eps.x = uniroot(f=f.eps, interval=c(1, 10), tol=1e-8)$root
  for (i in 1:D) if(is.na(eps[i])){eps[i] = eps.x}

  p = eps * p0
  z = qnorm(1-p)
  
  out = data.frame(cbind(p0, z0, p, z, eps, max.eps))
  return(out)
}
