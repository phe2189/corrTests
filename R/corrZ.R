#' Correlation Between Two Correlated Logrank Tests for Time-To-Event Endpoints
#' 
#' This function calculates the correlation between two correlated logrank tests 
#' Z_A(tau) and Z_B(tau') at analysis times tau and tau'. 
#' The correlation is due to the overlapping populations and timing.
#' #' 
#' @param e.strat Numbers of events for stratified analysis
#'                e.strat = list(AandB1, AnotB1, AandB2, BnotA2)
#'                AandB1: events for AnB at tau
#'                AnotB1: events for A and not B at tau
#'                AandB2: events for A and B at tau'
#'                BnotA2: events for B and not A at tau'
#'                
#' @param e.unstr Numbers of events for unstratified analysis
#'                e.unstrat = list(A1, B2, minAorB). 
#'                A1: events for A at tau 
#'                B2: events for B at tau'
#'                minAorB: events for AUB at min(tau, tau')
#' @param r.strat Proportions of experimental subjects for stratified 
#'                analysis in each of subgroup: 
#'                r.strat = list(AandB = 1/2, AnotB=0, BnotA=1/2).
#' @param r.unstr Proportions of experimental subjects in A, B and AandB 
#'                r.unstr = list(A=1/2, B=1/2, AandB=1/2). For simplicity, 
#'                no consideration of r change over time.
#' @param pAandB.unstr Proportion of subjects included in both A and B among all subjects in A or B. 
#' @param strat stratified analysis flag, "Y" or "N". Default, "Y". The stratified analysis 
#' means that testing HA is stratified by B, and vice versa. e.unstr, r.unstr, and 
#' gamma.unstr are required for unstratified analysis; 
#' e.strat and r.strat are required for stratified analysis. 
#' @param method Method for randomization ratio, i.e., proportion of subjects in experimental arm. 
#' "Approximation" or "Observed". "Approximation" method assumes the same randomization ratio in every stratum;
#' "Observed" method estimates the randomization ratio by the data. 
#' At study design stage without data, "Approximation" method should be used.
#' @return corr: Correlation of the two logrank test statistics
#'  
#'  @references 
#'  
#' @examples 
#'    #Example 2 in the manuscript.
#'    #1. Correlation between H1 and H2 at tau -- unstratified with approx.
#'    corrZ(e.strat = list(AandB1=120, AnotB1=100, AandB2=120, BnotA2=80),
#'    e.unstr = list(A1=220, B2=200, minAorB=300),                 
#'    r.strat = list(AandB = 1/2, AnotB=1/2, BnotA=1/2),
#'    r.unstr = list(A=1/2, B=1/2, AandB=1/2), pAandB.unstr = 0.5,
#'    strat = c("N"), method = c("Approximation"))
#'    
#'    #2. Correlation between H1 and H2  at tau -- unstratified Observed randomization ratios
#'    corrZ(e.strat = list(AandB1=120, AnotB1=100, AandB2=120, BnotA2=80),
#'    e.unstr = list(A1=220, B2=200, minAorB=300),                 
#'    r.strat = list(AandB = 1/2, AnotB=1/2, BnotA=1/2),
#'    r.unstr = list(A=1/2, B=1/2, AandB=1/2), pAandB.unstr = 0.5,
#'    strat = c("N"), method = c("Observed"))
#'    
#'    #3. Correlation between H1 and H1  at tau -- unstratified Observed randomization ratios
#'    corrZ(e.strat = list(AandB1=120, AnotB1=100, AandB2=120, BnotA2=80),
#'    e.unstr = list(A1=220, B2=220, minAorB=220),                 
#'    r.strat = list(AandB = 1/2, AnotB=1/2, BnotA=1/2),
#'    r.unstr = list(A=1/2, B=1/2, AandB=1/2), pAandB.unstr = 1,
#'    strat = c("N"), method = c("Approximation"))
#'    
#'    #4. Correlation between H1 and H1 at different times tau and tau' -- unstratified Observed randomization ratios
#'    corrZ(e.strat = list(AandB1=120, AnotB1=100, AandB2=120, BnotA2=80),
#'    e.unstr = list(A1=220, B2=280, minAorB=220),                 
#'    r.strat = list(AandB = 1/2, AnotB=1/2, BnotA=1/2),
#'    r.unstr = list(A=1/2, B=1/2, AandB=1/2), pAandB.unstr = 1,
#'    strat = c("N"), method = c("Approximation"))
#'    
#'    #5. Correlation between H1 and H2 at tau -- stratified with approx.
#'    corrZ(e.strat = list(AandB1=120, AnotB1=100, AandB2=120, BnotA2=80),
#'    e.unstr = list(A1=220, B2=200, minAorB=300),                 
#'    r.strat = list(AandB = 1/2, AnotB=1/2, BnotA=1/2),
#'    r.unstr = list(A=1/2, B=1/2, AandB=1/2), pAandB.unstr = 0.5,
#'    strat = c("Y"), method = c("Approximation"))
#'    
#'    #6. Correlation between H1 and H2 at tau -- stratified Observed randomization ratios
#'    corrZ(e.strat = list(AandB1=120, AnotB1=100, AandB2=120, BnotA2=80),
#'    e.unstr = list(A1=220, B2=200, minAorB=300),                 
#'    r.strat = list(AandB = 1/2, AnotB=1/2, BnotA=1/2),
#'    r.unstr = list(A=1/2, B=1/2, AandB=1/2), pAandB.unstr = 0.5,
#'    strat = c("Y"), method = c("Observed"))
#'    
#'    #7. Correlation between H1 and H3 at tau -- stratified Observed randomization ratios
#'    corrZ(e.strat = list(AandB1=220, AnotB1=0, AandB2=220, BnotA2=280),
#'    e.unstr = list(A1=220, B2=500, minAorB=500),                 
#'    r.strat = list(AandB = 1/2, AnotB=1/2, BnotA=1/2),
#'    r.unstr = list(A=1/2, B=1/2, AandB=1/2), pAandB.unstr = 0.5,
#'    strat = c("Y"), method = c("Observed"))    
#'    
#'    #8. Correlation between H2 and H3 at tau -- stratified Observed randomization ratios
#'    corrZ(e.strat = list(AandB1=200, AnotB1=0, AandB2=200, BnotA2=300),
#'    e.unstr = list(A1=200, B2=500, minAorB=500),                 
#'    r.strat = list(AandB = 1/2, AnotB=1/2, BnotA=1/2),
#'    r.unstr = list(A=1/2, B=1/2, AandB=1/2), pAandB.unstr = 0.5,
#'    strat = c("Y"), method = c("Observed"))    
#'    
#' @export
#' 
corrZ = function(e.strat = list(AandB1=120, AnotB1=100, AandB2=120, BnotA2=80),
                 e.unstr = list(A1=220, B2=200, minAorB=300),                 
                 r.strat = list(AandB = 1/2, AnotB=1/2, BnotA=1/2),
                 r.unstr = list(A=1/2, B=1/2, AandB=1/2),
                 pAandB.unstr = 0.5,
                 strat = c("N", "Y"),
                 method = c("Approximation", "Observed")){
  
  if (strat[1] == "N" && method[1] =="Approximation"){
    rho = pAandB.unstr * e.unstr$minAorB / sqrt(e.unstr$A1 * e.unstr$B2)
  } else if (strat[1] == "N" && method[1] =="Observed"){
    #nA * sigma2_A
    nA_S2A.unstr = r.unstr$A * (1 - r.unstr$A) * e.unstr$A1
    #nB * sigma2_B
    nB_S2B.unstr = r.unstr$B * (1 - r.unstr$B) * e.unstr$B2
    #nAB * sigma2_AB
    rr = r.unstr$A*r.unstr$B+(1-r.unstr$A-r.unstr$B)*r.unstr$AandB
    nAB_S2AB.unstr = rr*pAandB.unstr*e.unstr$minAorB
    #corr(ZA(tau), ZB(tau'))
    rho = nAB_S2AB.unstr / sqrt(nA_S2A.unstr*nB_S2B.unstr)
  } else if (strat[1] == "Y" && method[1] =="Approximation"){
    #events of A at tau
    eA1 = e.strat$AandB1 + e.strat$AnotB1
    #events of B at tau'
    eB2 = e.strat$AandB2 + e.strat$BnotA2
    #corr(ZA(tau), ZB(tau'))
    rho = min(e.strat$AandB1, e.strat$AandB2) / sqrt(eA1*eB2)
  } else if (strat[1] == "Y" && method[1] =="Observed"){
    #nAB * sigma2_AB at tau
    nAB_S2AB1 = r.strat$AandB * (1-r.strat$AandB) * e.strat$AandB1    
    #nAnotB * sigma2_AnotB at tau
    nAnotB_S2AnotB1 = r.strat$AnotB * (1 - r.strat$AnotB) * e.strat$AnotB1
    #nAB * sigma2_AB at tau'
    nAB_S2AB2 = r.strat$AandB * (1-r.strat$AandB) * e.strat$AandB2    
    #nBnotA * sigma2_BnotA at tau'
    nBnotA_S2BnotA2 = r.strat$BnotA * (1 - r.strat$BnotA) * e.strat$BnotA2

    #n * tilde Sigma2_A at tau
    n_S2A1 = nAB_S2AB1 + nAnotB_S2AnotB1
    #n * tilde Sigma2_B at tau'
    n_S2B2 = nAB_S2AB2 + nBnotA_S2BnotA2
    
    #corr(ZA(tau), ZB(tau'))
    rho = min(nAB_S2AB1, nAB_S2AB2) / sqrt(n_S2A1 * n_S2B2)
  } 
  
  return(rho)
}

