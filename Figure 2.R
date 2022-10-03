
################################
# Example 1: M+ subgroup 
# and overall population tests.
################################
#Control arm has exponential
#distribution with median 12 months in all strata. 
#1:1 randomization stratified by PD-L1+ status. 300 subjects in PD-L1+ and 
#450 total subjects in overall population. 3\% drop-off every year. Enrollment
#period is 18 months and weight 1.5.

library(corrTests)

#################
# h(t) and S(t)
#################

h0 = function(t){log(2)/12}
S0 = function(t){exp(-log(2)/12 * t)}
h0.5 = function(t){log(2)/12*0.5}
h0.65 = function(t){log(2)/12*0.65}
h0.8 = function(t){log(2)/12*0.8}
S0.5 = function(t){exp(-log(2)/12 * 0.5 * t)}
S0.65 = function(t){exp(-log(2)/12 * 0.65 * t)}
S0.8 = function(t){exp(-log(2)/12 * 0.8 * t)}

#Entry distribution: enrollment period 18 mo, acceleration weight 1.5.
Fe = function(t){(t/18)^1.5*as.numeric(t <= 18) + as.numeric(t > 18)}

#Drop-off distribution: 3\% drop-off every year.
G = function(t){1-exp(-0.03/12*t)}

#Plot correlation over time after enrollment complete
t = seq(18, 100, 5) #Analysis time, must be greater than enrollment period 18.
omega = matrix(NA, nrow=6, ncol=length(t))

for (i in 1:length(t)){
  #(1) Homogeneous regardless of PD-L1 status: HR = 0.65
  omega[1,i]=corrTime(T = t[i], n = list(AandB = 300, AnotB=0, BnotA=450), 
      r = list(AandB=1/2, AnotB =0, BnotA = 1/2), rA=1/2, rB=1/2, 
      h0=list(AandB=h0, AnotB=h0, BnotA=h0), 
      S0=list(AandB=S0, AnotB=S0, BnotA=S0),
      h1=list(AandB=h0.65, AnotB=NULL, BnotA=h0.65), 
      S1=list(AandB=S0.65, AnotB=NULL, BnotA=S0.65),
      F.entry = Fe, G.ltfu = G, strat.ana="Y")$corr
      
  omega[2,i]=corrTime(T = t[i], n = list(AandB = 300, AnotB=0, BnotA=450), 
      r = list(AandB=1/2, AnotB =0, BnotA = 1/2), rA=1/2, rB=1/2, 
      h0=list(AandB=h0, AnotB=h0, BnotA=h0), 
      S0=list(AandB=S0, AnotB=S0, BnotA=S0),
      h1=list(AandB=h0.65, AnotB=NULL, BnotA=h0.65), 
      S1=list(AandB=S0.65, AnotB=NULL, BnotA=S0.65),
      F.entry = Fe, G.ltfu = G, strat.ana="N")$corr
      
  #(2) Stronger effect in PD-L1+: HR = 0.65; PD-L1-: HR = 0.8
  omega[3,i]=corrTime(T = t[i], n = list(AandB = 300, AnotB=0, BnotA=450), 
      r = list(AandB=1/2, AnotB =0, BnotA = 1/2), rA=1/2, rB=1/2,
      h0=list(AandB=h0, AnotB=h0, BnotA=h0), 
      S0=list(AandB=S0, AnotB=S0, BnotA=S0),
      h1=list(AandB=h0.65, AnotB=NULL, BnotA=h0.8), 
      S1=list(AandB=S0.65, AnotB=NULL, BnotA=S0.8),
      F.entry = Fe, G.ltfu = G, strat.ana="Y")$corr
      
  omega[4,i]=corrTime(T = t[i], n = list(AandB = 300, AnotB=0, BnotA=450), 
      r = list(AandB=1/2, AnotB =0, BnotA = 1/2), rA=1/2, rB=1/2,
      h0=list(AandB=h0, AnotB=h0, BnotA=h0), 
      S0=list(AandB=S0, AnotB=S0, BnotA=S0),
      h1=list(AandB=h0.65, AnotB=NULL, BnotA=h0.8), 
      S1=list(AandB=S0.65, AnotB=NULL, BnotA=S0.8),
      F.entry = Fe, G.ltfu = G, strat.ana="N")$corr
      
  #(3) Weaker effect in PD-L1+: HR = 0.65; PD-L1-: HR = 0.5
  omega[5,i]=corrTime(T = t[i], n = list(AandB = 300, AnotB=0, BnotA=450), 
      r = list(AandB=1/2, AnotB =0, BnotA = 1/2), rA=1/2, rB=1/2, 
      h0=list(AandB=h0, AnotB=h0, BnotA=h0), 
      S0=list(AandB=S0, AnotB=S0, BnotA=S0),
      h1=list(AandB=h0.65, AnotB=NULL, BnotA=h0.5), 
      S1=list(AandB=S0.65, AnotB=NULL, BnotA=S0.5),
      F.entry = Fe, G.ltfu = G, strat.ana="Y")$corr
      
  omega[6,i]=corrTime(T = t[i], n = list(AandB = 300, AnotB=0, BnotA=450), 
      r = list(AandB=1/2, AnotB =0, BnotA = 1/2), rA=1/2, rB=1/2,
      h0=list(AandB=h0, AnotB=h0, BnotA=h0), 
      S0=list(AandB=S0, AnotB=S0, BnotA=S0),
      h1=list(AandB=h0.65, AnotB=NULL, BnotA=h0.5), 
      S1=list(AandB=S0.65, AnotB=NULL, BnotA=S0.5),
      F.entry = Fe, G.ltfu = G, strat.ana="N")$corr
}

#Plot the correlations vs time
tiff("Figure2M.tiff", width = 6, height = 4, units = 'in', res = 300)
par(mar=c(6, 4.1, 4.1, 2.1))
plot(t, omega[1,], type="n", ylim=range(omega),
     ylab="Correlation", cex.lab=0.8, xlab="")
title(xlab="Analysis Time (months)", line=3.2, cex.lab=0.8)
lines(t, omega[1,], lty=1, col=1, lwd=1, type="b", pch=1, cex=0.8)
lines(t, omega[2,], lty=2, col=1, lwd=1, type="b", pch=3, cex=0.8)
lines(t, omega[3,], lty=3, col=1, lwd=1, type="b", pch=2, cex=0.8)
lines(t, omega[4,], lty=4, col=1, lwd=1, type="b", pch=2, cex=0.8)
lines(t, omega[5,], lty=5, col=1, lwd=1, type="b", pch=4, cex=0.8)
lines(t, omega[6,], lty=6, col=1, lwd=1, type="b", pch=4, cex=0.8)
legend(x="bottomleft", inset=c(0, -.35), 
       c("M+/- HR 0.65/0.65: S", "M+/- HR 0.65/0.65: U"), xpd = TRUE,
       col=rep(1,2), lty=1:2, bty="n", cex=0.6, pch=c(1, 3)) 
legend(x="bottom", inset=c(0, -.35), 
       c("M+/- HR 0.65/0.8: S","M+/- HR 0.65/0.8: U"), xpd = TRUE,
       col=rep(1,2), lty=3:4, bty="n", cex=0.6, pch=c(2, 2)) 

legend(x="bottomright", inset=c(0, -.35),  
       c("M+/- HR 0.65/0.5: S","M+/- HR 0.65/0.5: U"),xpd = TRUE,
       col=rep(1,2), lty=5:6, bty="n", cex=0.6, pch=c(4, 4))
dev.off()

