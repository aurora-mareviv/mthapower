mthacases <- 
function(p0=p0, Nh=Nh, OR.cas.ctrl=OR.cas.ctrl, power=power, sig.level=sig.level){
      power <- power
      p0 <- p0
      Nh <- Nh
      ORcas.ctrl <- OR.cas.ctrl
      sig.level <- sig.level
      # Parameters related to sig.level, from [Table 2] of Samuels et al.
      # For 90% power and alpha = .05, Nscaled = 8.5
      if (sig.level == 0.05){
A <- -28 		# Parameter A for alpha=.05
x0 <- 2.6 		# Parameter x0 for alpha=.05
d <- 2.4 		# Parameter d for alpha=.05
      }
      if (sig.level == 0.01){
A <- -13 		# Parameter A for alpha=.01
x0 <- 5 		# Parameter x0 for alpha=.01
d <- 2.5 		# Parameter d for alpha=.01
      }
      if (sig.level == 0.001){
A <- -7 		# Parameter A for alpha=.001
x0 <- 7.4 		# Parameter x0 for alpha=.001
d <- 2.8 		# Parameter d for alpha=.001
      }
      cases.min <- NULL # initialize vector
      out.cas <- NULL # initialize vector
for(OR.cas.ctrl in OR.cas.ctrl){
    logparen <- (power-A)/(100-power) # 1. CALCULATES Nscaled FROM FIXED PARAMETERS AND DESIRED POWER
    Nscaled <- x0+d*log(logparen)
    ORctrl.cas <- 1 / ORcas.ctrl # 2. CALCULATES P1 FROM A PREDEFINED P0, AND A DESIRED OR
    OR <- ORctrl.cas  
    bracket.pw <- p0 / (OR - OR*p0) # obtained after isolating p1 in OR equation [3].
    p1 <- bracket.pw / (1 + bracket.pw)
      Nh037 <- Nh^0.37
      nump1 <- p1*(1-p1)
      nump0 <- p0*(1-p0)
      num <- nump1+nump0
      den <- (p1-p0)^2
      paren <- num/den
      cases.min <- Nscaled*Nh037*paren # 3. CALCULATES NCMIN	
} 	# Number of cases or controls required to detect a given OR at a desired power.
      out.cas <- data.frame(Nh, cases.min, p0, p1, ORctrl.cas, ORcas.ctrl, power, sig.level)
      out.cas 
      round(out.cas,3) # Results given with 3 decimals
}
