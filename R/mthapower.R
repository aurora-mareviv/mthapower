mthapower <-
function(n.cases=ncases, p0=p0, Nh=Nh, OR.cas.ctrl=OR.cas.ctrl, sig.level=sig.level){
      ncases <- n.cases
      p0 <- p0
      Nh <- Nh
      OR.cas.ctrl <- OR.cas.ctrl
      sig.level <- sig.level
      # Parameters related to sig.level, from [Table 2] of Samuels et al.
      # For 90% power and alpha = .05, Nscaled = 8.5
      if (sig.level == 0.05){
A <- -28 # Parameter A for alpha=.05
x0 <- 2.6 # Parameter x0 for alpha=.05
d <- 2.4 # Parameter d for alpha=.05
      }
      if (sig.level == 0.01){
A <- -13 # Parameter A for alpha=.01
x0 <- 5 # Parameter x0 for alpha=.01
d <- 2.5 # Parameter d for alpha=.01
      }
      if (sig.level == 0.001){
A <- -7 # Parameter A for alpha=.001
x0 <- 7.4 # Parameter x0 for alpha=.001
d <- 2.8 # Parameter d for alpha=.001
      }
      pow.ns <- NULL # initialize vector
      out.pow <- NULL # initialize vector
for(n.cases in n.cases){
    OR.ctrl.cas <- 1 / OR.cas.ctrl # 1. CALCULATE P1 FROM A PREDEFINED P0, AND A DESIRED OR
    OR <- OR.ctrl.cas  
    bracket.pw <- p0 / (OR - OR*p0) # obtained after isolating p1 in OR equation [3].
    p1 <- bracket.pw / (1 + bracket.pw)
    Nh037 <- Nh^0.37 # 2. CALCULATE NSCALED
    num.n <- ncases*((p1-p0)^2)
    den.n <- (p1*(1-p1) + p0*(1-p0))*Nh037
    Nscaled.data <- num.n/den.n
      num.power <- A - 100 # 3. CALCULATE POWER
      den.power <- 1 + exp((Nscaled.data - x0)/d)
      power <- 100 + (num.power/den.power) # The power I have to detect a given OR with my data, at a given alpha
}

      out.pow <- data.frame(Nh,ncases, p0, p1, OR.ctrl.cas, OR.cas.ctrl, power, sig.level)
      out.pow 
      round(out.pow,3) # Results given with 3 decimals
}
