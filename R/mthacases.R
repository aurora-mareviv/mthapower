#' Sample size calculations - mtDNA haplogroups

#' Determine the minimum number of cases (Ncmin), required to detect: either a change from p0 to p1, or a given OR, with a predefined confidence interval, in a study with Nh haplogroups.
#' Note: I assume that case-control equations are valid for a cohort.
#' This function may not be generalizable for all studies.

#' USAGE:
#' @param p0 the frequency of the haplogroup in the control population, the controls among exposed. It depends on haplogroup baseline frequency.
#' @param Nh number of categories for haplogroups. Usually 10 haplogroups plus one category for rare haplogrs.: Nh <- 11
#' @param OR.cas.ctrl (p1 / (1-p1)) / (p0 / (1-p0)) the OR you want to detect with your data. It can be either a single value, or a sequence: OR.cas.ctrl <- 2; OR.cas.ctrl <- seq(1.25,3 by=0.5)
#' @param power the power I want to detect any OR in my study.
#' @param sig.level error alpha accepted. Can take 3 possible values: 0.05, 0.01 and 0.001 (see [Table 2] of Samuels et al.)
#' @param cases.min number of cases or controls that I need to recruit.
#'
#' Gives the result in a data frame, easy to print in a plot.
#' @return None
#'
#' @examples
#' mydata <- mthacases(p0=0.445, Nh=11, OR.cas.ctrl=c(2), power=80, sig.level=0.05) # Baudouin
#' mydata <- mthacases(p0=0.445, Nh=11, OR.cas.ctrl=c(1.25,1.5,1.75,2,2.25,2.5,2.75,3), power=80, sig.level=0.05)
#' mydata <- mydata[c(2,6)]
#' mydata
#' plot(mydata)
#'
#' @export
mthacases <-
  function(p0=p0, Nh=Nh, OR.cas.ctrl=OR.cas.ctrl, power=power, sig.level=sig.level)
  { ... }


#'
#' @author Aurora Baluja, \email{mariauror@gmail.com}
#' @references \url{http://www.cell.com/ajhg/fulltext/S0002-9297(07)63709-4}
#' @references \url{https://github.com/aurora-mareviv/mthapower}

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
