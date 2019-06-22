#' Power calculations - mtDNA haplogroups
#'
#' For a given study size, determine the minimum effect size that can be detected with the desired power and significance level, in a study with \code{Nh} haplogroups.
#' Note: I assume that case-control equations are valid for cohorts with a balanced number of cases and controls.
#' This function may not be generalizable for all studies involving mtDNA haplogroups.
#
#' @param n.cases number of cases or controls from the study. It can be either a single value, or a sequence: \code{n.cases <- 300}; \code{n.cases <- seq(50,500 by=10)}.
#' @param p0 the frequency of the haplogroup in the control population. It depends on haplogroup baseline frequency.
#' @param Nh number of categories for haplogroups. Usually 10 haplogroups plus one category for rare haplogroups: \code{Nh <- 11}.
#' @param OR.cas.ctrl (p1 / (1-p1)) / (p0 / (1-p0)) the OR you want to detect with your data.
#' @param sig.level the alpha error accepted. Can take 3 possible values: \code{0.05}, \code{0.01} and \code{0.001} (see [Table 2] of Samuels et al).
#'
#'
#' @author Author and maintainer: Aurora Baluja. Email: \email{mariauror@gmail.com}
#' @references 1. DC Samuels, AD Carothers, R Horton, PF Chinnery. The Power to Detect Disease Associations with Mitochondrial DNA Haplogroups. AJHG, 2006. 78(4):713-720. \href{https://www.ncbi.nlm.nih.gov/pmc/PMC1424681}{DOI:10.1086/502682}.
#' @references 2. Source code: \href{https://github.com/aurora-mareviv/mthapower}{github.com/aurora-mareviv/mthapower}.
#' @references 3. Shiny app: \href{https://aurora.shinyapps.io/mtDNA_power_calc}{aurora.shinyapps.io/mtDNA_power_calc}.
#'
#'
#' @return Calculates power given the number of cases and other parameters. The output is an object of class \code{data.frame}, ready to plot.
#'
#' @examples
#' # Example 1:
#' pow <- mthapower(n.cases=203, p0=0.443, Nh=13, OR.cas.ctrl=2.33, sig.level=0.05)
#'
#' # Example 2:
#' # Create data frames
#' pow.H150 <- mthapower(n.cases=seq(50,1000,by=50), p0=0.433, Nh=11,
#'                       OR.cas.ctrl=1.5, sig.level=0.05)
#' pow.H175 <- mthapower(n.cases=seq(50,1000,by=50), p0=0.433, Nh=11,
#'                       OR.cas.ctrl=1.75, sig.level=0.05)
#' pow.H200 <- mthapower(n.cases=seq(50,1000,by=50), p0=0.433, Nh=11,
#'                       OR.cas.ctrl=2, sig.level=0.05)
#' pow.H250 <- mthapower(n.cases=seq(50,1000,by=50), p0=0.433, Nh=11,
#'                       OR.cas.ctrl=2.5, sig.level=0.05)
#'
#' # Bind the three data frames:
#' bindata <- rbind(pow.H150,pow.H175,pow.H200,pow.H250)
#' # Adds column OR to binded data frame:
#' bindata$OR <- rep(factor(c(1.50,1.75,2,2.5)),
#'               times = c(nrow(pow.H150),
#'                         nrow(pow.H175),
#'                         nrow(pow.H200),
#'                         nrow(pow.H250)))
#' # Create plot:
#' # install.packages("car")
#' library(car)
#' scatterplot(power~ncases | OR, regLine=FALSE,
#'             smooth=FALSE,
#'             boxplots=FALSE,  by.groups=TRUE,
#'             data=bindata)
#'
#' @export
mthapower <-
  function(n.cases=ncases, p0=p0, Nh=Nh, OR.cas.ctrl=OR.cas.ctrl, sig.level=sig.level)
  { ... }


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
