
rm(list = ls())
set.seed(1)

setwd("") ## Insert your working directory where you need to have a "data" folder
datadir <- "data/" ## Here the simulation output will be saved

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Simulation parameters
B <- 2000 # Number of iterations
N_boot <- 500 # Number of iterations for bootstrapping   
LEVEL <- .95 # confidence level
PROBS <- c( (1-LEVEL)/2,1-(1-LEVEL)/2 )
PARlist <- list(
  n_vec = c(10,50,100,300,500,1000),
  delta_vec = c(0, 2),
  sigma_vec = c(1, 5),
  alpha_vec = c(0, 10)
)

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
### Functions
snpar <- function(xi=0,omega=1,alpha=0) {
  delta <- alpha/sqrt(1+alpha^2)
  mu <- xi + omega * delta * sqrt( 2/pi )
  sigma2 <- omega^2 * ( 1 - (2*delta^2)/pi )
  return(list(mu = mu, sigma = sqrt(sigma2)))
}

##
sninvpar <- function( mu=0, sigma=1, xi=NULL, omega=NULL, alpha=0 ) {
  if (is.null(omega)) {
    delta <- alpha/sqrt(1+alpha^2)
    omega2 <- sigma^2 / ( 1 - (2*delta^2) / pi )
    omega <- sqrt( omega2 )
  }
  if (is.null(xi)) {
    delta <- alpha/sqrt(1+alpha^2)
    xi <- mu - omega * delta * sqrt( 2/pi )
  }
  return( list( xi = xi, omega = omega, alpha = alpha ) )
}

##
min_dskew_normal <- function( x = seq( -5, 5, by = .01 ), xi = 0, omega = 1, alpha = 0, 
                              return.all = FALSE ) {
  
  require( sn )
  y1 <- dsn( x, xi = 0, omega = 1, alpha = 0 )
  y2 <- dsn( x, alpha = alpha, xi = xi, omega = omega )
  dy <- ifelse( y1 < y2, y1, y2 )
  gData <- data.frame( x, y1, y2, dy )  
  
  if (return.all) {
    return( list( gData = gData ) )
  } else {
    return( dy )  
  }
}

##
true_cohen <- function( mu1, mu2, sigma1, sigma2 ) {
  pool.s <- sqrt( (sigma1^2 + sigma2^2) / 2 )
  (mu1-mu2) / pool.s 
}

##
relative_mean_bias <- function(theta, true_theta, eps = 1e-6, use_abs_for_zero = TRUE) {
  # theta: estimated values
  # true_theta: true values of effect size
  # eps: small value to avoid deviding by zero
  # use_abs_for_zero: if TRUE, when true_theta=0 returns absolute bias
  
  bias <- theta - true_theta
  
  rmb <- ifelse(
    true_theta == 0,
    if (use_abs_for_zero) bias else bias / eps,   # deals with cases of zero
    bias / true_theta
  )
  
  return(rmb)
}

# Function to calculate the CLES index and it's confidence interval
CLES_IC <- function(x1, x2, alpha = .9 ) {
  # Calculates CLES index
  n1 <- length(x1)
  n2 <- length(x2)
  
  # Calculates number of comparisons
  confronti <- sum(outer(x1, x2, ">"))
  
  # Calculates CLES
  CLES <- confronti / (n1 * n2)
  
  # Calculates variance
  var_CLES <- (CLES * (1 - CLES)) / (n1 * n2)
  
  # Calculates the confidence interval
  Z <-  qnorm(1 - (1 - alpha) / 2) 
  errore_standard <- sqrt(var_CLES)
  
  IC_lower <- CLES - Z * errore_standard
  IC_upper <- CLES + Z * errore_standard
  
  # <returns the results
  return(list(CLES = CLES, IC = c(IC_lower, IC_upper)))
}

true_CLES_OV <- function( true_d ) {
  true_CLES <- pnorm( abs(true_d) / sqrt(2) )
  true_OV <- 2*pnorm( (abs(true_d)*(-1)) / 2 )
  
  return(list(true_CLES = true_CLES, true_OV = true_OV ))
} 


## 
core_sim <- function( DESIGN, B = 10 ) {  
  
  require(overlapping)
  require(sn)
  require(effectsize)
  # set simulation parameters
  n <- DESIGN$n
  delta <- DESIGN$xi 
  alpha <- DESIGN$alpha 
  omega <- DESIGN$omega
  condition <- DESIGN$K
  
  # calculates mu and sigma
  MUSI <- snpar(xi = delta, omega = omega, alpha = alpha)
  
  cat( paste0( "condition ",condition, ": n = ",n,"; mu = ",round(MUSI$mu,1), 
               "; sigma = ",round(MUSI$sigma,2),"; alpha = ",
               round(alpha,2) ),"\n")
  
  # true population overlap
  true_overlap <- integrate( min_dskew_normal, -Inf, Inf, 
                                xi = delta, alpha = alpha, 
                                omega = omega )$value
  
  # true population Cohen's d
  true_d <- true_cohen( mu1 = 0, mu2 = MUSI$mu, sigma1 = 1, 
                        sigma2 = MUSI$sigma )
  
  # true CLES and OV
  true_cles_ov <- true_CLES_OV( true_d )
  
  # data simulation
  t( sapply(1:B, function(b){
    
    y1 <- rsn(n,xi=0, omega=1, alpha=0)
    y2 <- rsn(n,xi=delta, omega=omega,alpha=alpha)
    
    mx1 <- mean(y1); sx1 <- sd(y1)
    mx2 <- mean(y2); sx2 <- sd(y2)
    
    ETA1 <- overlap(list(y1,y2))$OV
    ETA2 <- overlap(list(y1,y2), type = "2" )$OV
    COHEN <- cohens_d(y1,y2, ci = LEVEL )
    OV <- p_overlap( y1, y2, ci = LEVEL )
    
    ETA_boot <- boot.overlap( list(y1,y2), B = N_boot )
    ETA1_CI <- quantile( ETA_boot$OVboot_dist, probs = PROBS )
    ETA_boot <- boot.overlap( list(y1,y2), B = N_boot, type = "2" )
    ETA2_CI <- quantile( ETA_boot$OVboot_dist, probs = PROBS )
    cles_ic <- CLES_IC( y1, y2, alpha = LEVEL )
    
    output <- c( mx1, sx1, mx2, sx2, ETA1, ETA2, n, 
                 delta, alpha, omega, 
                 true_overlap, true_d, COHEN$Cohens_d, OV$Overlap, 
                 COHEN$CI_low, COHEN$CI_high, ETA1_CI[1], ETA1_CI[2], 
                 ETA2_CI[1], ETA2_CI[2], 
                 cles_ic$CLES, cles_ic$IC[1], cles_ic$IC[2],
                 OV$CI_low, OV$CI_high)
    names(output) <- c( "mx1", "sx1", "mx2", "sx2", "eta1", "eta2", "n",
                        "xi","alpha","omega", 
                        "true_overlap", "true_d","Cohens_d","OV",
                        "Cohen_low","Cohen_high",
                        "eta1_low", "eta1_high","eta2_low", "eta2_high",
                        "Cles","Cles_low","Cles_high","OV_low","OV_high")
    return(output)  
  }) )
}

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Prepare Design
MUSI <- with( PARlist, expand.grid(mu=delta_vec, sigma=sigma_vec, alpha = alpha_vec) )

DESIGN <- NULL 
for (i in 1:nrow(MUSI)) {
  DESIGN <- rbind( DESIGN, 
                   unlist( with( MUSI, sninvpar(mu[i],sigma[i],alpha = alpha[i]) )) )
}
DD <- data.frame(DESIGN)

DESIGN <- NULL
for (n in PARlist$n_vec) {
  DD$n <- n
  DESIGN <- rbind(DESIGN,DD)
}
DESIGN$K <- 1:nrow(DESIGN)
DESIGN <- apply(DESIGN, 1, as.list)

## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Simulation

library( parallel )
INIZIO <- Sys.time()
OUT <- mclapply(DESIGN, core_sim, B=B, mc.cores=detectCores())
FINE <- Sys.time()

ALL <- data.frame( do.call(rbind, OUT) )
ALL$mu <- with(ALL,snpar(xi,omega,alpha))$mu
ALL$sigma <- with(ALL,snpar(xi,omega,alpha))$sigma
ALL$bias_d <- with(ALL, Cohens_d - true_d )
ALL$bias_eta1 <- with(ALL, eta1 - true_overlap ) 
ALL$bias_eta2 <- with(ALL, eta2 - true_overlap ) 
ALL$rbias_eta1 <- with( ALL, (eta1 - true_overlap)/true_overlap )
ALL$rbias_eta2 <- with( ALL, (eta2 - true_overlap)/true_overlap )
ALL$rbias_d <- with( ALL, 
                     relative_mean_bias(Cohens_d, true_d ))

save(ALL, INIZIO, B, PARlist, FINE,  file = paste0(datadir, "exp2601.rda") )  

