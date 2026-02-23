
rm(list = ls())
set.seed(1)

#setwd("/mnt/data/users/massimiliano.pastore")
setwd("/home/cox/MEGA/lavori/overlapping/")
datadir <- "data/"


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
  # theta: valori stimati
  # true_theta: valori veri
  # eps: piccolo valore per evitare divisione per zero
  # use_abs_for_zero: se TRUE, quando true_theta=0 ritorna il bias assoluto
  
  bias <- theta - true_theta
  
  rmb <- ifelse(
    true_theta == 0,
    if (use_abs_for_zero) bias else bias / eps,   # gestione dei casi con 0
    bias / true_theta
  )
  
  return(rmb)
}

# Funzione per calcolare l'indice CLES e il suo intervallo di confidenza
CLES_IC <- function(x1, x2, alpha = .95 ) {
  # Calcola l'indice CLES
  n1 <- length(x1)
  n2 <- length(x2)
  
  # Calcola il numero di confronti
  confronti <- sum(outer(x1, x2, ">"))
  
  # Calcola CLES
  CLES <- confronti / (n1 * n2)
  
  # Calcola la varianza
  var_CLES <- (CLES * (1 - CLES)) / (n1 * n2)
  
  # Calcola l'intervallo di confidenza
  Z <-  qnorm(1 - (1 - alpha) / 2) 
  errore_standard <- sqrt(var_CLES)
  
  IC_lower <- CLES - Z * errore_standard
  IC_upper <- CLES + Z * errore_standard
  
  # Restituisce i risultati
  return(list(CLES = CLES, IC = c(IC_lower, IC_upper)))
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
  
  # calcolo mu e sigma
  MUSI <- snpar(xi = delta, omega = omega, alpha = alpha)
  
  cat( paste0( "n = ",n,"; mu = ",round(MUSI$mu,1), 
               "; sigma = ",round(MUSI$sigma,2),"; alpha = ",
               round(alpha,2) ),"\n")
  
  # true population overlap
  true_overlap <- integrate( min_dskew_normal, -Inf, Inf, xi = delta, alpha = alpha, 
                             omega = omega )$value
  
  # true population Cohen's d
  true_d <- true_cohen( mu1 = 0, mu2 = MUSI$mu, sigma1 = 1, 
                        sigma2 = MUSI$sigma )
  
  # data simulation
  t( sapply(1:B, function(b){
    
    y1 <- rsn(n,xi=0, omega=1, alpha=0)
    y2 <- rsn(n,xi=delta, omega=omega,alpha=alpha)
    
    mx1 <- mean(y1); sx1 <- sd(y1)
    mx2 <- mean(y2); sx2 <- sd(y2)
    
    ETA1 <- overlap(list(y1,y2))$OV
    ETA2 <- overlap(list(y1,y2), type = "2" )$OV
    COHEN <- cohens_d(y1,y2, ci = .95)
    CLES <- d_to_overlap(COHEN$Cohens_d)
    COHEN_CI_low <- COHEN$CI_low
    COHEN_CI_high <- COHEN$CI_high
    ETA_boot <- boot.overlap( list(y1,y2), B = 200 )
    ETA_CI <- quantile( ETA_boot$OVboot_dist, probs = c(.025,.975))
    
    output <- c( mx1, sx1, mx2, sx2, ETA1, ETA2, n, 
                 delta, alpha, omega, 
                 true_overlap, true_d, COHEN$Cohens_d, CLES, 
                 COHEN_CI_low, COHEN_CI_high, ETA_CI[1], ETA_CI[2])
    names(output) <- c( "mx1", "sx1", "mx2", "sx2", "eta1", "eta2", "n",
                        "xi","alpha","omega", 
                        "true_overlap", "true_d","Cohens_d","cles",
                        "Cohen_low","Cohen_high",
                        "eta_low", "eta_high")
    return(output)  
  }) )
}


## +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
## Parameters
B <- 500
PARlist <- list(
  n_vec = c(10,50,100,300,500,1000),
  delta_vec = c(0, 2),
  sigma_vec = c(1, 5),
  alpha_vec = c(0, 10)
)

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

save(ALL, INIZIO, FINE,  file = paste0(datadir, "ALL_95.rda") )  

