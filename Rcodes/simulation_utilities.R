## Utilities for the simulation

# +++++++++++++++++++++++++++++++++++++++++++++
#' @name mindskew_normal_brms
#' @description Calcola il minimo tra due densità Skew-Normal
#' @note La prima densità ha i parametri fissati \code{xi = 0}, 
#' \code{omega = 1}, \code{alpha = 0}, per cui di fatto è normale standard
#' @param x = x vector
#' @param xi = location parameter
#' @param omega = scale parameter
#' @param alpha = skewness parameter
#' @param plot = logical, if \code{TRUE} produce la rappresentazione grafica delle densità 
#' e dell'area di sovrapposizione
#' @param return.all = logical, if \code{TRUE} restituisce il data set completo delle 
#' densità delle due distribuzioni
min_dskew_normal_brms <- function( x = seq( -5, 5, by = .01 ), xi = 0, omega = 1, alpha = 0, 
                              plot = FALSE, return.all = FALSE ) {
  
  require( brms )
  require( ggplot2 )
  y1 <- dskew_normal( x, xi = 0, omega = 1, alpha = 0 )
  y2 <- dskew_normal( x, alpha = alpha, xi = xi, omega = omega )
  dy <- ifelse( y1 < y2, y1, y2 )
  gData <- data.frame( x, y1, y2, dy )  
  
  OVplot <- ggplot(gData, aes(x,dy) ) + theme_bw() +
    geom_ribbon( aes(x, ymin = 0, ymax = dy), fill = "orange") + 
    geom_line() + geom_line( aes(x,y1), color="red" ) + geom_line( aes(x,y2), color="blue" )
  
  if (plot) {  
    print( OVplot  )
  }
  
  if (return.all) {
    return( list( gData = gData, OVplot = OVplot ) )
  } else {
    return( dy )  
  }
}
# END

# +++++++++++++++++++++++++++++++++++++++++++++
#' @name min_dskew_normal
#' @description Calcola il minimo tra due densità Skew-Normal
#' @note La prima densità ha i parametri fissati \code{xi = 0}, 
#' \code{omega = 1}, \code{alpha = 0}, per cui di fatto è normale standard
#' @param x = x vector
#' @param xi = location parameter
#' @param omega = scale parameter
#' @param alpha = skewness parameter
#' @param plot = logical, if \code{TRUE} produce la rappresentazione grafica delle densità 
#' e dell'area di sovrapposizione
#' @param return.all = logical, if \code{TRUE} restituisce il data set completo delle 
#' densità delle due distribuzioni
min_dskew_normal <- function( x = seq( -5, 5, by = .01 ), xi = 0, omega = 1, alpha = 0, 
                              plot = FALSE, return.all = FALSE ) {
  
  require( sn )
  require( ggplot2 )
  y1 <- dsn( x, xi = 0, omega = 1, alpha = 0 )
  y2 <- dsn( x, alpha = alpha, xi = xi, omega = omega )
  dy <- ifelse( y1 < y2, y1, y2 )
  gData <- data.frame( x, y1, y2, dy )  
  
  OVplot <- ggplot(gData, aes(x,dy) ) + theme_bw() +
    geom_ribbon( aes(x, ymin = 0, ymax = dy), fill = "orange") + 
    geom_line() + geom_line( aes(x,y1), color="red" ) + geom_line( aes(x,y2), color="blue" )
  
  if (plot) {  
    print( OVplot  )
  }
  
  if (return.all) {
    return( list( gData = gData, OVplot = OVplot ) )
  } else {
    return( dy )  
  }
}
# END

# +++++++++++++++++++++++++++++++++++++++++++++
#' @name maxdskew_normal_brms
#' @description Calcola il massimo tra due densità Skew-Normal
#' @note La prima densità ha i parametri fissati \code{xi = 0}, 
#' \code{omega = 1}, \code{alpha = 0}, per cui di fatto è normale standard
#' @param x = x vector
#' @param xi = location parameter
#' @param omega = scale parameter
#' @param alpha = skewness parameter
#' @param plot = logical, if \code{TRUE} produce la rappresentazione grafica delle densità 
#' e dell'area di sovrapposizione
#' @param return.all = logical, if \code{TRUE} restituisce il data set completo delle 
#' densità delle due distribuzioni
max_dskew_normal_bmrs <- function( x = seq( -5, 5, by = .01 ), xi = 0, omega = 1, alpha = 0, 
                              plot = FALSE, return.all = FALSE ) {
  
  require( brms )
  require( ggplot2 )
  y1 <- dskew_normal( x, xi = 0, omega = 1, alpha = 0 )
  y2 <- dskew_normal( x, alpha = alpha, xi = xi, omega = omega )
  dy <- ifelse( y1 > y2, y1, y2 )
  gData <- data.frame( x, y1, y2, dy )  
  
  OVplot <- ggplot(gData, aes(x,dy) ) + theme_bw() +
    geom_ribbon( aes(x, ymin = 0, ymax = dy), fill = "orange") + 
    geom_line() + geom_line( aes(x,y1), color="red" ) + geom_line( aes(x,y2), color="blue" )
  
  if (plot) {  
    print( OVplot  )
  }
  
  if (return.all) {
    return( list( gData = gData, OVplot = OVplot ) )
  } else {
    return( dy )  
  }
}
# END

# +++++++++++++++++++++++++++++++++++++++++++++
#' @name max_dskew_normal
#' @description Calcola il massimo tra due densità Skew-Normal
#' @note La prima densità ha i parametri fissati \code{xi = 0}, 
#' \code{omega = 1}, \code{alpha = 0}, per cui di fatto è normale standard
#' @param x = x vector
#' @param xi = location parameter
#' @param omega = scale parameter
#' @param alpha = skewness parameter
#' @param plot = logical, if \code{TRUE} produce la rappresentazione grafica delle densità 
#' e dell'area di sovrapposizione
#' @param return.all = logical, if \code{TRUE} restituisce il data set completo delle 
#' densità delle due distribuzioni
max_dskew_normal <- function( x = seq( -5, 5, by = .01 ), xi = 0, omega = 1, alpha = 0, 
                              plot = FALSE, return.all = FALSE ) {
  
  require( sn )
  require( ggplot2 )
  y1 <- dsn( x, xi = 0, omega = 1, alpha = 0 )
  y2 <- dsn( x, alpha = alpha, xi = xi, omega = omega )
  dy <- ifelse( y1 > y2, y1, y2 )
  gData <- data.frame( x, y1, y2, dy )  
  
  OVplot <- ggplot(gData, aes(x,dy) ) + theme_bw() +
    geom_ribbon( aes(x, ymin = 0, ymax = dy), fill = "orange") + 
    geom_line() + geom_line( aes(x,y1), color="red" ) + geom_line( aes(x,y2), color="blue" )
  
  if (plot) {  
    print( OVplot  )
  }
  
  if (return.all) {
    return( list( gData = gData, OVplot = OVplot ) )
  } else {
    return( dy )  
  }
}
# END

# +++++++++++++++++++++++++++++++++++++++++++++
#' @name true_cohen
true_cohen <- function( mu1, mu2, sigma1, sigma2 ) {
  pool.s <- sqrt( (sigma1^2 + sigma2^2) / 2 )
  mu1-mu2 / pool.s 
}
# END

# +++++++++++++++++++++++++++++++++++++++++++++
#' @name pars_skew_normal
pars_skew_normal <- function( xi = 0, omega = 1, alpha = 0) {
  delta <- alpha / sqrt( 1 + alpha^2 )
  mx <- xi + omega*delta*sqrt(2/pi)
  vx <- omega^2 * ( 1 - 2*delta^2 / pi )
  return( list( mean = mx, variance = vx ) )
}
# END

# +++++++++++++++++++++++++++
#' @name inv_pars_skew_normal
inv_pars_skew_normal <- function( mu=0, sigma=1, xi=NULL, omega=NULL, alpha=0 ) {
  
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
# END
#' @description Restituisce i parametri della SkewNormal, xi e omega a partire da media, dev.st. e alpha 

# +++++++++++++++++++++++++++++++++++++++++++++
#' @title core_sim
#' @description Simulation core. Simula i dati 
#' su due gruppi e calcola gli indici. Il primo gruppo
#' è campionato da una normale standard il secondo da una 
#' SkewNormal(delta,alpha,omega) 
#' @param n = sample size
#' @param delta = media della SkewNormal
#' @param alpha = skewness della SkewNormal
#' @param omega = variabilità della SkewNormal
#' @param B = number of replicates
core_sim_brms <- function(n, delta, alpha, omega, B) {
  
  true_overlap <- integrate( min_dskew_normal, -Inf, Inf, xi = delta, alpha = alpha, omega = omega )$value
  
  SN_par <- pars_skew_normal(xi=delta, omega=omega, alpha=alpha)
  true_d <- true_cohen( mu1 = 0, mu2 = SN_par$mean, sigma1 = 1, sigma2 = sqrt(SN_par$variance) )
  
  OUT <- t(sapply(1:B, function(b){
    y1 <- rskew_normal(n,xi=0,alpha=0, omega=1)
    y2 <- rskew_normal(n,xi=delta,alpha=alpha, omega=omega)
    
    COHEN <- cohens_d(y1,y2)
    HEDGES <- hedges_g(y1,y2)
    CLES <- d_to_overlap(COHEN$Cohens_d)
    U3 <- d_to_u3(COHEN$Cohens_d)
    PSUP <- d_to_p_superiority(COHEN$Cohens_d)
    
    ETA1 <- overlap(list(y1,y2))$OV
    ETA2 <- overlap(list(y1,y2), type="2")$OV
    
    output <- c( with( COHEN, c(Cohens_d,CI_low,CI_high) ),
                 with( HEDGES, c(Hedges_g,CI_low,CI_high) ),
                 CLES, U3, PSUP, 
                 ETA1, ETA2, n, delta, alpha, omega, true_overlap, true_d )
    
    names(output) <- c("Cohens_d","d_CI_low","d_CI_high","Hedges_g","g_CI_low","g_CI_high","CLES","U3","Prob_sup","eta1","eta2", "n", "delta","alpha","omega", "true_overlap", "true_d" )
    output
  }))
  
  return(OUT)
  
}
# END

# +++++++++++++++++++++++++++++++++++++++++++++
#' @title core_sim
#' @description Simulation core. Simula i dati 
#' su due gruppi e calcola gli indici. Il primo gruppo
#' è campionato da una normale standard il secondo da una 
#' SkewNormal(delta,alpha,omega) 
#' @param n = sample size
#' @param delta = media della SkewNormal
#' @param alpha = skewness della SkewNormal
#' @param omega = variabilità della SkewNormal
#' @param B = number of replicates
core_sim <- function(n, delta, alpha, omega, B) {
  
  true_overlap <- integrate( min_dskew_normal, -Inf, Inf, xi = delta, alpha = alpha, omega = omega )$value
  
  SN_par <- pars_skew_normal(xi=delta, omega=omega, alpha=alpha)
  true_d <- true_cohen( mu1 = 0, mu2 = SN_par$mean, sigma1 = 1, sigma2 = sqrt(SN_par$variance) )
  
  OUT <- t(sapply(1:B, function(b){
    y1 <- rsn(n,xi=0,alpha=0, omega=1)
    y2 <- rsn(n,xi=delta,alpha=alpha, omega=omega)
    
    COHEN <- cohens_d(y1,y2)
    HEDGES <- hedges_g(y1,y2)
    CLES <- d_to_overlap(COHEN$Cohens_d)
    U3 <- d_to_u3(COHEN$Cohens_d)
    PSUP <- d_to_p_superiority(COHEN$Cohens_d)
    
    ETA1 <- overlap(list(y1,y2))$OV
    ETA2 <- overlap(list(y1,y2), type="2")$OV
    
    output <- c( with( COHEN, c(Cohens_d,CI_low,CI_high) ),
                 with( HEDGES, c(Hedges_g,CI_low,CI_high) ),
                 CLES, U3, PSUP, 
                 ETA1, ETA2, n, delta, alpha, omega, true_overlap, true_d )
    
    names(output) <- c("Cohens_d","d_CI_low","d_CI_high","Hedges_g","g_CI_low","g_CI_high","CLES","U3","Prob_sup","eta1","eta2", "n", "delta","alpha","omega", "true_overlap", "true_d" )
    output
  }))
  
  return(OUT)
  
}
# END

