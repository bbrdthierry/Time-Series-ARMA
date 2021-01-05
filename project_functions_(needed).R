

##########################################################
# Graphique d'afc avec seuil pour un test global
##########################################################


acfG <- function(x, lag.max=NULL) {
  n <- length(x)
  if (is.null(lag.max)) 
    lag.max <- floor(10*log10(n))
  lag.max <- as.integer(min(lag.max, n-1L))
  rho <- acf(x, lag.max=lag.max, plot=F)$acf
  if ( all.equal(rho[1],1)==TRUE )
    rho <- rho[-1]
  nh <- length(rho) 
  out <- sum( abs(rho) > 1.96/sqrt(n) ) 
  alpha.global <- 0.05
  alpha <- 1 - (1-alpha.global)^(1/nh)
  seuil.global <- qnorm(1-alpha/2, sd=1) / sqrt(n)
  acf(x, lag.max=lag.max)
  abline(h=seuil.global, lty=3, col="red")
  abline(h=-seuil.global, lty=3, col="red")
  #
  testG <- binom.test(out, n=nh, p=0.05)
  return(testG)
}


#########################################################
# Test d'indï¿½pendance des points tournants
#########################################################



PtTourn.test <- function (vec)
{
  n <- length(vec)
  nT <- 0
  for (i in 2:(n-1)){
    if (((vec[i] > vec[i-1]) & (vec[i] > vec[i+1])) | ( (vec[i] < vec[i-1])
                                                        & (vec[i] < vec[i+1]))) 
    {nT <- nT+1}
    
  }
  Tobs <- ( nT - 2*(n-2)/3 ) / sqrt( (16*n-29)/90 )
  p <- 2 * (1- pnorm(abs(Tobs)))
  res <- c(n,nT,Tobs,p)
  names(res) <- c("n","nT","stat"," p-valeur")
  return(res)
}

periodog <- function(x, graph = T) {
  freqs <- (2 * pi * (0:floor(length(x)/2)))/length(x)
  periodo <- (Mod(fft(x - mean(x)))^2/(2 * pi * length(x)))[1:(floor(length(x)/2) +
                                                                 1)]
  if (graph == T) {
    plot(freqs, periodo, type = "l", main = "Periodogram")
  }
  else {
    return(data.frame(periodo = periodo, freqs = freqs))
  }
}

Saison.test <- function(vec, d) {
  n <- length(vec)
  ns <- n%/%d
  n <- ns * d
  vec <- vec[1:n]
  rangs.vec <- matrix(rep(0, n), nrow = d, ncol = ns)
  for (j in 1:ns) {
    saisonj <- vec[(d * (j - 1) + 1):(d * j)]
    rangs.vec[order(saisonj), j] <- 1:d
  }
  rangs <- apply(rangs.vec, 1, sum)
  stat <- 12 * sum((rangs - ns * (d + 1)/2)^2)/(ns * d * (d + 1))
  pval <- 1 - pchisq(stat, df = d - 1)
  res <- c(d, d - 1, stat, pval)
  names(res) <- c("Periode", "d.f.", "Tobs", "p-valeur")
  return(res)
}


PtMont.test <- function (vec)
{
  n <- length(vec)
  nM <- 0
  for (i in 1:(n-1)){
    if (vec[i] < vec[i+1]) {
      nM <- nM+1         
    }}
  Tobs <- ( nM - (n-1)/2 ) / sqrt( (n+1)/12 )
  p <- 2 * (1- pnorm(abs(Tobs)))
  res <- c(n,nM,Tobs,p)
  names(res) <- c("n","nM","stat"," p-valeur")
  return(res)
}


PtDisc.test <- function (vec)
{
  n <- length(vec)
  res <- cor.test(1:n, vec, method="kendall")
  nD=floor((1-res$estimate)/4*n*(n-1))
  Tobs <- ( 1 - 4*nD/(n*(n-1)) ) / sqrt( 2*(2*n+5)/(9*n*(n-1)) )
  p <- 2 * (1- pnorm(abs(Tobs)))
  res <- c(n,nD,Tobs,p)
  names(res) <- c("n","nD","stat"," p-valeur")
  return(res)
}

acfMA <- function(x) {
  n <- length(x)
  rho <- acf(x, plot = F)$acf[-1]
  nlags <- length(acf(x, plot = F)$lag)
  wh <- c(1, 1 + 2 * cumsum(rho^2))
  seuilsSup <- 1.96 * sqrt(wh/n)
  seuilsSup <- as.ts(seuilsSup)
  seuilsSup <- ts(seuilsSup, start = 1/frequency(x),
                  end = nlags/frequency(x), frequency = frequency(x))
  acf(x)
  lines(seuilsSup, lty = 2, col = "red",
        lwd = 2)
  lines(-seuilsSup, lty = 2, col = "red",
        lwd = 2)
}



##########################################################
# Test d'indépendance de Portmanteau
##########################################################



Portmant.test <- function(vec)
{
  n <- length(vec)
  autocor <- acf(vec, plot = FALSE, na.action = na.pass)
  H <- floor(n/4)
  H <- min( length(autocor$acf) -1 , H)   
  stat <- n * sum(autocor$acf[2:(H+1)] ^2 )
  pval <- 1 - pchisq(stat, df=H)
  res <- c(H,stat,pval)
  names(res) <- c("d.f.","stat"," p-value")
  return(res)
}


PortmantLB.test <- function(vec)
{
  n <- length(vec)
  autocor <- acf(vec, plot = FALSE, na.action = na.pass)
  H <- floor(n/4)
  H <- min( length(autocor$acf) -1 , H) 
  h <- 1:H  
  stat <- n * (n+2) * sum(autocor$acf[2:(H+1)] ^2 / (n-h) )
  pval <- 1 - pchisq(stat, df=H)
  res <- c(H,stat,pval)
  names(res) <- c("d.f.","stat"," p-value")
  return(res)
}
