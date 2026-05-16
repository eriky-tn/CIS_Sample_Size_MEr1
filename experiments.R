################################################################################
# Title: A Bayesian Approach to Determine Sample Size in Erlang Single-Server Queues
#
# Reference:
# Singh, L., Gomes, E. S., & Cruz, F. R. B. (2025).
# A Bayesian approach to determine sample size in Erlang single-server queues.
# Communications in Statistics - Theory and Methods, 54(12), 3692-3710.
# https://doi.org/10.1080/03610926.2024.2401602
#
# Authors:
# Eriky S. Gomes
# Frederico R. B. Cruz
#
# Affiliation:
# Department of Statistics
# Universidade Federal de Minas Gerais (UFMG)
#
# Contact:
# eriky-tn@ufmg.br
# fcruz@est.ufmg.br
#
# Copyright (c) 2025 Gomes & Cruz
# Version: v2026
#
# Description:
# R script developed for computational experiments of the article
################################################################################
# Auxiliar functions
################################################################################

# optimizer by golden-section algorithm
MaxFuncGS <- function(func, ...){
  TOL = 1e-3                            
  GOLDEN = (sqrt(5) - 1) / 2  # inverse of golden number
  iterMax = ceiling(log(TOL) / log(GOLDEN))
  a <- TOL
  b <- 1 - TOL
  x1 <- a + (1 - GOLDEN) * (b - a)
  x2 <- a + GOLDEN * (b - a)
  fx1 <- func(x1, ...)
  fx2 <- func(x2, ...)
  for(i in 1:iterMax){
    if(fx1 < fx2){
      a <- x1
      x1 <- x2
      fx1 <- fx2
      x2 <- a + GOLDEN * (b - a)
      fx2 <- func(x2, ...)
    }
    else{
      b <- x2
      x2 <- x1
      fx2 <- fx1
      x1 <- a + (1 - GOLDEN) * (b - a)
      fx1<-func(x1,...)
    }
  }
  xMax <- (a + b) / 2
  fMax <- func(xMax, ...)
  return(list(xMax = xMax, fMax = fMax))
}

# likelihood function
Lik <- function(p, x, r){
  y <- sum(x)
  n <- length(x)
  const <- r^(-y)
  for(i in 1:n){
    const = const * choose(r + x[i] - 1, x[i])
  }
  return(const * p^y * (1 + p / r)^(-y - n * r))
}

# log-likelihood function
LogLik <- function(p, x, r){
  y <- sum(x)
  n <- length(x)
  const <- (-y) * log(r)
  for(i in 1:n){
    const = const + log(choose(r + x[i] - 1, x[i]))
  }
  #cat('const = ', const,'\n')
  return(const + y * log(r) + (-y - n * r) * log((1 + p / r)))
}

# maximum likelihood estimator
MLE <-function(x, r){
  tol <- 1e-5
  y <- sum(x)
  n <- length(x)
  if(y / n >= 1){
    return(1 - tol)
  }
  else{
    return(y / n)
  }
}

GaussHyp <- function(a, b, c, z) {
  if (any(c(b, c) <= 0) || c <= b || z < -1) stop('invalid parameters')
  
  AuxF <- function(u, a, b, c, z) {
    auxOut <- u^(b - 1) * (1 - u)^(c - b - 1) * (1 - z * u)^(-a)
    return(auxOut)
  }
  return(integrate(AuxF, 0, 1, a , b, c, z)[[1]] / 
           beta(b, c - b))
}


# mean of distribution function
MeanF <- function(dtb, ...){
  AuxF<-function(p, ...){
    auxOut <- p * dtb(p, ...)
    return(auxOut)
  }
  return(integrate(AuxF, 0, 1, ...)[[1]])
}


Ixi <- function(a, b, xi){
  AuxF <- function(u){
    return(u^(a - 1) * (1 - u)^(b - 1))
  }
  return(integrate(AuxF, 0, xi)[[1]] / beta(a, b))
}

# bissection root-finding
BissectionRF <- function(f, lInf, lSup, ...){
  TOL = 1e-3
  k <- 0
  while(abs(lSup - lInf) > TOL){
    f1 <- f(lInf, ...)
    f2 <- f(lSup, ...)
    x <- (lInf + lSup) / 2
    fx <- f(x, ...)
    ifelse(f1 * fx < 0, lSup <- x, lInf <- x)
    k <- k + 1
  }
  root <- (lInf + lSup) / 2
  #cat(k,'iterations')
  return(root)
}

# random number generator by acceptance-rejection method
RandNumAR <- function(n, dtb, ...){
  TOL = 1e-2
  randNum <- numeric(n)
  fmax <- MaxFuncGS(dtb, ...)$fMax
  i <- 0
  while(i <= n){
    u <- runif(1, TOL, 1 - TOL)
    y <- runif(1, TOL, 1 - TOL)
    if(u < dtb(y, ...) / fmax){
      randNum[i] <- y
      i <- i + 1
    }
  }
  return(randNum)
}

################################################################################
# Beta functions
################################################################################

# beta prior
BetaPrior <- function(p, r, a, b){ # r is dummy
  result <- p^(a - 1) * (1 - p)^(b - 1) / beta(a, b)
  if(any(!is.finite(result))) stop('non-finite result')
  return(result)
}

# beta posterior
BetaPosterior <- function(p, x, r, a, b){
  y <- sum(x)
  n <- length(x)
  result <- p^(y + a - 1) * (1 - p)^(b - 1) * (1 + p / r)^(-y - n * r) /
           beta(y + a, b) / GaussHyp(y + n * r, y + a, y + a + b, -1 / r)
  if(any(!is.finite(result))) stop('non-finite result')
  return(result)
}

# beta self estimator
BetaSelf <- function(x, r, a, b){
  y <- sum(x)
  n <- length(x)
  result <- (y + a) / (y + a + b) *
    GaussHyp(y + n * r, y + a + 1, y + a + 1 + b, -1 / r) / 
    GaussHyp(y + n * r, y + a, y + a + b, -1 / r)
  if(any(!is.finite(result))) stop('non-finite result')
  return(result)
}

################################################################################
# Gauss hypergeometric functions
################################################################################

# gauss hypergeometric prior
GaussHypPrior <- function(p, r, a, b, c){
  result <- p^(b - 1) * (1 - p)^(c - b - 1) * (1 + p / r)^(-a) /
    beta(b, c - b) / GaussHyp(a, b, c, -1 / r)
  if(any(!is.finite(result))) stop('non-finite result')
  return(result)
}

# gauss hypergeometric posterior
GaussHypPosterior <- function(p, x, r, a, b, c){
  y <- sum(x)
  n <- length(x)
  result <- p^(b - 1 + y) * (1 - p)^(c - b - 1) * (1 + p / r)^(-a - y - n * r) /
    beta(b + y, c - b) / GaussHyp(a + y + n * r, b + y, c + y, -1 / r)
  if(any(!is.finite(result))) stop('non-finite result')
  return(result)
}

# gauss hypergeometric self estimator
GaussHypSelf <- function(x, r, a, b, c){
  y <- sum(x)
  n <- length(x)
  result <- (b + y) / (c + y) * 
    GaussHyp(a + y + n * r, b + y + 1, c + y + 1, -1 / r) /
    GaussHyp(a + y + n * r, b + y, c + y, -1 / r)
  if(any(!is.finite(result))) stop('non-finite result')
  return(result)
}

################################################################################
# Jeffreys functions
################################################################################

# jeffreys prior
JefPrior <- function(p, r){
  result <- r ^ (1 / 2) * p ^ (-1 / 2) *(r + p) ^(-1 / 2)
  if(any(!is.finite(result))) stop('non-finite result')
  return(result)
}

JefPosterior <- function(p, x, r){
  y <- sum(x)
  n <- length(x)
  result <- r^(-1) * (p / r)^(y - 1 / 2) * (1 + p / r)^(-y - n * r - 1 / 2) /
    beta(y + 1 / 2, n * r) /
           Ixi(y + 1 / 2, n * r, 1 / (1 + r))
  if(any(!is.finite(result))) stop('non-finite result')
  return(result)
}

# jeffreys self estimator
JefSelf <- function(x, r){
  y <- sum(x)
  n <- length(x)
  result <- r * (y + 1 / 2) / (n * r - 1) *
           Ixi(y + 3 / 2, n * r - 1, 1 / (1 + r)) /
           Ixi(y + 1 / 2, n * r, 1 / (1 + r))
  if(any(!is.finite(result))) stop('non-finite result')
  return(result)
           
}


JefSelf2 <- function(x, r) {
  y <- sum(x)
  n <- length(x)
  result <-
    2 * r ^ (-y  - 1 / 2) / (2 * y + 3) * 
    GaussHyp(y + 3 / 2, n * r + y + 1 / 2, y + 5 / 2,-1 / r) /
    beta(y + 1 / 2, n * r) /
    Ixi(y + 1 / 2, n * r, 1 / (1 + r))
  if (any(!is.finite(result)))
    stop('non-finite result')
  return(result)
}

################################################################################
# Credible intervals
################################################################################

# equal-tailed credible interval
ETI <- function(postF, ...){
  DELTA = 1e-3
  ALPHA = 0.05
  p = seq(DELTA, 1 - DELTA, by = DELTA)
  len = length(p)
  f <- postF(p, ...)
  dArea <- 0
  cover <- 0
  
  areaLeft <- 0
  i <- 1
  while(areaLeft < ALPHA / 2 & i < len){
    dArea <- DELTA * (f[i] + f[i + 1]) / 2
    areaLeft <- areaLeft + dArea
    i <- i + 1
  }
  lInf <- p[i - 1]
  areaLeft <- areaLeft - dArea # area correction
  
  areaRight <- 0
  j <- len
  while(areaRight < ALPHA / 2 & j > i){
    dArea <- DELTA * (f[j] + f[j - 1]) / 2
    areaRight <- areaRight + dArea
    j <- j - 1
  }
  lSup <- p[j + 1]
  areaRight <- areaRight - dArea # area correction
  
  cover <- 1 - areaLeft - areaRight
  return(list(lInf = lInf, lSup = lSup, cover = cover))
}

# symmetrical credible interval (from pontual estimate)
SCI <- function(DELTA, pEstF, postF, ...){
  pEst <- pEstF(...)
  lInf <- max(0, pEst - DELTA)
  lSup <- min(1, pEst + DELTA)
  cover <- integrate(postF, lInf, lSup, subdivisions = 1000, ...)[[1]]
  return(list(lInf = lInf, lSup = lSup, cover = cover))
}

################################################################################
# Sample size criteria
################################################################################

# average coverage criterion
ACC <- function(DELTA, r, priori, post, pEstF, ...){
  ALPHA = 0.05
  k = 1000 # samples drawn for each n
  nMin = 2
  nMax = 300
  p <- RandNumAR(k, priori, r, ...)

  while(nMax - nMin > 1){
    An <- 0
    n <- round((nMin + nMax) / 2, 0)
    i <- 1
    while(i <= k){
      x <- rnbinom(n, r, 1 / (1 + p[i] / r))
      interv <- SCI(DELTA, pEstF, post, x, r, ...)
      cover <- interv$cover
      An <- An + cover
      i <- i + 1
    }
    An <- An / k
    ifelse(An >= 1 - ALPHA , nMax <- n, nMin <- n)
    #cat("n=", n, "\n")
  }
  n <- nMax
  if(n >= 295){
    #cat('n = ', Inf, '\n')
    return(Inf)
  }
  return(n)
}

# average length criterion
ALC <- function(WIDTH, r, priori, post, ...){
  k = 1000 # samples drawn for each n
  nMin = 2
  nMax = 300
  p <- RandNumAR(k, priori, r, ...)
  
  while(nMax - nMin > 1){
    Wn <- 0
    n <- round((nMax + nMin) / 2, 0)
    i <- 1
    while(i <= k){
      x <- rnbinom(n, r, 1 / (1 + p[i] / r))
      interv <- ETI(post, x, r, ...)
      Wn <- Wn + (interv$lSup - interv$lInf)
      i <- i + 1
    }
    Wn <- Wn / k
    ifelse(Wn > WIDTH, nMin <- n, nMax <- n)
    #cat("n=", n, "\n")
  }
  n <- nMax
  if(n >= 295){
    cat('n = ', Inf, '\n')
    return(Inf)
  }
  return(n)
}


################################################################################
# Table functions
################################################################################

# table function for average coverage criterion
TabACC <- function(DELTA, r, prior, posterior, pEstF, ...){
  tab <- numeric(length(DELTA))
  for(i in 1:length(DELTA)){
    set.seed(2023)
    tab[i] <- ACC(DELTA[i], r, prior, posterior, pEstF, ...)
    cat(round(i / length(DELTA) * 100, 0), "%\n")
  }
  return(tab)
}

# table function for average length criterion
TabALC <- function(WIDTH, r, prior, posterior, ...){
  tab <- numeric(length(WIDTH))
  for(i in 1:length(WIDTH)){
    set.seed(2023)
    tab[i] <- ALC(WIDTH[i], r, prior, posterior, ...)
    cat(round(i / length(WIDTH) * 100, 0), "%\n")
  }
  return(tab)
}

################################################################################
# Tables for abacuses
################################################################################

# equal-tailed credible interval table
TabETI <- function(r, priori, post, size, ...){
  set.seed(2023)
  rep <- 1000
  p <- RandNumAR(rep, priori, r, ...)
  width <- numeric(rep)
  tab <- matrix(nrow = length(size), ncol = 1, dimnames = list(size, 'width'))
  
  for(i in 1:length(size)){
    for(j in 1:rep){
      x <- rnbinom(size[i], r, 1 / (1 + p[j] / r))
      interv <- ETI(post, x, r, ...)
      width[j] <- interv$lSup - interv$lInf
    }
    tab[i, 1] <- mean(width)
    cat(round(i / length(size) * 100, 0), "%\n")
  }
  return(tab)
}

# symmetrical credible interval table
TabSCI <- function(delta, r, prior, post, pEstF, size, ...){
  set.seed(2023)
  rep <- 1000
  p <- RandNumAR(rep, prior, r, ...)
  cover <- numeric(rep)
  tab <- matrix(nrow = length(size), ncol = length(delta),
                dimnames = list(size, delta))
  
  for(i in 1:length(size)){
    for(j in 1:length(delta)){
      for(k in 1:rep){
        x <- rnbinom(size[i], r, 1 / (1 + p[k] / r))
        cover[k] <- SCI(delta[j], pEstF, post, x, r, ...)$cover
      }
      tab[i, j] <- mean(cover)
    }
    cat(round(i / length(size) * 100, 0), "%\n")
  }
  return(tab)
}

################################################################################
# Simulation study: tables for abacuses
################################################################################

# r = 2 ------------------------------------------------------------------------
size <- seq(10, 400, 10)
DELTA <- c(0.08, 0.10, 0.12, 0.15)
r <- 2

# beta prior
a <- c(1, 1, 2, 2)
b <- c(1, 2, 1, 2)

tabB1_ETI <- TabETI(r, BetaPrior, BetaPosterior, size, a[1], b[1])
tabB2_ETI <- TabETI(r, BetaPrior, BetaPosterior, size, a[2], b[2])
tabB3_ETI <- TabETI(r, BetaPrior, BetaPosterior, size, a[3], b[3])
tabB4_ETI <- TabETI(r, BetaPrior, BetaPosterior, size, a[4], b[4])

write.table(cbind(tabB1_ETI,tabB2_ETI,tabB3_ETI,tabB4_ETI),'tabBeta_ETI_2.txt',
          sep = '\t', quote = FALSE)
rm(tabB1_ETI, tabB2_ETI, tabB3_ETI, tabB4_ETI)

tabB1_SCI <- TabSCI(DELTA, r, BetaPrior, BetaPosterior, BetaSelf, size, a[1], b[1])
tabB2_SCI <- TabSCI(DELTA, r, BetaPrior, BetaPosterior, BetaSelf, size, a[2], b[2])
tabB3_SCI <- TabSCI(DELTA, r, BetaPrior, BetaPosterior, BetaSelf, size, a[3], b[3])
tabB4_SCI <- TabSCI(DELTA, r, BetaPrior, BetaPosterior, BetaSelf, size, a[4], b[4])

write.table(cbind(tabB1_SCI,tabB2_SCI,tabB3_SCI,tabB4_SCI),'tabBeta_SCI_2.txt',
            sep = '\t', quote = FALSE)
rm(a, b, tabB1_SCI, tabB2_SCI, tabB3_SCI, tabB4_SCI)


# gauss hypergeometric prior
a <- c(1, 1, 1, 1)
b <- c(1, 2, 3, 4)
c <- c(3, 6, 4, 6)

tabG1_ETI <- TabETI(r, GaussHypPrior, GaussHypPosterior, size, a[1], b[1], c[1])
tabG2_ETI <- TabETI(r, GaussHypPrior, GaussHypPosterior, size, a[2], b[2], c[2])
tabG3_ETI <- TabETI(r, GaussHypPrior, GaussHypPosterior, size, a[3], b[3], c[3])
tabG4_ETI <- TabETI(r, GaussHypPrior, GaussHypPosterior, size, a[4], b[4], c[4])

write.table(cbind(tabG1_ETI,tabG2_ETI,tabG3_ETI,tabG4_ETI),'tabGau_ETI_2.txt',
            sep = '\t', quote = FALSE)
rm(tabG1_ETI, tabG2_ETI, tabG3_ETI, tabG4_ETI)

tabG1_SCI <- TabSCI(DELTA, r, GaussHypPrior, GaussHypPosterior, GaussHypSelf, size, a[1], b[1], c[1])
tabG2_SCI <- TabSCI(DELTA, r, GaussHypPrior, GaussHypPosterior, GaussHypSelf, size, a[2], b[2], c[2])
tabG3_SCI <- TabSCI(DELTA, r, GaussHypPrior, GaussHypPosterior, GaussHypSelf, size, a[3], b[3], c[3])
tabG4_SCI <- TabSCI(DELTA, r, GaussHypPrior, GaussHypPosterior, GaussHypSelf, size, a[4], b[4], c[4])

write.table(cbind(tabG1_SCI,tabG2_SCI,tabG3_SCI,tabG4_SCI),'tabGau_SCI_2.txt',
            sep = '\t', quote = FALSE)
rm(a, b, c, tabG1_SCI, tabG2_SCI, tabG3_SCI, tabG4_SCI)

# jeffreys prior
size <- seq(10, 250, 10)
tabJ_ETI <- TabETI(r, JefPrior, JefPosterior, size)
write.table(tabJ_ETI,'tabJef_ETI_2.txt', sep = '\t', quote = FALSE)
rm(tabJ_ETI)

size <- seq(10, 300, 10)
tabJ_SCI <- TabSCI(DELTA, r, JefPrior, JefPosterior, JefSelf, size)
write.table(tabJ_SCI,'tabJef_SCI_2.txt', sep = '\t', quote = FALSE)
rm(tabJ_SCI)

# r = 4 ------------------------------------------------------------------------
size <- seq(10,400,10)
DELTA <- c(0.08, 0.10, 0.12, 0.15)
r <- 4

# beta prior
a <- c(1, 1, 2, 2)
b <- c(1, 2, 1, 2)

tabB1_ETI <- TabETI(r, BetaPrior, BetaPosterior, size, a[1], b[1])
tabB2_ETI <- TabETI(r, BetaPrior, BetaPosterior, size, a[2], b[2])
tabB3_ETI <- TabETI(r, BetaPrior, BetaPosterior, size, a[3], b[3])
tabB4_ETI <- TabETI(r, BetaPrior, BetaPosterior, size, a[4], b[4])

write.table(cbind(tabB1_ETI,tabB2_ETI,tabB3_ETI,tabB4_ETI),'tabBeta_ETI_4.txt',
            sep = '\t', quote = FALSE)
rm(tabB1_ETI, tabB2_ETI, tabB3_ETI, tabB4_ETI)

tabB1_SCI <- TabSCI(DELTA, r, BetaPrior, BetaPosterior, BetaSelf, size, a[1], b[1])
tabB2_SCI <- TabSCI(DELTA, r, BetaPrior, BetaPosterior, BetaSelf, size, a[2], b[2])
tabB3_SCI <- TabSCI(DELTA, r, BetaPrior, BetaPosterior, BetaSelf, size, a[3], b[3])
tabB4_SCI <- TabSCI(DELTA, r, BetaPrior, BetaPosterior, BetaSelf, size, a[4], b[4])

write.table(cbind(tabB1_SCI,tabB2_SCI,tabB3_SCI,tabB4_SCI),'tabBeta_SCI_4.txt',
            sep = '\t', quote = FALSE)
rm(a, b, tabB1_SCI, tabB2_SCI, tabB3_SCI, tabB4_SCI)


# gauss hypergeometric prior
a <- c(1, 1, 1, 1)
b <- c(1, 2, 3, 4)
c <- c(3, 6, 4, 6)

tabG1_ETI <- TabETI(r, GaussHypPrior, GaussHypPosterior, size, a[1], b[1], c[1])
tabG2_ETI <- TabETI(r, GaussHypPrior, GaussHypPosterior, size, a[2], b[2], c[2])
tabG3_ETI <- TabETI(r, GaussHypPrior, GaussHypPosterior, size, a[3], b[3], c[3])
tabG4_ETI <- TabETI(r, GaussHypPrior, GaussHypPosterior, size, a[4], b[4], c[4])

write.table(cbind(tabG1_ETI,tabG2_ETI,tabG3_ETI,tabG4_ETI),'tabGau_ETI_4.txt',
            sep = '\t', quote = FALSE)
rm(tabG1_ETI, tabG2_ETI, tabG3_ETI, tabG4_ETI)

tabG1_SCI <- TabSCI(DELTA, r, GaussHypPrior, GaussHypPosterior, GaussHypSelf, size, a[1], b[1], c[1])
tabG2_SCI <- TabSCI(DELTA, r, GaussHypPrior, GaussHypPosterior, GaussHypSelf, size, a[2], b[2], c[2])
tabG3_SCI <- TabSCI(DELTA, r, GaussHypPrior, GaussHypPosterior, GaussHypSelf, size, a[3], b[3], c[3])
tabG4_SCI <- TabSCI(DELTA, r, GaussHypPrior, GaussHypPosterior, GaussHypSelf, size, a[4], b[4], c[4])

write.table(cbind(tabG1_SCI,tabG2_SCI,tabG3_SCI,tabG4_SCI),'tabGau_SCI_4.txt',
            sep = '\t', quote = FALSE)
rm(a, b, c, tabG1_SCI, tabG2_SCI, tabG3_SCI, tabG4_SCI)

# jeffreys prior
size <- seq(10, 250, 10)
tabJ_ETI <- TabETI(r, JefPrior, JefPosterior, size)
write.table(tabJ_ETI,'tabJef_ETI_4.txt', sep = '\t', quote = FALSE)
rm(tabJ_ETI)

size <- seq(10, 250, 10)
tabJ_SCI <- TabSCI(DELTA, r, JefPrior, JefPosterior, JefSelf, size)
write.table(tabJ_SCI,'tabJef_SCI_4.txt', sep = '\t', quote = FALSE)
rm(tabJ_SCI)

################################################################################
# Simulation study: sample size tables
################################################################################

# r = 2 ------------------------------------------------------------------------
size <- seq(10, 400,10)
DELTA <- c(0.08, 0.10, 0.12, 0.15)
WIDTH <- c(0.16, 0.20, 0.24, 0.30)
r <- 2

# beta prior
matrixACC_Beta <- matrix(NA, nrow = 0, ncol = length(DELTA))
matrixALC_Beta <- matrix(NA, nrow = 0, ncol = length(WIDTH))
colnames(matrixACC_Beta) = DELTA
colnames(matrixALC_Beta) = WIDTH

a <- c(1, 1, 2, 2)
b <- c(1, 2, 1, 2)

for(i in 1:length(a)){
  matrixACC_Beta <- rbind(matrixACC_Beta,
    TabACC(DELTA, r, BetaPrior, BetaPosterior, BetaSelf, a[i], b[i]))
  matrixALC_Beta <- rbind(matrixALC_Beta,
    TabALC(WIDTH, r, BetaPrior, BetaPosterior, a[i], b[i]))
}
rownames(matrixACC_Beta) = paste0('Beta(', a, ',', b, ')')
rownames(matrixALC_Beta) = paste0('Beta(', a, ',', b, ')')

# gauss hypergeometric prior
matrixACC_Gauss <- matrix(NA, nrow = 0, ncol = length(DELTA))
matrixALC_Gauss <- matrix(NA, nrow = 0, ncol = length(WIDTH))
colnames(matrixACC_Gauss) = DELTA
colnames(matrixALC_Gauss) = WIDTH

a <- c(1, 1, 1, 1)
b <- c(1, 2, 3, 4)
c <- c(3, 6, 4, 6)

for(i in 1:length(a)){
  matrixACC_Gauss <- rbind(matrixACC_Gauss,
                     TabACC(DELTA, r, GaussHypPrior, GaussHypPosterior,
                            GaussHypSelf, a[i], b[i], c[i]))
  matrixALC_Gauss <- rbind(matrixALC_Gauss,
                     TabALC(WIDTH, r, GaussHypPrior, GaussHypPosterior,
                            a[i], b[i], c[i]))
}
rownames(matrixACC_Gauss) = paste0('Gauss(', a, ',', b, ',', c, ')')
rownames(matrixALC_Gauss) = paste0('Gauss(', a, ',', b, ',', c, ')')

# jeffreys prior
size <- seq(10,300,10)
matrixACC_Jeffreys <- matrix(NA, nrow = 0, ncol = length(DELTA))
matrixALC_Jeffreys <- matrix(NA, nrow = 0, ncol = length(WIDTH))
colnames(matrixACC_Jeffreys) = DELTA
colnames(matrixALC_Jeffreys) = WIDTH

matrixACC_Jeffreys <- rbind(matrixACC_Jeffreys,
                   TabACC(DELTA, r, JefPrior, JefPosterior, JefSelf))
matrixALC_Jeffreys <- rbind(matrixALC_Jeffreys,
                   TabALC(WIDTH, r, JefPrior, JefPosterior))

rownames(matrixACC_Jeffreys) = paste0('Jeffreys')
rownames(matrixALC_Jeffreys) = paste0('Jeffreys')

# writing results
write.table(rbind(matrixACC_Beta,matrixACC_Gauss,matrixACC_Jeffreys),
      file = 'TableACC_2.txt', sep = '\t', quote = FALSE, row.names = TRUE)
write.table(rbind(matrixALC_Beta,matrixALC_Gauss,matrixALC_Jeffreys),
            file = 'TableALC_2.txt', sep = '\t', quote = FALSE, row.names = TRUE)
rm(matrixACC_Beta, matrixACC_Gauss,matrixACC_Jeffreys)
rm(matrixALC_Beta, matrixALC_Gauss,matrixALC_Jeffreys)
rm(r, a, b, c, DELTA, WIDTH)

# r = 4 ------------------------------------------------------------------------
DELTA <- c(0.08, 0.10, 0.12, 0.15)
WIDTH <- c(0.16, 0.20, 0.24, 0.30)
r <- 4

# beta prior
matrixACC_Beta <- matrix(NA, nrow = 0, ncol = length(DELTA))
matrixALC_Beta <- matrix(NA, nrow = 0, ncol = length(WIDTH))
colnames(matrixACC_Beta) = DELTA
colnames(matrixALC_Beta) = WIDTH

a <- c(1, 1, 2, 2)
b <- c(1, 2, 1, 2)

for(i in 1:length(a)){
  matrixACC_Beta <- rbind(matrixACC_Beta,
                          TabACC(DELTA, r, BetaPrior, BetaPosterior, BetaSelf, a[i], b[i]))
  matrixALC_Beta <- rbind(matrixALC_Beta,
                          TabALC(WIDTH, r, BetaPrior, BetaPosterior, a[i], b[i]))
}
rownames(matrixACC_Beta) = paste0('Beta(', a, ',', b, ')')
rownames(matrixALC_Beta) = paste0('Beta(', a, ',', b, ')')

# gauss hypergeometric prior
matrixACC_Gauss <- matrix(NA, nrow = 0, ncol = length(DELTA))
matrixALC_Gauss <- matrix(NA, nrow = 0, ncol = length(WIDTH))
colnames(matrixACC_Gauss) = DELTA
colnames(matrixALC_Gauss) = WIDTH

a <- c(1, 1, 1, 1)
b <- c(1, 2, 3, 4)
c <- c(3, 6, 4, 6)

for(i in 1:length(a)){
  matrixACC_Gauss <- rbind(matrixACC_Gauss,
                           TabACC(DELTA, r, GaussHypPrior, GaussHypPosterior,
                                  GaussHypSelf, a[i], b[i], c[i]))
  matrixALC_Gauss <- rbind(matrixALC_Gauss,
                           TabALC(WIDTH, r, GaussHypPrior, GaussHypPosterior,
                                  a[i], b[i], c[i]))
}
rownames(matrixACC_Gauss) = paste0('Gauss(', a, ',', b, ',', c, ')')
rownames(matrixALC_Gauss) = paste0('Gauss(', a, ',', b, ',', c, ')')

# jeffreys prior 
matrixACC_Jeffreys <- matrix(NA, nrow = 0, ncol = length(DELTA))
matrixALC_Jeffreys <- matrix(NA, nrow = 0, ncol = length(WIDTH))
colnames(matrixACC_Jeffreys) = DELTA
colnames(matrixALC_Jeffreys) = WIDTH

# (must bound ACC and ALC on nMax = 300)
matrixACC_Jeffreys <- rbind(matrixACC_Jeffreys,
                            TabACC(DELTA, r, JefPrior, JefPosterior, JefSelf))
matrixALC_Jeffreys <- rbind(matrixALC_Jeffreys,
                            TabALC(WIDTH, r, JefPrior, JefPosterior))

rownames(matrixACC_Jeffreys) = paste0('Jeffreys')
rownames(matrixALC_Jeffreys) = paste0('Jeffreys')

# writing results
write.table(rbind(matrixACC_Beta,matrixACC_Gauss,matrixACC_Jeffreys),
            file = 'TableACC_4.txt', sep = '\t', quote = FALSE, row.names = TRUE)
write.table(rbind(matrixALC_Beta,matrixALC_Gauss,matrixALC_Jeffreys),
            file = 'TableALC_4.txt', sep = '\t', quote = FALSE, row.names = TRUE)
rm(matrixACC_Beta, matrixACC_Gauss,matrixACC_Jeffreys)
rm(matrixALC_Beta, matrixALC_Gauss,matrixALC_Jeffreys)
rm(r, a, b, c, DELTA, WIDTH)

#save.image('mer1-ss.rdata')
#load('mer1-ss.rdata')
################################################################################
# TESTS
################################################################################

# testing optimizer by golden-section algorithm
rep <- 1000
a <- runif(rep, 1, 5)
b <- runif(rep, 1, 5)
r <- 1
maxGS <- numeric(rep)
error <- numeric(rep)
for(i in 1: rep){
  error[i] <- MaxFuncGS(BetaPrior, r, a[i], b[i])$xMax - 
    (a[i] - 1) / (a[i] + b[i] - 2)
}
all(is.finite(error)) # must be TRUE
max(abs(error)) # maximum error
rm(rep, a, b, r, maxGS, error, i)

# testing likelihood function and mle
n = 50; r = 1; p = 0.5
x <- rnbinom(n, r, 1 / (1 + p / r))
MLE(x, r) - MaxFuncGS(Lik, x, r)$xMax # error
MLE(x,r)
MaxFuncGS(Lik, x, r)[1]
p <- seq(0, 1, 1e-3)
plot(p, Lik(p, x, r),'l')
abline(v = MLE(x,r))

n = 100; r = 2; p = 0.7
x <- rnbinom(n, r, 1 / (1 + p / r))
MLE(x, r) - MaxFuncGS(Lik, x, r)$xMax # error
p <- seq(0,1,1e-3)
plot(p, Lik(p, x, r),'l')
abline(v = MLE(x,r))

n = 200; r = 4; p = 0.2
x <- rnbinom(n, r, 1 / (1 + p / r))
MLE(x, r) - MaxFuncGS(Lik, x, r)$xMax # error
p <- seq(0,1,1e-3)
plot(p, Lik(p, x, r),'l')
abline(v = MLE(x,r))

rm(MaxFuncGS, x, n, r, p)

n=1; r=1; p=0.5 # extreme values
x <- rnbinom(n, r, 1 / (1 + p / r))
p <- seq(0,1,1e-2)
plot(p, Lik(p, x, r), 'l')

dev.off()


# testing gaussian hypergeometric function
rep <- 1000
max <- 500
a <- runif(rep, 1, max)
b <- runif(rep, 1, max)
c <- b + runif(rep, 1, max)
z <- rep(-1 / 1, rep)
f <- numeric(rep)
for(i in 1:rep){
  f[i] =  GaussHyp(a[i], b[i], c[i], z[i])
}
all(!is.nan(f) & !is.na(f)) # must be TRUE
rm(rep, max, a, b, c, z, f, i)

# testing Ixi
a = 1; b = 1 # uniform
xi <- seq(0, 1, 1e-2)
error <- numeric(length(xi))
for(k in 1:length(xi)){
  error[k] <- Ixi(a, b, xi[k]) - xi[k]
}
max(abs(error))

a = 1; b = 2 # left triangular
xi <- seq(0, 1, 1e-2)
error <- numeric(length(xi))
for(k in 1:length(xi)){
  error[k] <- Ixi(a, b, xi[k]) - (1 - (1 - xi[k])^2)
}
max(abs(error))

a = 2; b = 1 # right triangular
xi <- seq(0, 1, 1e-2)
error <- numeric(length(xi))
for(k in 1:length(xi)){
  error[k] <- Ixi(a, b, xi[k]) - (xi[k]^2)
}
max(abs(error))
rm(a, b, xi, error)

# testing bissection root-finding
rep = 10000
root <- runif(rep, 0, 1)
error <- numeric(rep)
f2GrauTeste <- function(x, root){
  # auxiliar root fixed at x = -10
  return(x^2 - x * (10 + root) + 10 * root)
}
for(i in 1: rep){
  error[i] <- BissectionRF(f2GrauTeste, 0, 1, root[i]) - root[i] 
}
all(is.finite(error)) # must be TRUE
max(abs(error)) # maximum error
rm(rep, root, error, f2GrauTeste, i)

# testing random number generator by acceptance-rejection method
dUnif <- function(p){
  return(1)
}
rep <- 10000
p <- seq(0, 1, by = 1e-3)
randNum <- RandNumAR(rep, dUnif)
f <- numeric(length(p))
for(i in 1:length(p)){
  f[i] <- dUnif(p[i])
}
hist(randNum, prob = TRUE,) # must be uniform
lines(p, f, lwd=3, col = 'red')
rm(randNum, f)
dev.off()

a = runif(3,1,5)
b = runif(3,1,5)
r = 1
randNum <- RandNumAR(rep, BetaPrior,r, a[1],b[1])
hist(randNum, prob = TRUE, ylim=c(0,4))
lines(p, BetaPrior(p, r, a[1], b[1]), lwd=3, col = 'red')
rm(randNum)
dev.off()

randNum <- RandNumAR(rep, BetaPrior,r, a[2],b[2])
hist(randNum, prob = TRUE, ylim = c(0,4))
lines(p, BetaPrior(p,r, a[2], b[2]), lwd=3, col = 'red')
rm(randNum)
dev.off()

randNum <- RandNumAR(rep, BetaPrior,r, a[3],b[3])
hist(randNum, prob = TRUE, ylim=c(0,4))
lines(p, BetaPrior(p,r, a[3], b[3]), lwd=3, col = 'red')
rm(randNum)
dev.off()

rm(rep, randNum,p, a, b)

# testing beta prior
p <- seq(0, 1, 1e-3)
a <- 1; b <- 1
r <- 1
plot(p, BetaPrior(p,r,a,b),'l') # uniform
integrate(BetaPrior, 0, 1,r, a, b)[[1]]
a <- 1; b <- 2
plot(p, BetaPrior(p,r,a,b),'l') # left triangular
integrate(BetaPrior, 0, 1,r, a, b)[[1]]
a <- 2; b <- 1
plot(p, BetaPrior(p,r,a,b),'l') # right triangular
integrate(BetaPrior, 0, 1,r, a, b)[[1]]
a <- 2; b <- 2
plot(p, BetaPrior(p,r,a,b),'l') # parabola
integrate(BetaPrior, 0, 1,r, a, b)[[1]]
a <- 1; b <- 3
plot(p, BetaPrior(p,r,a,b),'l') # left curve
integrate(BetaPrior, 0, 1,r, a, b)[[1]]
a <- 3; b <- 1
plot(p, BetaPrior(p,r,a,b),'l') # right curve
integrate(BetaPrior, 0, 1, r, a, b)[[1]]
dev.off()
rm(p, a, b)

# testing beta posterior
Literal_BetaPosterior <- function(p, x, r, a, b){
  AuxF <- function(p){
    return(Lik(p, x, r) * BetaPrior(p,r, a, b))
  }
  return(AuxF(p) / integrate(AuxF, 0, 1)[[1]])
}
n = 50; r = 1; p = 0.5
a = 1; b = 2
x <- rnbinom(n, r, 1 / (1 + p / r))
integrate(BetaPosterior, 0, 1, x, r, a, b)[[1]] # must be 1
p <- seq(0, 1, 1e-3)
max(abs(BetaPosterior(p, x, r, a, b) - 
          Literal_BetaPosterior(p, x, r, a, b))) # absolute error
plot(p , BetaPosterior(p, x, r, a, b), 'l', col='blue')
lines(p, Literal_BetaPosterior(p, x, r, a, b), 'l', col='red')
dev.off()
rm(Literal_BetaPosterior, n, r, p, a, b, x)

# testing beta self estimator
n = 20; r = 1; p = 0.2
a = 1; b = 2
x <- rnbinom(n, r, 1 / (1 + p / r))
abs(BetaSelf(x, r, a, b) - MeanF(BetaPosterior, x, r, a, b))
n = 50; r = 1; p = 0.5
a = 2; b = 2
x <- rnbinom(n, r, 1 / (1 + p / r))
abs(BetaSelf(x, r, a, b) - MeanF(BetaPosterior, x, r, a, b))
n = 100; r = 4; p = 0.8
a = 3; b = 2
x <- rnbinom(n, r, 1 / (1 + p / r))
abs(BetaSelf(x, r, a, b) - MeanF(BetaPosterior, x, r, a, b))

rm(n, r, p, x)

# testing gauss hypergeometric prior
p <- seq(0, 1, 1e-3)
r <- 1

plot(NA, xlim = c(0, 1), ylim = c(0, 3)) # left curves
a <- 1; b <- 1; c <- 2
lines(p, GaussHypPrior(p, r, a, b, c),'l',col='black')
integrate(GaussHypPrior, 0, 1, r, a, b, c)[[1]]
a <- 1; b <- 1; c <- 3
lines(p, GaussHypPrior(p, r, a, b, c),'l',col='red')
integrate(GaussHypPrior, 0, 1, r, a, b, c)[[1]]
a <- 1; b <- 1; c <- 4
lines(p, GaussHypPrior(p, r, a, b, c),'l',col='blue')
integrate(GaussHypPrior, 0, 1, r, a, b, c)[[1]]
a <- 1; b <- 2; c <- 5
lines(p, GaussHypPrior(p, r, a, b, c),'l',col='green')
integrate(GaussHypPrior, 0, 1, r, a, b, c)[[1]]
a <- 1; b <- 2; c <- 6
lines(p, GaussHypPrior(p, r, a, b, c),'l',col='orange')
integrate(GaussHypPrior, 0, 1, r, a, b, c)[[1]]
dev.off()

plot(NA, xlim = c(0, 1), ylim = c(0, 3)) # right curves
a <- 1; b <- 2; c <- 3
lines(p, GaussHypPrior(p, r, a, b, c),'l', col = 'black')
integrate(GaussHypPrior, 0, 1, r, a, b, c)[[1]]
a <- 1; b <- 3; c <- 4
lines(p, GaussHypPrior(p, r, a, b, c),'l',col = 'red')
integrate(GaussHypPrior, 0, 1, r, a, b, c)[[1]]
a <- 1; b <- 3; c <- 5
lines(p, GaussHypPrior(p, r, a, b, c),'l', col = 'blue')
integrate(GaussHypPrior, 0, 1, r, a, b, c)[[1]]
a <- 1; b <- 4; c <- 6
lines(p, GaussHypPrior(p, r, a, b, c),'l', col = 'green')
integrate(GaussHypPrior, 0, 1, r, a, b, c)[[1]]
a <- 1; b <- 5; c <- 6
dev.off()
rm(p, a, b, c)

# testing gauss hypergeometric posterior
Literal_GaussHypPosterior <- function(p, x, r, a, b, c){
  AuxF <- function(p){
    return(Lik(p, x, r) * GaussHypPrior(p, r, a, b, c))
  }
  return(AuxF(p) / integrate(AuxF, 0, 1)[[1]])
}
n = 50; r = 1; p = 0.5
a = 1; b = 2; c = 3
x <- rnbinom(n, r, 1 / (1 + p / r))
integrate(GaussHypPosterior, 0, 1, x, r, a, b, c)[[1]] # must be 1
p <- seq(0, 1, 1e-3)
max(abs(GaussHypPosterior(p, x, r, a, b, c) - 
          Literal_GaussHypPosterior(p, x, r, a, b, c))) # absolute error
plot(p , GaussHypPosterior(p, x, r, a, b, c), 'l', col='blue')
lines(p, Literal_GaussHypPosterior(p, x, r, a, b, c), 'l', col='red')
dev.off()
rm(Literal_GaussHypPosterior, n, r, p, a, b, c, x)

# testing gauss hypergeometric self estimator
n = 20; r = 1; p = 0.2
a = 1; b = 2; c = 3
x <- rnbinom(n, r, 1 / (1 + p / r))
abs(GaussHypSelf(x, r, a, b, c) - MeanF(GaussHypPosterior, x, r, a, b, c))
n = 50; r = 2; p = 0.5
a = 1; b = 2; c = 4
x <- rnbinom(n, r, 1 / (1 + p / r))
abs(GaussHypSelf(x, r, a, b, c) - MeanF(GaussHypPosterior, x, r, a, b, c))
n = 100; r = 4; p = 0.8
a = 1; b = 3; c = 5
x <- rnbinom(n, r, 1 / (1 + p / r))
abs(GaussHypSelf(x, r, a, b, c) - MeanF(GaussHypPosterior, x, r, a, b, c))
n = 1; r = 1; p = 0.8
a = 1; b = 3; c = 5
x <- rnbinom(n, r, 1 / (1 + p / r))
abs(GaussHypSelf(x, r, a, b, c) - MeanF(GaussHypPosterior, x, r, a, b, c))

rm(n, r, p, x, a, b, c)

# testing jeffeys prior
p <- seq(1e-3, 1, 1e-3)
r = 1
plot(NA, xlim = c(0, 1), ylim = c(0, 4))
lines(p, JefPrior(p, r) / integrate(JefPrior, 0, 1, r)[[1]],'l', col = 'black')
integrate(JefPrior, 0, 1, r)[[1]] # not normalized
r = 2
lines(p, JefPrior(p, r) / integrate(JefPrior, 0, 1, r)[[1]],'l', col = 'blue')
integrate(JefPrior, 0, 1, r)[[1]] # not normalized
r = 3
lines(p, JefPrior(p, r) / integrate(JefPrior, 0, 1, r)[[1]],'l', col = 'red')
integrate(JefPrior, 0, 1, r)[[1]] # not normalized
dev.off()
rm(p, r)

# testing jeffreys posterior
Literal_JefPosterior <- function(p, x, r){
  AuxF <- function(p){
    return(Lik(p, x, r) * JefPrior(p,r))
  }
  return(AuxF(p) / integrate(AuxF, 0, 1)[[1]])
}
n = 50; r = 1; p = 0.5
a = 1; b = 2
x <- rnbinom(n, r, 1 / (1 + p / r))
integrate(JefPosterior, 0, 1, x, r)[[1]] # must be 1
p <- seq(1e-3, 1, 1e-3)
max(abs(JefPosterior(p, x, r) - 
          Literal_JefPosterior(p, x, r))) # absolute error
plot(p, JefPosterior(p, x, r), 'l', col='blue')
lines(p, Literal_JefPosterior(p, x, r), 'l', col='red')
dev.off()
n = 100; r = 8; p = 0.8
a = 2; b = 2
x <- rnbinom(n, r, 1 / (1 + p / r))
integrate(JefPosterior, 0, 1, x, r)[[1]] # must be 1
p <- seq(1e-3, 1, 1e-3)
max(abs(JefPosterior(p, x, r) - 
          Literal_JefPosterior(p, x, r))) # absolute error
plot(p, JefPosterior(p, x, r), 'l', col='blue')
lines(p, Literal_JefPosterior(p, x, r), 'l', col='red')
dev.off()
rm(Literal_BetaPosterior, n, r, p, x)

# testing jeffreys self estimator
n = 20; r = 1; p = 0.2
x <- rnbinom(n, r, 1 / (1 + p / r))
abs(JefSelf(x, r) - MeanF(JefPosterior, x, r)) # must be 0
n = 50; r = 2; p = 0.5
x <- rnbinom(n, r, 1 / (1 + p / r))
abs(JefSelf(x, r) - MeanF(JefPosterior, x, r)) # must be 0
n = 100; r = 4; p = 0.8
x <- rnbinom(n, r, 1 / (1 + p / r))
abs(JefSelf(x, r) - MeanF(JefPosterior, x, r)) # must be 0
n = 1; r = 1; p = 0.8 # extreme values
x <- rnbinom(n, r, 1 / (1 + p / r))
MeanF(JefPosterior, x, r)
abs(JefSelf(x, r) - MeanF(JefPosterior, x, r)) # the estimator is undefined
rm(n, r, p, x)

# testing new jeffreys self estimator
set.seed(2024)
ns = c(10,20,50,100,200)
rs = c(1,2,3)
ps = c(0.01, 0.1,0.5,0.7,0.9,0.99)
i <- 1
err <- numeric(length=length(ns)*length(rs)*length(ps))
for (n in ns) {
  for (r in rs) {
    for (p in ps) {
      x <- rnbinom(n, r, 1 / (1 + p / r))
      err[i] <- abs(JefSelf2(x, r) - MeanF(JefPosterior, x, r))
      i+1
    }
  }
}
max(err) # must be 0
rm(ns, rs, ps, i, err, x)


# testing ETI
rep <- 1000
size <- 20
r <- 1
p <- runif(rep, 0, 1)
a <- runif(rep, 1, 5)
b <- runif(rep, 1, 5)
lInf <- numeric(rep)
lSup <- numeric(rep)
cover <- numeric(rep)
for(i in 1:rep){
  x <- rnbinom(size, r, 1 / (1 + p / r))
  interv <- ETI(BetaPosterior, x, r, a[i], b[i])
  lInf[i] <- interv$lInf
  lSup[i] <- interv$lSup
  cover[i] <- interv$cover
}
all(!is.nan(cover)) # must be true
all(!is.nan(lInf))  # must be true
all(!is.nan(lSup))  # must be true
c(mean = mean(lInf), min = min(lInf), max = max(lInf)) # lInf
c(mean = mean(lSup), min = min(lSup), max = max(lSup)) # lMax
c(mean = mean(cover), min = min(cover), max = max(cover)) # coverage
plot(cover,ylim=c(0.94,0.96),pch=16,cex=0.1)
dev.off()
rm(rep, size, r, p, a, b, lInf, lSup, cover, x, interv)

rep <- 1000
size <- 50
r <- 2
p <- runif(rep, 0, 1)
a <- runif(rep, 1, 5)
b <- runif(rep, 1, 5)
lInf <- numeric(rep)
lSup <- numeric(rep)
cover <- numeric(rep)
for(i in 1:rep){
  x <- rnbinom(size, r, 1 / (1 + p / r))
  interv <- ETI(BetaPosterior, x, r, a[i], b[i])
  lInf[i] <- interv$lInf
  lSup[i] <- interv$lSup
  cover[i] <- interv$cover
}
all(!is.nan(cover)) # must be true
all(!is.nan(lInf))  # must be true
all(!is.nan(lSup))  # must be true
c(mean = mean(lInf), min = min(lInf), max = max(lInf)) # lInf
c(mean = mean(lSup), min = min(lSup), max = max(lSup)) # lMax
c(mean = mean(cover), min = min(cover), max = max(cover)) # coverage
plot(cover,ylim=c(0.94,0.96),pch=16,cex=0.1)
dev.off()
rm(rep, size, r, p, a, b, lInf, lSup, cover, x, interv)

rep <- 1000
size <- 200
r <- 4
p <- runif(rep, 0, 1)
a <- runif(rep, 1, 5)
b <- runif(rep, 1, 5)
lInf <- numeric(rep)
lSup <- numeric(rep)
cover <- numeric(rep)
for(i in 1:rep){
  x <- rnbinom(size, r, 1 / (1 + p / r))
  interv <- ETI(BetaPosterior, x, r, a[i], b[i])
  lInf[i] <- interv$lInf
  lSup[i] <- interv$lSup
  cover[i] <- interv$cover
}
all(!is.nan(cover)) # must be true
all(!is.nan(lInf))  # must be true
all(!is.nan(lSup))  # must be true
c(mean = mean(lInf), min = min(lInf), max = max(lInf)) # lInf
c(mean = mean(lSup), min = min(lSup), max = max(lSup)) # lMax
c(mean = mean(cover), min = min(cover), max = max(cover)) # coverage
plot(cover,ylim=c(0.94,0.96),pch=16,cex=0.1)
dev.off()
rm(rep, size, r, p, a, b, lInf, lSup, cover, x, interv)


# testing SCI
DELTA <- 0.10
rep <- 1000
size <- 20
r <- 1
p <- runif(rep, 0, 1)
a <- runif(rep, 1, 5)
b <- runif(rep, 1, 5)
lInf <- numeric(rep)
lSup <- numeric(rep)
cover <- numeric(rep)
for(i in 1:rep){
  x <- rnbinom(size, r, 1 / (1 + p / r))
  interv <- SCI(DELTA, BetaSelf, BetaPosterior, x, r, a[i], b[i])
  lInf[i] <- interv$lInf
  lSup[i] <- interv$lSup
  cover[i] <- interv$cover
}
all(!is.nan(cover)) # must be true
all(!is.nan(lInf))  # must be true
all(!is.nan(lSup))  # must be true
c(mean = mean(lInf), min = min(lInf), max = max(lInf)) # lInf
c(mean = mean(lSup), min = min(lSup), max = max(lSup)) # lMax
c(mean = mean(cover), min = min(cover), max = max(cover)) # coverage
plot(seq(1,rep), cover, pch = 16, ylim = c(0, 1), cex = 0.1, col = 'red')
rm(DELTA, rep, size, r, p, a, b, lInf, lSup, cover, x, interv)

DELTA <- 0.10
rep <- 1000
size <- 50
r <- 1
p <- runif(rep, 0, 1)
a <- runif(rep, 1, 5)
b <- runif(rep, 1, 5)
lInf <- numeric(rep)
lSup <- numeric(rep)
cover <- numeric(rep)
for(i in 1:rep){
  x <- rnbinom(size, r, 1 / (1 + p / r))
  interv <- SCI(DELTA, BetaSelf, BetaPosterior, x, r, a[i], b[i])
  lInf[i] <- interv$lInf
  lSup[i] <- interv$lSup
  cover[i] <- interv$cover
}
all(!is.nan(cover)) # must be true
all(!is.nan(lInf))  # must be true
all(!is.nan(lSup))  # must be true
c(mean = mean(lInf), min = min(lInf), max = max(lInf)) # lInf
c(mean = mean(lSup), min = min(lSup), max = max(lSup)) # lMax
c(mean = mean(cover), min = min(cover), max = max(cover)) # coverage
points(seq(1,rep), cover, pch = 16,cex = 0.1, col = 'blue')
rm(DELTA, rep, size, r, p, a, b, lInf, lSup, cover, x, interv)

DELTA <- 0.10
rep <- 1000
size <- 100
r <- 1
p <- runif(rep, 0, 1)
a <- runif(rep, 1, 5)
b <- runif(rep, 1, 5)
lInf <- numeric(rep)
lSup <- numeric(rep)
cover <- numeric(rep)
for(i in 1:rep){
  x <- rnbinom(size, r, 1 / (1 + p / r))
  interv <- SCI(DELTA, BetaSelf, BetaPosterior, x, r, a[i], b[i])
  lInf[i] <- interv$lInf
  lSup[i] <- interv$lSup
  cover[i] <- interv$cover
}
all(!is.nan(cover)) # must be true
all(!is.nan(lInf))  # must be true
all(!is.nan(lSup))  # must be true
c(mean = mean(lInf), min = min(lInf), max = max(lInf)) # lInf
c(mean = mean(lSup), min = min(lSup), max = max(lSup)) # lMax
c(mean = mean(cover), min = min(cover), max = max(cover)) # coverage
points(seq(1, rep), cover, pch = 16,cex = 0.1, col = 'green')
legend('bottomright',legend = c('size = 20', 'size = 50', 'size = 100'),
       col=c('red','blue','green'),pch=c(16,16,16))
rm(DELTA, rep, size, r, p, a, b, lInf, lSup, cover, x, interv)
dev.off()

DELTA <- 0.05
rep <- 1000
size <- 200
r <- 4
p <- runif(rep, 0, 1)
a <- runif(rep, 1, 5)
b <- runif(rep, 1, 5)
lInf <- numeric(rep)
lSup <- numeric(rep)
cover <- numeric(rep)
for(i in 1:rep){
  x <- rnbinom(size, r, 1 / (1 + p / r))
  interv <- SCI(DELTA, BetaSelf, BetaPosterior, x, r, a[i], b[i])
  lInf[i] <- interv$lInf
  lSup[i] <- interv$lSup
  cover[i] <- interv$cover
}
all(!is.nan(cover)) # must be true
all(!is.nan(lInf))  # must be true
all(!is.nan(lSup))  # must be true
c(mean = mean(lInf), min = min(lInf), max = max(lInf)) # lInf
c(mean = mean(lSup), min = min(lSup), max = max(lSup)) # lMax
c(mean = mean(cover), min = min(cover), max = max(cover)) # coverage
plot(seq(1, rep), cover, ylim = c(0.7, 1), pch = 16, cex = 0.1, col = 'red')
rm(DELTA,rep, size, r, p, a, b, lInf, lSup, cover, x, interv)

DELTA <- 0.15
rep <- 1000
size <- 200
r <- 4
p <- runif(rep, 0, 1)
a <- runif(rep, 1, 5)
b <- runif(rep, 1, 5)
lInf <- numeric(rep)
lSup <- numeric(rep)
cover <- numeric(rep)
for(i in 1:rep){
  x <- rnbinom(size, r, 1 / (1 + p / r))
  interv <- SCI(DELTA, BetaSelf, BetaPosterior, x, r, a[i], b[i])
  lInf[i] <- interv$lInf
  lSup[i] <- interv$lSup
  cover[i] <- interv$cover
}
all(!is.nan(cover)) # must be true
all(!is.nan(lInf))  # must be true
all(!is.nan(lSup))  # must be true
c(mean = mean(lInf), min = min(lInf), max = max(lInf)) # lInf
c(mean = mean(lSup), min = min(lSup), max = max(lSup)) # lMax
c(mean = mean(cover), min = min(cover), max = max(cover)) # coverage
points(seq(1, rep), cover, pch = 16,cex = 0.1, col = 'blue')
legend('bottomright',legend = c('DELTA = 0.05', 'DELTA = 0.15'),
       col=c('red','blue'),pch=c(16,16))
rm(DELTA,rep, size, r, p, a, b, lInf, lSup, cover, x, interv)
dev.off()

# testing average coverage criterion
rep = 10000
delta = 0.10
r = 1
# a = runif(1, 1, 5);  b = runif(1, 1, 5)
a = 2; b = 1
p <- RandNumAR(rep, BetaPrior, r, a, b)
cover <- numeric(rep)
start = Sys.time()
nEst <- ACC(delta, r, BetaPrior, BetaPosterior, BetaSelf, a, b)
end = Sys.time()
cat('elapsed time: ', end - start, ' s')
for(i in 1:rep){
  x <- rnbinom(nEst, r, 1 / (1 + p[i] / r))
  cover[i] = SCI(delta, BetaSelf, BetaPosterior, x, r, a, b)$cover
}
mean(cover)

rep = 10000
delta = 0.10
r = 2
a = 2
b = 1
p <- RandNumAR(rep, BetaPrior, r, a, b)
cover <- numeric(rep)
start = Sys.time()
nEst <- ACC(delta, r, BetaPrior, BetaPosterior, BetaSelf, a, b)
end = Sys.time()
cat('elapsed time: ', end - start, ' s')
for(i in 1:rep){
  x <- rnbinom(nEst, r, 1 / (1 + p[i] / r))
  cover[i] = SCI(delta, BetaSelf, BetaPosterior, x, r, a, b)$cover
}
mean(cover)

rep = 10000
delta = 0.15
r = 4
a = 1
b = 2
p <- RandNumAR(rep, BetaPrior, r, a, b)
cover <- numeric(rep)
start = Sys.time()
nEst <- ACC(delta, r, BetaPrior, BetaPosterior, BetaSelf, a, b)
end = Sys.time()
cat('elapsed time: ', end - start, ' s')
for(i in 1:rep){
  x <- rnbinom(nEst, r, 1 / (1 + p[i] / r))
  cover[i] = SCI(delta, BetaSelf, BetaPosterior, x, r, a, b)$cover
}
mean(cover)
rm(start, end, rep, delta, r, a, b, p, cover, nEst, x)

# testing average length criterion
rep = 10000
w = 0.16
r = 1
# a = runif(1, 1, 5);  b = runif(1, 1, 5)
a = 2; b = 1
p <- RandNumAR(rep, BetaPrior, r, a, b)
width <- numeric(rep)
start = Sys.time()
nEst <- ALC(w, r, BetaPrior, BetaPosterior, a, b)
end = Sys.time()
cat('elapsed time: ', end - start, ' s')
for(i in 1:rep){
  x <- rnbinom(nEst, r, 1 / (1 + p[i] / r))
  interv = ETI(BetaPosterior, x, r, a, b)
  width[i] <- interv$lSup - interv$lInf
}
mean(width)

rep = 10000
w = 0.20
r = 2
a = 2
b = 1
p <- RandNumAR(rep, BetaPrior, r, a, b)
width <- numeric(rep)
start = Sys.time()
nEst <- ALC(w, r, BetaPrior, BetaPosterior, a, b)
end = Sys.time()
cat('elapsed time: ', end - start, ' s')
for(i in 1:rep){
  x <- rnbinom(nEst, r, 1 / (1 + p[i] / r))
  interv = ETI(BetaPosterior, x, r, a, b)
  width[i] <- interv$lSup - interv$lInf
}
mean(width)

rep = 10000
w = 0.30
r = 4
a = 1
b = 2
p <- RandNumAR(rep, BetaPrior, r, a, b)
width <- numeric(rep)
start = Sys.time()
nEst <- ALC(w, r, BetaPrior, BetaPosterior, a, b)
end = Sys.time()
cat('elapsed time: ', end - start, ' s')
for(i in 1:rep){
  x <- rnbinom(nEst, r, 1 / (1 + p[i] / r))
  interv = ETI(BetaPosterior, x, r, a, b)
  width[i] <- interv$lSup - interv$lInf
}
mean(width)
rm(start, end, rep, interv, width, r, a, b, p, nEst, x)

rep = 10000
w = 0.20
r = 2
a = 2
b = 1
p <- RandNumAR(rep, BetaPrior, r, a, b)
width <- numeric(rep)
start = Sys.time()
nEst <- ALC(w, r, BetaPrior, BetaPosterior, a, b)
end = Sys.time()
cat('elapsed time: ', end - start, ' s')
for(i in 1:rep){
  x <- rnbinom(nEst, r, 1 / (1 + p[i] / r))
  interv = ETI(BetaPosterior, x, r, a, b)
  width[i] <- interv$lSup - interv$lInf
}
mean(width)

rep = 10000
w = 0.30
r = 4
a = 1
b = 2
p <- RandNumAR(rep, BetaPrior, r, a, b)
width <- numeric(rep)
start = Sys.time()
nEst <- ALC(w, r, BetaPrior, BetaPosterior, a, b)
end = Sys.time()
cat('elapsed time: ', end - start, ' s')
for(i in 1:rep){
  x <- rnbinom(nEst, r, 1 / (1 + p[i] / r))
  interv = ETI(BetaPosterior, x, r, a, b)
  width[i] <- interv$lSup - interv$lInf
}
mean(width)
rm(start, end, rep, interv, width, r, a, b, p, nEst, x)
