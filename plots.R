#-------------------------------------------------------------------------------
# A Bayesian Approach to Determine Sample Size in Erlang Single Server Queues
# Plotting charts
#-------------------------------------------------------------------------------
#
# PROGRAMMED BY:
# Eriky S. Gomes & Frederico R. B. Cruz
# E-mail: eriky-tn@ufmg.br
# E-mail: fcruz@est.ufmg.br
#
# AFFILIATION: 
# Universidade Federal de Minas Gerais
#
# @2023 Gomes & Cruz
# v.2023.12.18
#
#-------------------------------------------------------------------------------
# Directory
#-------------------------------------------------------------------------------
setwd('C:/Users/eriky/Box/Projetos/5MM1/Er-1 ss/results 2023-12-20')
#-------------------------------------------------------------------------------
# Packages
#-------------------------------------------------------------------------------

library(ggplot2) # for plotting
library(dplyr)  # for analysis
library(reshape2) # for analysys
library(Cairo)

#-------------------------------------------------------------------------------
# basic plotting layout 
#-------------------------------------------------------------------------------

# some manual customization variables ------------------------------------------

#cpManual<-c('black','blue','red2','green4','orange3','purple2','navyblue','cyan3',
#            'salmon4','gold','violet','limegreen','springgreen','slateblue4')
cpManual <- c('black','red2', 'blue2', 'green4')
#ltManual<-c('solid','longdash','dashed','dotted','dotdash')
ltManual <- c('solid', 'longdash', 'dashed', 'dotted')
shapeManual <- c(15, 16, 17, 18, 4, 8)

# basic plot -------------------------------------------------------------------
ggBase <-
  ggplot() + theme_minimal() + labs(color = '',
                                    linetype = '',
                                    shape = '') +
  theme(
    #panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.7, linetype='solid', color='gray80'),
    panel.grid.minor = element_line(linewidth = 0.2, linetype='solid', color='gray95'),
    panel.background = element_rect(fill = 'white', color = 'black'),
    axis.ticks = element_line(),
    axis.ticks.length = unit(c(-0, 1), 'mm'),
    #legend.box.spacing = unit(0,'pt'),
    legend.margin = margin(
      t = 1,
      r = 1,
      b = 1,
      l = 1,
      unit = 'mm'
    ),
    plot.margin = margin(10, 10, 10, 10, "mm"),
    panel.border = element_rect(fill = NA, linetype = "solid", color = "black")
  ) +
  scale_linetype_manual(values = ltManual) +
  scale_color_manual(values = cpManual) +
  scale_shape_manual(values = shapeManual)

rm(cpManual, ltManual, shapeManual)
# ------------------------------------------------------------------------------
# Sample size plotting function
# ------------------------------------------------------------------------------

# plotting function for ACC
PlotSS_ACC <- function(tabPlot, namePlot){
  myPlot <- ggBase + 
    xlab('size') + 
    ylab('average coverage') +
    geom_line(data = tabPlot, 
              aes(x = size, y = value,
                  color = variable, linetype = variable)) +
    geom_point(data = tabPlot, 
               aes(x = size, y = value, 
                   color = variable, shape = variable), size = 1.5) +
    theme(legend.position = c(0.98, 0), legend.justification = c(1, 0)) +
    scale_x_continuous(minor_breaks = seq(0, 500, 10)) +
    scale_y_continuous(minor_breaks = seq(0, 1, 0.05))
  
  #myPlot
  ggsave(
    file = paste0(namePlot),
    plot = myPlot,
    width = 6,
    height = 4,
    device = cairo_pdf
  )
}

# plotting function for ALC
PlotSS_ALC <- function(tabPlot, namePlot){
  myPlot <- ggBase + 
    xlab('size') + 
    ylab('average length') +
    geom_line(data = tabPlot, 
              aes(x = size, y = value,
                  color = variable, linetype = variable)) +
    geom_point(data = tabPlot, 
               aes(x = size, y = value, 
                   color = variable, shape = variable), size = 1.5) +
    theme(legend.position = c(1.0, 1.1), legend.justification = c(1, 1)) +
    scale_x_continuous(minor_breaks = seq(0, 500, 10)) +
    scale_y_continuous(minor_breaks = seq(0, 1, 0.05))
  
  #myPlot
  ggsave(
    file = paste0(namePlot),
    plot = myPlot,
    width = 6,
    height = 4,
    device = cairo_pdf
  )
}

#-------------------------------------------------------------------------------
# Average coverage criterion (SCI)
# Abacuses for r = 2
#-------------------------------------------------------------------------------

size <- seq(10,450,10)
DELTA <- c(0.08, 0.10, 0.12, 0.15)

# Beta -------------------------------------------------------------------------
a <- c(1, 1, 2, 2)
b <- c(1, 2, 1, 2)

tabBeta <- read.table('tabBeta_SCI_2.txt', header=TRUE)
col <- paste0('B(',rep(a, each = 4),',', rep(b, each = 4), ') ', '\u03B4 = ' ,
              rep(DELTA, times = 4))
colnames(tabBeta) = col
tabBeta[,] <- lapply(tabBeta[,], function(x) as.numeric(x))
tabBeta <- cbind(size = as.numeric(row.names(tabBeta)), tabBeta)


tabBeta1 <- melt(tabBeta[, c(1, 2:5)], id.vars = 'size')
tabBeta2 <- melt(tabBeta[, c(1, 6:9)], id.vars = 'size')
tabBeta3 <- melt(tabBeta[, c(1, 10:13)], id.vars = 'size')
tabBeta4 <- melt(tabBeta[, c(1, 14:17)], id.vars = 'size')

PlotSS_ACC(tabBeta1, 'abacBeta1_ACC_2.pdf')
PlotSS_ACC(tabBeta2, 'abacBeta2_ACC_2.pdf')
PlotSS_ACC(tabBeta3, 'abacBeta3_ACC_2.pdf')
PlotSS_ACC(tabBeta4, 'abacBeta4_ACC_2.pdf')

rm(a, b, tabBeta, col, tabBeta1, tabBeta2, tabBeta3, tabBeta4)

# Gaussian ---------------------------------------------------------------------
a <- c(1, 1, 1, 1)
b <- c(1, 2, 3, 4)
c <- c(3, 6, 4, 6)

tabGau <- read.table('tabGau_SCI_2.txt', header=TRUE)
col <- paste0('G(',rep(a, each = 4),',', rep(b, each = 4),',', rep(c, each = 4),
              ') ', '\u03B4 = ', rep(DELTA, times = 4))
colnames(tabGau) = col
tabGau[,] <- lapply(tabGau[,], function(x) as.numeric(x))
tabGau <- cbind(size = as.numeric(row.names(tabGau)), tabGau)

tabGau1 <- melt(tabGau[, c(1, 2:5)], id.vars = 'size')
tabGau2 <- melt(tabGau[, c(1, 6:9)], id.vars = 'size')
tabGau3 <- melt(tabGau[, c(1, 10:13)], id.vars = 'size')
tabGau4 <- melt(tabGau[, c(1, 14:17)], id.vars = 'size')

PlotSS_ACC(tabGau1, 'abacGau1_ACC_2.pdf')
PlotSS_ACC(tabGau2, 'abacGau2_ACC_2.pdf')
PlotSS_ACC(tabGau3, 'abacGau3_ACC_2.pdf')
PlotSS_ACC(tabGau4, 'abacGau4_ACC_2.pdf')
rm(a, b, c, tabGau, col, tabGau1, tabGau2, tabGau3, tabGau4)

# Jeffreys ---------------------------------------------------------------------
tabJef <- read.table('tabJef_SCI_2.txt', header=TRUE)
col <- paste0('Jefreys',' \u03B4 = ', DELTA)
colnames(tabJef) = col
tabJef[,] <- lapply(tabJef[,], function(x) as.numeric(x))
tabJef <- cbind(size = as.numeric(row.names(tabJef)), tabJef)

tabJef1 <- melt(tabJef[, c(1, 2)], id.vars = 'size')
tabJef2 <- melt(tabJef[, c(1, 3)], id.vars = 'size')
tabJef3 <- melt(tabJef[, c(1, 4)], id.vars = 'size')
tabJef4 <- melt(tabJef[, c(1, 5)], id.vars = 'size')

PlotSS_ACC(tabJef1, 'abacJef1_ACC_2.pdf')
PlotSS_ACC(tabJef2, 'abacJef2_ACC_2.pdf')
PlotSS_ACC(tabJef3, 'abacJef3_ACC_2.pdf')
PlotSS_ACC(tabJef4, 'abacJef4_ACC_2.pdf')

rm(tabJef, col, tabJef1, tabJef2, tabJef3, tabJef4)


#-------------------------------------------------------------------------------
# Average coverage criterion (SCI)
# Abacuses for r = 4
#-------------------------------------------------------------------------------

size <- seq(10, 450, 10)
DELTA <- c(0.08, 0.10, 0.12, 0.15)

# Beta -------------------------------------------------------------------------
a <- c(1, 1, 2, 2)
b <- c(1, 2, 1, 2)

tabBeta <- read.table('tabBeta_SCI_4.txt', header=TRUE)
col <- paste0('B(',rep(a, each = 4),',', rep(b, each = 4), ') ', '\u03B4 = ' ,
              rep(DELTA, times = 4))
colnames(tabBeta) = col
tabBeta[,] <- lapply(tabBeta[,], function(x) as.numeric(x))
tabBeta <- cbind(size = as.numeric(row.names(tabBeta)), tabBeta)


tabBeta1 <- melt(tabBeta[, c(1, 2:5)], id.vars = 'size')
tabBeta2 <- melt(tabBeta[, c(1, 6:9)], id.vars = 'size')
tabBeta3 <- melt(tabBeta[, c(1, 10:13)], id.vars = 'size')
tabBeta4 <- melt(tabBeta[, c(1, 14:17)], id.vars = 'size')

PlotSS_ACC(tabBeta1, 'abacBeta1_ACC_4.pdf')
PlotSS_ACC(tabBeta2, 'abacBeta2_ACC_4.pdf')
PlotSS_ACC(tabBeta3, 'abacBeta3_ACC_4.pdf')
PlotSS_ACC(tabBeta4, 'abacBeta4_ACC_4.pdf')

rm(a,b, tabBeta, col, tabBeta1, tabBeta2, tabBeta3, tabBeta4)

# Gaussian ---------------------------------------------------------------------
a <- c(1, 1, 1, 1)
b <- c(1, 2, 3, 4)
c <- c(3, 6, 4, 6)

tabGau <- read.table('tabGau_SCI_4.txt', header=TRUE)
col <- paste0('G(',rep(a, each = 4),',', rep(b, each = 4),',', rep(c, each = 4),
              ') ', '\u03B4 = ', rep(DELTA, times = 4))
colnames(tabGau) = col
tabGau[,] <- lapply(tabGau[,], function(x) as.numeric(x))
tabGau <- cbind(size = as.numeric(row.names(tabGau)), tabGau)


tabGau1 <- melt(tabGau[, c(1, 2:5)], id.vars = 'size')
tabGau2 <- melt(tabGau[, c(1, 6:9)], id.vars = 'size')
tabGau3 <- melt(tabGau[, c(1, 10:13)], id.vars = 'size')
tabGau4 <- melt(tabGau[, c(1, 14:17)], id.vars = 'size')

PlotSS_ACC(tabGau1, 'abacGau1_ACC_4.pdf')
PlotSS_ACC(tabGau2, 'abacGau2_ACC_4.pdf')
PlotSS_ACC(tabGau3, 'abacGau3_ACC_4.pdf')
PlotSS_ACC(tabGau4, 'abacGau4_ACC_4.pdf')

rm(a, b, c, tabGau, col, tabGau1, tabGau2, tabGau3, tabGau4)

# Jeffreys ---------------------------------------------------------------------
tabJef <- read.table('tabJef_SCI_4.txt', header=TRUE)
col <- paste0('Jefreys',' \u03B4 = ', DELTA)
colnames(tabJef) = col
tabJef[,] <- lapply(tabJef[,], function(x) as.numeric(x))
tabJef <- cbind(size = as.numeric(row.names(tabJef)), tabJef)

tabJef1 <- melt(tabJef[, c(1, 2)], id.vars = 'size')
tabJef2 <- melt(tabJef[, c(1, 3)], id.vars = 'size')
tabJef3 <- melt(tabJef[, c(1, 4)], id.vars = 'size')
tabJef4 <- melt(tabJef[, c(1, 5)], id.vars = 'size')

PlotSS_ACC(tabJef1, 'abacJef1_ACC_4.pdf')
PlotSS_ACC(tabJef2, 'abacJef2_ACC_4.pdf')
PlotSS_ACC(tabJef3, 'abacJef3_ACC_4.pdf')
PlotSS_ACC(tabJef4, 'abacJef4_ACC_4.pdf')

rm(tabJef, col, tabJef1, tabJef2, tabJef3, tabJef4)

#-------------------------------------------------------------------------------
# Average length criterion (ETI)
# Abacuses for r = 2
#-------------------------------------------------------------------------------

size <- seq(10,450,10)
WIDTH <- c(0.16, 0.20, 0.24, 0.30)

# Beta -------------------------------------------------------------------------
a <- c(1, 1, 2, 2)
b <- c(1, 2, 1, 2)

tabBeta <- read.table('tabBeta_ETI_2.txt', header=TRUE)
col <- paste0('B(', a,',', b, ') ')
colnames(tabBeta) = col
tabBeta[,] <- lapply(tabBeta[,], function(x) as.numeric(x))
tabBeta <- cbind(size = as.numeric(row.names(tabBeta)), tabBeta)

tabBeta <- melt(tabBeta, id.vars = 'size')
PlotSS_ALC(tabBeta, 'abacBeta_ALC_2.pdf')

rm(a,b, tabBeta, col)

# Gaussian ---------------------------------------------------------------------
a <- c(1, 1, 1, 1)
b <- c(1, 2, 3, 4)
c <- c(3, 6, 4, 6)

tabGau <- read.table('tabGau_ETI_2.txt', header=TRUE)
col <- paste0('G(', a,',', b,',', c, ') ')
colnames(tabGau) = col
tabGau[,] <- lapply(tabGau[,], function(x) as.numeric(x))
tabGau <- cbind(size = as.numeric(row.names(tabGau)), tabGau)

tabGau <- melt(tabGau, id.vars = 'size')
PlotSS_ALC(tabGau, 'abacGau_ALC_2.pdf')

rm(a, b, c, tabGau, col)

# Jeffreys ---------------------------------------------------------------------
tabJef <- read.table('tabJef_ETI_2.txt', header=TRUE)
col <- paste0('Jefreys')
colnames(tabJef) = col
tabJef[,] <- as.numeric(tabJef[,])
tabJef <- cbind(size = as.numeric(row.names(tabJef)), tabJef)

tabJef <- melt(tabJef, id.vars = 'size')
PlotSS_ALC(tabJef, 'abacJef_ALC_2.pdf')

rm(tabJef, col)

#-------------------------------------------------------------------------------
# Average length criterion (ETI)
# Abacuses for r = 4
#-------------------------------------------------------------------------------

size <- seq(10,450,10)
WIDTH <- c(0.16, 0.20, 0.24, 0.30)

# Beta -------------------------------------------------------------------------
a <- c(1, 1, 2, 2)
b <- c(1, 2, 1, 2)

tabBeta <- read.table('tabBeta_ETI_4.txt', header=TRUE)
col <- paste0('B(', a,',', b, ') ')
colnames(tabBeta) = col
tabBeta[,] <- lapply(tabBeta[,], function(x) as.numeric(x))
tabBeta <- cbind(size = as.numeric(row.names(tabBeta)), tabBeta)

tabBeta <- melt(tabBeta, id.vars = 'size')
PlotSS_ALC(tabBeta, 'abacBeta_ALC_4.pdf')

rm(a,b, tabBeta, col)

# Gaussian ---------------------------------------------------------------------
a <- c(1, 1, 1, 1)
b <- c(1, 2, 3, 4)
c <- c(3, 6, 4, 6)

tabGau <- read.table('tabGau_ETI_4.txt', header=TRUE)
col <- paste0('G(', a,',', b,',', c, ') ')
colnames(tabGau) = col
tabGau[,] <- lapply(tabGau[,], function(x) as.numeric(x))
tabGau <- cbind(size = as.numeric(row.names(tabGau)), tabGau)

tabGau <- melt(tabGau, id.vars = 'size')
PlotSS_ALC(tabGau, 'abacGau_ALC_4.pdf')

rm(a, b, c, tabGau, col)

# Jeffreys ---------------------------------------------------------------------
tabJef <- read.table('tabJef_ETI_4.txt', header=TRUE)
col <- paste0('Jefreys')
colnames(tabJef) = col
tabJef[,] <- as.numeric(tabJef[,])
tabJef <- cbind(size = as.numeric(row.names(tabJef)), tabJef)

tabJef <- melt(tabJef, id.vars = 'size')
PlotSS_ALC(tabJef, 'abacJef_ALC_4.pdf')

rm(tabJef, col)

#-------------------------------------------------------------------------------
# Plotting prior distributions
#-------------------------------------------------------------------------------

PlotDensityRho <- function(tabPlot, namePlot, 
                           legend_pos = c(0, 0), legend_justif = c(0, 0),
                           nrowLegend = 2){
  myPlot <- ggBase + 
    xlab(expression(rho)) + # rho unicode
    ylab('density') +
    geom_line(data = tabPlot, 
              aes(x = p, y = value,
                  color = variable, linetype = variable)) +
    geom_point(data = filter(tabPlot, row_number() %% 3 == 1), # gap between points 
               aes(x = p, y = value, 
                   color = variable, shape = variable), size = 1.5) +
    theme(legend.position = legend_pos, 
          legend.justification = legend_justif,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    guides(linetype = guide_legend(nrow = nrowLegend))
    
  #myPlot
  ggsave(
    file = paste0(namePlot),
    plot = myPlot,
    width = 6,
    height = 4,
    device = cairo_pdf
  )
}

p <- seq(0, 1, 1e-2)

# Beta prior  ------------------------------------------------------------------
a <- c(1, 1, 2, 2)
b <- c(1, 2, 1, 2)
r = 2 # dummy

betap1 <- BetaPrior(p, r, a[1], b[1])
betap2 <- BetaPrior(p, r, a[2], b[2])
betap3 <- BetaPrior(p, r, a[3], b[3])
betap4 <- BetaPrior(p, r, a[4], b[4])

betap <- as.data.frame(cbind(betap1, betap2, betap3, betap4))
col <- paste0('B(', a,',', b, ') ')
colnames(betap) <- col
betap <- cbind(p = p, betap)
betap <- melt(betap, id.vars = 'p')

PlotDensityRho(betap, 'betaPrior.pdf', 
               legend_pos = c(0.5, 0), legend_justif = c(0.5, 0))
rm(a, b, r, betap1, betap2, betap3, betap4, betap, col)

# Gaussian prior ---------------------------------------------------------------
a <- c(1, 1, 1, 1)
b <- c(1, 2, 3, 4)
c <- c(3, 6, 4, 6)

# r = 2
r = 2

gaup1 <- GaussHypPrior(p, r, a[1], b[1], c[1])
gaup2 <- GaussHypPrior(p, r, a[2], b[2], c[2])
gaup3 <- GaussHypPrior(p, r, a[3], b[3], c[3])
gaup4 <- GaussHypPrior(p, r, a[4], b[4], c[4])

gaup <- as.data.frame(cbind(gaup1, gaup2, gaup3, gaup4))
col <- paste0('G(', a,',', b,',', c,',','r=', r,') ')
colnames(gaup) <- col
gaup <- cbind(p = p, gaup)
gaup <- melt(gaup, id.vars = 'p')

PlotDensityRho(gaup, 'gauPrior-r2.pdf', 
               legend_pos = c(0.5, 1.07), legend_justif = c(0.5, 1))
rm(r, gaup1, gaup2, gaup3, gaup4, col, gaup)

# r = 4
r = 4

gaup1 <- GaussHypPrior(p, r, a[1], b[1], c[1])
gaup2 <- GaussHypPrior(p, r, a[2], b[2], c[2])
gaup3 <- GaussHypPrior(p, r, a[3], b[3], c[3])
gaup4 <- GaussHypPrior(p, r, a[4], b[4], c[4])

gaup <- as.data.frame(cbind(gaup1, gaup2, gaup3, gaup4))
col <- paste0('G(', a,',', b,',', c,',','r=', r,') ')
colnames(gaup) <- col
gaup <- cbind(p = p, gaup)
gaup <- melt(gaup, id.vars = 'p')

PlotDensityRho(gaup, 'gauPrior-r4.pdf', 
               legend_pos = c(0.5, 1.07), legend_justif = c(0.5, 1))

rm(a, b, c, r, gaup1, gaup2, gaup3, gaup4, col, gaup)

# Jeffreys prior  --------------------------------------------------------------
tol = 1e-2
p <- seq(tol, 1 - tol, tol)
r = c(2, 4)

jefp1 <- JefPrior(p, r[1])
jefp2 <- JefPrior(p, r[2])

jefp <- as.data.frame(cbind(jefp1, jefp2))
col <- paste0('Jeffreys, r=', r)
colnames(jefp) <- col
jefp <- cbind(p = p, jefp)
jefp <- melt(jefp, id.vars = 'p')

PlotDensityRho(jefp, 'jefPrior.pdf', 
               legend_pos = c(0.98, 1.05), legend_justif = c(1, 1))

rm(jefp1, jefp2, jefp, col)

