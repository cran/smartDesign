## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, tidy.opts=list(width.cutoff=80), comment=NA)
require(smartDesign)

## ----smartDTR13---------------------------------------------------------------
mumat13 <- cbind(G1=c(30,25), G0=c(20,22))
mumat13
varmat13 <- cbind(G1=c(16,16),G0=c(16,16))
varmat13
dtr13 <- smartDTR(mu_Barm=mumat13, sigsq_Barm=varmat13,
                 Barm=c(1,3), nsubject=252, pG_A1=0.8, pran_Barm=c(.5,.5))
dtr13$dtrdat  # last row shoul dhave sigsq just larger than to input
dtr13$true_mumix
dtr13$true_sigmix

## ----plotmeanDTR13------------------------------------------------------------
plot.smartDTR(dtr13, metric="mean", xtype="spec", mar=c(4,4,4,6), legend.inset=c(-.2, 0))
plot(dtr13, metric="mean", relativeBias=TRUE, xtype="sens", mar=c(4,4,4,6), legend.inset=c(-.2, 0), ylim=c(0,10), xlab="Sens")


## ----plotDTRB-----------------------------------------------------------------
plot.smartDTR(dtr13, metric="variance", xtype="spec", mar=c(4,4,4,6), legend.inset=c(-.2, 0))
plot.smartDTR(dtr13, metric="variance", xtype="sens", relativeBias=TRUE)

## ----plotDTRbias--------------------------------------------------------------
plot.smartDTR(dtr13, metric="mean", relativeBias=TRUE, xtype="sens", mar=c(4,4,4,6), legend.inset=c(-.2, 0))
plot.smartDTR(dtr13, metric="variance", relativeBias=TRUE, xtype="sens", mar=c(4,4,4,6), legend.inset=c(-.2, 0))

## ----dtrpower-----------------------------------------------------------------
mumat13 <- cbind(G1=c(30,35), G0=c(20,28))
varmat13 <- cbind(G1=c(100,100),G0=c(100,100))
#mumat13 <- cbind(G1=c(30,25), G0=c(20,22))
#varmat13 <- cbind(G1=c(16,16),G0=c(16,16))

dtr13 <- smartDTR(mu_Barm=mumat13, sigsq_Barm=varmat13,
                 Barm=c(1,3), nsubject=252)
plot(dtr13, metric="mean", xtype="spec", relativeBias=TRUE)
plot(dtr13, metric="variance", xtype="sens", relativeBias=FALSE)
mumat24 <- cbind(G1=c(25,32), G0=c(18,23))
varmat24 <- cbind(G1=c(100,100),G0=c(100,100))
dtr24 <- smartDTR(mu_Barm=mumat24, sigsq_Barm=varmat24,
                 Barm=c(2,4), nsubject=252)
plot(dtr24)
pdtr13vs24 <- powerDTR(dtr13, dtr24)
names(pdtr13vs24)
#source("../R/powerDTR.R")
plot(pdtr13vs24,  mar=c(4,4,4,6), legend.inset=c(-.2, 0), hline=0.8)
plot(pdtr13vs24, alpha=0.05, xtype="sens", cex.lab=.5, cex.axis=.5)

## ----dtrpower2----------------------------------------------------------------
mumat13 <- cbind(G1=c(30,35), G0=c(20,10))
varmat13 <- cbind(G1=c(100,100),G0=c(100,100))
dtr13 <- smartDTR(mu_Barm=mumat13, sigsq_Barm=varmat13,
                 Barm=c(1,3), nsubject=252)

mumat57 <- cbind(G1=c(29,32), G0=c(18,9))
varmat57 <- cbind(G1=c(100,100),G0=c(100,100))
dtr57 <- smartDTR(mu_Barm=mumat57, sigsq_Barm=varmat57,
                 Barm=c(5,7), nsubject=252)

pdtr13vs57 <- powerDTR(dtr13, dtr57)
names(pdtr13vs57)

## Plot DTR Power

plot(pdtr13vs57, mar=c(4,4,4,6), cex.axis=.5, cex.lab=.5, legend.inset=c(-.2, 0), hline=0.8)

