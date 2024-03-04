## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, tidy.opts=list(width.cutoff=80), comment=NA)
require(smartDesign)

## ----smartSST-----------------------------------------------------------------
sst1 <- smartSST(mu_Barm=c(G1=30, G0=20), sigsq_Barm=c(G1=16,G0=16),
        Barm=1, sens=seq(.5, 1, by=.1),  spec=seq(.5, 1, by=.1),
        nsubject=252, pG_A1=0.8, pran_A1 = 0.5, pran_Barm = 0.5)
names(sst1$sstdat)
sst1

## ----sst3---------------------------------------------------------------------
sst3 <- smartSST(mu_Barm=c(G1=25, G0=20), sigsq_Barm=c(G1=16,G0=16), Barm=3, sens=seq(.5, 1, by=.1), nsubject=252, spec=seq(.5, 1, by=.1), pG_A1=0.8)
print(sst3, digits=3)

## ----b1means------------------------------------------------------------------
#plot(sst1, metric="mean", xtype="spec",cex.lab=1.5, cex.axis=1.2, 
#     mar=c(5,5,2,9), legend.inset=c(-0.45,0.1))
plot.smartSST(sst1, metric="mean", xtype="spec", cex.lab=1.5, cex.axis=1.2, 
     mar=c(5,5,2,9), legend.inset=c(-0.45,0.1), relativeBias=FALSE)



## ----b1variance---------------------------------------------------------------
plot.smartSST(sst1, metric="variance", xtype="spec", cex.lab=1.5, cex.axis=1.2, 
     mar=c(5,5,2,9), legend.inset=c(-0.45,0.1))

## ----b3meanplots--------------------------------------------------------------
plot.smartSST(sst3, metric="mean", xtype="spec", ylim=c(-30,10),
mar=c(4,4,4,6), relativeBias=TRUE)

## ----b3varplots---------------------------------------------------------------
plot.smartSST(sst3, metric="variance", xtype="sens",
  mar=c(4,4,4,6), cex.axis=.5,  cex.lab=.5)

## ----b1bias-------------------------------------------------------------------
plot.smartSST(sst3, metric="mean", xtype="sens", cex.axis=.5,cex.lab=.5, mar=c(4,4,4,6), relativeBias=TRUE)

## ----powerspec----------------------------------------------------------------
## need to add warnigns that power does not work if sens/spec are different
sst1 <- smartSST(mu_Barm=c(G1=30, G0=20), sigsq_Barm=c(G1=100,G0=100), Barm=1, sens=seq(.5, 1, by=.1), nsubject=316, spec=seq(.5, 1, by=.1), pG_A1=0.8)
sst2 <- smartSST(mu_Barm=c(G1=25, G0=11), sigsq_Barm=c(G1=100,G0=100), Barm=2, sens=seq(.5, 1, by=.1), nsubject=316, spec=seq(.5, 1, by=.1), pG_A1=0.8)

powsst1vs2 <- powerSST(sst1, sst2, alpha=0.05)
## access the full power data.frame/matrix
powsst1vs2$powerdat ## same as: print.powerSST(powsst1vs2)
plot.powerSST(powsst1vs2, xtype="sens", hline=0.8, ylim=c(0,1), mar=c(4,4,4,6), cex.axis=.5, cex.lab=.5) 
plot.powerSST(powsst1vs2, xtype="sens", alpha=0.10, ylim=c(.5, 1), hline=.9, mar=c(4,4,4,6), cex.axis=.5,cex.lab=.5)

## ----powerspec2---------------------------------------------------------------
## need to add warnigns that power does not work if sens/spec are different
sst3 <- smartSST(mu_Barm=c(G1=35, G0=28), sigsq_Barm=c(G1=100,G0=100), Barm=3, sens=seq(.5, 1, by=.1), nsubject=1256, spec=seq(.5, 1, by=.1), pG_A1=0.8)
sst4 <- smartSST(mu_Barm=c(G1=30, G0=23), sigsq_Barm=c(G1=100,G0=100), Barm=4, sens=seq(.5, 1, by=.1), nsubject=1256, spec=seq(.5, 1, by=.1), pG_A1=0.8)

powsst3vs4 <- powerSST(sst3, sst4, alpha=0.05)
## access the full power data.frame/matrix
powsst3vs4$powerdat ## same as: print.powerSST(powsst3vs4)
plot.powerSST(powsst3vs4, xtype="sens", hline=0.8, ylim=c(0,1), cex.axis=.5, cex.lab=.5, mar=c(4,4,4,6), legend.inset=c(-.2,0))

