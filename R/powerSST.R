
## Power of SST

powerSST <- function(sst1, sst2, pG_A1 = 0.8, pG_A2 = 0.8, alpha=0.05) {
    ## future argument for power of other stats: power="response", 
    ##     Part III: SST-- Plot power for A1 ##
    Blev1 <- sst1$Barm
    Blev2 <- sst2$Barm
   
    ##Generate mu's and sigsq's dataframes
    mu1_dat <- sst1$sstdat[,c("sens","spec","mu_Barm")]
                                  #cbind(sst1$sens.dat, sst1$spec.dat, sst1$mu_B_matrix)
    mu2_dat <- sst2$sstdat[,c("sens","spec","mu_Barm")]
                                   #cbind(sst2$sens.dat, sst2$spec.dat, sst2$mu_B_matrix)
    
    sigsq1_dat <- sst1$sstdat[,c("sens","spec","sigsq_Barm", "n_Barm")]
            #cbind(sst1$sens.dat, sst1$spec.dat, sst1$sigsq_B_matrix)
    sigsq2_dat <- sst2$sstdat[,c("sens","spec", "sigsq_Barm", "n_Barm")]
            #cbind(sst2$sens.dat, sst2$spec.dat, sst2$sigsq_B_matrix)

    zstat <- (mu1_dat[,"mu_Barm"] - mu2_dat[,"mu_Barm"]) /
        sqrt(sigsq1_dat[,"sigsq_Barm"]/sigsq1_dat$n_Barm +
             sigsq2_dat[,"sigsq_Barm"]/sigsq2_dat$n_Barm)
                                             
    pval <- 2*pnorm(-abs(zstat))
    ##p_B1B2<-2*pnorm(-abs(z_B1B2))       # p-value of B1 vs. B2
   
    cvalue <- round(qnorm(1-alpha/2), 2)  #critical value from normal dist, alpha level, 2-sided
    power_Barms <- 1 - pnorm(cvalue - zstat) + pnorm(-1*cvalue - zstat) 

    powerdf <- round(cbind.data.frame(mu1_dat, sigsq1_dat$sigsq_Barm, mu2_dat$mu_Barm, sigsq2_dat$sigsq_Barm, zstat, pval, power_Barms), 5)
    colnames(powerdf) <- c("sens", "spec", paste0("mu_B", Blev1), paste0("sigsq_B", Blev1), paste0("mu_B",Blev2), paste0("sigsq_B", Blev2),"z_Barm", "pval", "power")
#ztest_A1_dat_B1B2 <- round(cbind(pR_A1_dat$pie1_A1, pR_A1_dat$pie2_A1,
#                               mu_B1_dat$mu_B1,mu_B2_dat$mu_B2,z_B1B2,p_B1B2,power_B1B2),3)
#colnames(ztest_A1_dat_B1B2) <- c("pie1","pie2","mu_B1","mu_B2","z_B1B2","p_B1B2","power_B1B2")
    out <- list(powerdat=powerdf, Barm=c(Blev1, Blev2), alpha=alpha) 
    class(out) <- c("powerSST", "list")
    return(out)
}

print.powerSST <- function(x, ...) {
   cat("\n B-levels: ", paste(x$Barm, collapse=", "), "; alpha = ", x$alpha, "\n\n")   
   print(x$powerdat, ...)
   invisible(x$powerdat)
}
plot.powerSST <- function(x, xtype="spec",  metric=c("power", "zstat"),
             hline=0.8, mar=c(5.1, 4.1, 4.1, 8.5),
             ylim=c(0,1), ylab=NULL, xlab=NULL, legend.inset= c(-.3, 0),
             alpha=NULL, cex.lab=1.0,cex.axis=1.0, ...) {
    ## JPS 2/28/21 items to fix up
    ## lty and col better then 1:10
    ## xlab and ylab allow user to pass, must pull off from dots
    usrpar <- par()
    on.exit(par(usrpar))
  
    sensuniq <- unique(x$powerdat$sens)
    specuniq <- unique(x$powerdat$spec)
    ## add Barms to plot

    if(!is.null(alpha)) {
      x$alpha <- alpha
      cvalue <- round(qnorm(1-alpha/2), 2)  
      x$powerdat$power <- 1 - pnorm(cvalue - x$powerdat$z_Barm) +
          pnorm(-1*cvalue - x$powerdat$z_Barm)
    }
   
    
    ## could have metric of power vs z_Barm
    ## ask Jessie if this is wanted
    if(any(casefold(metric) %in% c("power","pow"))) {    
    ## put horizontal line at 80% power or allow user to set where that line is with "hline"
        ## let ylim override power range so we can see the full range
        ## only if ylim not valid, set to range of power
      if(is.null(ylim) | min(c(ylim,0),na.rm=TRUE) < 0 | max(c(1,ylim),na.rm=TRUE) > 1.0) {
          ylim <- range(x$powerdat$power, na.rm=TRUE)
      }
      if(xtype == "spec") {
        ##power specificity
        xline <- specuniq
        xlab <- ifelse(!is.null(xlab), xlab, expression(paste("Specificity",(pi[2]))))
        ylab <- ifelse(!is.null(ylab), ylab,paste0("Power: B", x$Barm[1]," vs. B",x$Barm[2]))
        par(mar=mar, xpd=TRUE) #c(5.1, 4.5, 4.1, 8.1), xpd=FALSE)
        interaction.plot(x$powerdat$spec, x$powerdat$sens, x$powerdat$power,
                         xtick=FALSE,  col=c(1:10), lwd=2, legend=FALSE,
                         xlab=xlab, ylab=ylab, ylim=ylim, 
                         cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.lab*1.1, ...) 
        legend("topright", inset=legend.inset, # c(-0.42,0.1),
               legend = sensuniq, 
               lty=c(length(sensuniq):1), col=1:length(sensuniq),
               lwd=2,title=expression(paste("Sensitivity",(pi[1]))),bty="n",cex=cex.axis)
      } else if(xtype == "sens") {
        xline <- sensuniq
        par(mar=mar, xpd=TRUE)

        xlab <- ifelse(!is.null(xlab), xlab,expression(paste("Sensitivity",(pi[1]))))
        ylab <- ifelse(!is.null(ylab), ylab,paste0("Power: B", x$Barm[1]," vs. B",x$Barm[2]))
        interaction.plot(x$powerdat$sens, x$powerdat$spec, x$powerdat$power,
                       xtick=FALSE, col=c(1:10), legend=FALSE,
                       xlab=xlab, ylab=ylab, ylim=ylim, 
                       cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.lab*1.1, ...) 
        #par(mar=c(5.1, 5.1, 4.1, 8.1), xpd=TRUE)
        legend("topright", inset=legend.inset, legend = specuniq, 
               lty=c(length(specuniq):1), col=1:length(specuniq),
               lwd=2,title=expression(paste("Specificity",(pi[2]))),bty="n",cex=cex.axis)

      }
      ## ABE, please work on this line not spanning full x range (xline)
      xline <- c(xline, 6) ## this line should not be needed, but JPS hard-coded
      lines(x=c(min(xline), max(xline)), y=c(hline, hline), lty=2, col="gray")  
    }
    if(any(casefold(metric[1]) %in% c("zstat", "zratio", "z_stat", "z_barm"))) {
      ylim=range(x$powerdat$`z_Barm`, na.rm=TRUE)
      if(xtype == "spec") {
        ##power specificity
        xline <- specuniq
        xlab <- ifelse(!is.null(xlab), xlab, expression(paste("Specificity",(pi[2]))))
        ylab <- ifelse(!is.null(ylab), ylab,paste0("Z-stat: B", x$Barm[1]," vs. B",x$Barm[2]))
        par(mar=mar, xpd=TRUE) #c(5.1, 4.5, 4.1, 8.1), xpd=FALSE)
        interaction.plot(x$powerdat$spec, x$powerdat$sens, x$powerdat$`z_Barm`,
                         xtick=FALSE,  col=c(1:10), lwd=2, legend=FALSE,
                         xlab=xlab, ylab=ylab, ylim=ylim, 
                         cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.lab*1.1, ...) 
        legend("topright", inset=legend.inset, legend = sensuniq, 
               lty=c(length(sensuniq):1), col=1:length(sensuniq),
               lwd=2,title=expression(paste("Sensitivity",(pi[1]))),bty="n",cex=cex.axis)
      } else if(xtype == "sens") {
        xline <- sensuniq
        par(mar=mar, xpd=TRUE)
        xlab <- ifelse(!is.null(xlab), xlab,expression(paste("Sensitivity",(pi[1]))))
        ylab <- ifelse(!is.null(ylab), ylab,paste0("Z-stat: B", x$Barm[1]," vs. B",x$Barm[2]))
        interaction.plot(x$powerdat$sens, x$powerdat$spec, x$powerdat$`z_Barm`,
                       xtick=FALSE, col=c(1:10), legend=FALSE,
                       xlab=xlab, ylab=ylab, ylim=ylim, 
                       cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.lab*1.1, ...) 
        #par(mar=c(5.1, 5.1, 4.1, 8.1), xpd=TRUE)
        legend("topright", inset=legend.inset, legend = specuniq, 
               lty=c(length(specuniq):1), col=1:length(specuniq),
               lwd=2,title=expression(paste("Specificity",(pi[2]))),bty="n",cex=cex.axis)

      }
    }
           
    invisible()
}


