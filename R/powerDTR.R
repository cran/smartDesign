
## Power of DTR

powerDTR <- function(dtr1, dtr2, pG_A1 = 0.8, pG_A2 = 0.8, alpha=0.05) {
    
## future argument for other power calculations: power="response"
    
    Barm1 <- dtr1$Barm
    Barm2 <- dtr2$Barm
 
    M1 <- dtr1$sst1$sstdat$n_Barm /
        (dtr1$sst1$sstdat$n_Barm + dtr1$sst2$sstdat$n_Barm)
    M3 <- dtr1$sst2$sstdat$n_Barm /
        (dtr1$sst1$sstdat$n_Barm + dtr1$sst2$sstdat$n_Barm)
    M2 <- dtr2$sst1$sstdat$n_Barm /
        (dtr2$sst1$sstdat$n_Barm + dtr2$sst2$sstdat$n_Barm)
    M4 <- dtr2$sst2$sstdat$n_Barm /
        (dtr2$sst1$sstdat$n_Barm + dtr2$sst2$sstdat$n_Barm)
    Num1 <- M1*dtr1$sst1$sstdat$mu_Barm +
            M3*dtr1$sst2$sstdat$mu_Barm -
            M2*dtr2$sst1$sstdat$mu_Barm -
            M4*dtr2$sst2$sstdat$mu_Barm
 
   ## Num1=M1*mu_B1_dat$mu_B1-M2*mu_B2_dat$mu_B2+M3*mu_B3_dat$mu_B3-M4*mu_B4_dat$mu_B4
    Den1 <- sqrt(M1^2 * dtr1$sst1$sstdat$sigsq_Barm/dtr1$sst1$sstdat$n_Barm +
                M2^2 * dtr2$sst1$sstdat$sigsq_Barm/dtr2$sst1$sstdat$n_Barm +
                M3^2 * dtr1$sst2$sstdat$sigsq_Barm/dtr1$sst2$sstdat$n_Barm +
                M4^2 * dtr2$sst2$sstdat$sigsq_Barm/dtr2$sst2$sstdat$n_Barm)
    zstat <- Num1/Den1
    ##    z_B1B3vsB2B4=Num1/Den1
    pval <- 2*pnorm(-abs(zstat))
    ##p_B1B3vsB2B4<-2*pnorm(-abs(z_B1B3vsB2B4))
    c_value <- round(qnorm(1-alpha/2), 2) #1.96 #c_value is critical value
    powerstat <- 1-pnorm(c_value-zstat) + pnorm(-c_value-zstat)
    
    powerdf <- cbind.data.frame(dtr1$dtrdat[,c("sens","spec","mu_Barm")],
                          dtr2$dtrdat[,c("mu_Barm")], 
                          z_stat=zstat, p_value=pval, power=powerstat)
    colnames(powerdf)[c(3,4)] <- c(paste0("mu_B", paste(dtr1$Barm, collapse="")),
                                   paste0("mu_B", paste(dtr2$Barm, collapse="")))

    out <- list(powerdat=powerdf, Barms=cbind(Barm1, Barm2), alpha=alpha) 
    class(out) <- c("powerDTR", "list")
    return(out)
}
    
print.powerDTR <- function(x, ...) {
## cat("\n B-group: ", paste(x$Barms[,1], collapse=", "), "; alpha = ", x$alpha, "\n\n")   
   print(x$powerdat, ...)
   invisible(x$powerdat)
}

plot.powerDTR <- function(x, xtype="spec", metric="power", legend.inset=c(-.3, 0),
                          mar=c(5.1, 4.1, 4.1, 8.1), ylim=NULL, hline=0.8,
                          cex.axis=1.0, cex.lab=1.0, alpha=NULL, ...) {

    usrpar <- par()
    on.exit(par(usrpar))
   
    sensuniq <- sort(unique(x$powerdat$sens), decreasing=FALSE)
    specuniq <- sort(unique(x$powerdat$spec), decreasing=FALSE)
    if(xtype == "spec") {
        seqend <- max(sensuniq)
        seqst <- min(sensuniq)
        seqtick <- sensuniq[2] - sensuniq[1]
        sequniq <- sensuniq
    }
    if(xtype == "sens") {
        seqend <- max(specuniq)
        seqst <- min(specuniq)
        seqtick <- specuniq[2] - specuniq[1]
        sequniq <- specuniq
    }

    if(!is.null(alpha)) {
      x$alpha <- alpha
      cvalue <- round(qnorm(1-alpha/2), 2)
     # browser()
    #  x$powerdat$power <- 1 - pnorm(cvalue - x$powerdat$z_Barm) +
    #      pnorm(-1*cvalue - x$powerdat$z_Barm)
    #  c_value <- round(qnorm(1-alpha/2), 2) #1.96 #c_value is critical value
      x$powerdat$power <- 1-pnorm(cvalue - x$powerdat$z_stat) +
                                   pnorm(-cvalue - x$powerdat$z_stat)
      
    }
    if(any(casefold(metric) == c("pow","power"))) {

      if(is.null(ylim)) ylim=c(0,1)

      if(xtype == "sens") {
        xline <- specuniq
        par(mar=mar, xpd=TRUE)
        interaction.plot(x$powerdat$sens, x$powerdat$spec, x$powerdat$power,
                 xlab=expression(paste("Sensitivity",(pi[1]))),
                 ylab="Power (B1B3 vs. B2B4)",col=c(1:10), lwd=2,
                 legend=FALSE,cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.axis, ylim=ylim) 
        legend("topright", inset=legend.inset,
          legend = c(seq(seqst,seqend,by=seqtick)),
          lty=c(nrow(x$powerdat):1),col=c(1:nrow(x$powerdat)),
          lwd=2,title=expression(paste("Specificity",(pi[2]))),bty="n",cex=cex.axis)
      }
      if(xtype == "spec") {
        xline <- specuniq
        par(mar=mar, xpd=TRUE)
        interaction.plot(x$powerdat$spec, x$powerdat$sens, x$powerdat$power,
                 xlab=expression(paste("Specificity",(pi[2]))),
                 ylab="Power (B1B3 vs. B2B4)",col=c(1:10), lwd=2,
                 legend=FALSE,cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.lab, ylim=ylim) 
        legend("topright", inset=legend.inset,
          legend = c(seq(seqst,seqend,by=seqtick)),
          lty=c(nrow(x$powerdat):1),col=c(1:nrow(x$powerdat)),
          lwd=2,title=expression(paste("Specificity",(pi[2]))),bty="n",cex=cex.axis)
      }
      xline <- c(xline, 6) ## this line should not be needed, but JPS hard-coded
      lines(x=c(min(xline), max(xline)), y=c(hline, hline), lty=2, col="gray")
    }
    if(any(casefold(metric) %in% c("zstat", "z_stat", "zratio"))) {
      if(xtype == "sens") {
        par(mar=mar, xpd=TRUE)  
        interaction.plot(x$powerdat$sens, x$powerdat$spec, x$powerdat$`z_stat`,
                 xlab=expression(paste("Sensitivity",(pi[1]))),
                 ylab="Z Ratio (B1B3 vs. B2B4)",col=c(1:10), lwd=2,
                 legend=FALSE,cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.axis,ylim=ylim) 
        legend("topright", inset=legend.inset,
          legend = c(seq(seqst,seqend,by=seqtick)),
          lty=c(nrow(x$powerdat):1),col=c(1:nrow(x$powerdat)),
          lwd=2,title=expression(paste("Specificity",(pi[2]))),bty="n",cex=cex.axis)
      }
      if(xtype == "spec") {
        par(mar=mar, xpd=TRUE)  
        interaction.plot(x$powerdat$spec, x$powerdat$sens, x$powerdat$`z_stat`,
                 xlab=expression(paste("Sensitivity",(pi[1]))),
                 ylab="Z Ratio (B1B3 vs. B2B4)",col=c(1:10), lwd=2,
                 legend=FALSE,cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.lab, ylim=ylim) 
        legend("topright", inset=legend.inset,
          legend = c(seq(seqst,seqend,by=seqtick)),
          lty=c(nrow(x$powerdat):1),col=c(1:nrow(x$powerdat)),
          lwd=2,title=expression(paste("Specificity",(pi[2]))),bty="n",cex=cex.axis)
      }
    }
    invisible()
}


