
## plot smart DTR object

plot.smartDTR <- function(x, xtype="sens", metric="mean",
                     ylim=NULL, ylab=NULL, xlab=NULL, lwd=2,
                     relativeBias=FALSE, mar=c(4, 4, 4, 6),
                     legend.inset=c(-.2,0), cex.axis = .7, cex.lab=.9, ...) {
    
    if(!("smartDTR" %in% class(x))) {
        stop("Object not a smartDTR class")
    }
    if(any(mar < 0) | length(mar) < 4) {
        stop("invalid vector for mar")
    } 
    if(xtype != "sens" & xtype != "spec"){
        stop("Please enter valid xtype (spec or sens): ", xtype)
    }
    if(metric != "mean" & metric != "variance") {
        stop("Please enter valid metric (mean, variance): ", metric)
    }
    usrpar <- par()
    on.exit(par(usrpar))
    
    ## ticks for legend, so if xtype is sens, then legend is spec
    if(xtype == "sens") {
      ## AEC: u-sens and u-spec to be used below in legend
      ## rather than nrow(dat)
      uspec <- unique(x$dtrdat$spec)
      seqst <- min(uspec)
      seqend <- max(uspec)
      seqtick <- uspec[2] - uspec[1]
      legend.col <- legend.lty <- 1:length(uspec)
    } else {
      usens <- unique(x$dtrdat$sens)
      seqst <- min(usens)
      seqend <- max(usens)
      seqtick <- usens[2] - usens[1]
      legend.col <- legend.lty <- 1:length(usens)
    }  

    if(metric == "mean" & relativeBias ==FALSE) {
    #Plot mixed mean
      mudat <- x$dtrdat[,c("sens", "spec", "mu_Barm")]
      Blabel <- paste0("B",x$Barm[1], "B", x$Barm[2])        
      ylab <- ifelse(is.null(ylab), paste0("Mean of ", Blabel), ylab)
      if(is.null(ylim)) {
          ylim <- range(mudat[,3], na.rm=TRUE)
      }
      if(xtype == "spec") { 
        xlab <- ifelse(is.null(xlab), expression(paste("Specificity",(pi[2]))), xlab) 
        par(mar=mar, xpd=TRUE)
        interaction.plot(mudat$spec, mudat$sens, mudat$mu_Barm,
           xlab=xlab, ylab=ylab, ylim=ylim,
           col=legend.col, lty=legend.lty, lwd=lwd,xtick=FALSE,
           legend=FALSE,cex.lab=cex.lab, cex.axis=cex.axis, cex.main=1.1*cex.lab, ...)
        legend("topright", inset=legend.inset,
           legend = c(seq(seqst,seqend,by=seqtick)),
           lty=legend.lty,col=legend.col,
           lwd=lwd,title=expression(paste("Sensitivity",(pi[1]))),bty="n",cex=cex.axis) 
      }  ## if spec
      if(xtype == "sens") {
        xlab <- ifelse(is.null(xlab), expression(paste("Sensitivity",(pi[1]))), xlab)  
        par(mar=mar, xpd=TRUE)
        interaction.plot(mudat$sens, mudat$spec, mudat$mu_Barm,
           xlab=xlab, ylab=ylab,ylim=ylim,
           col=legend.col, lty=legend.lty, lwd=lwd, xtick=FALSE,
           legend=FALSE,cex.lab=cex.lab, cex.axis=cex.axis,
           cex.main=1.1*cex.lab, ...)

        legend("topright", inset=legend.inset,
           legend = c(seq(seqst,seqend,by=seqtick)),
           lty=legend.lty,col=legend.col,
           lwd=lwd,title=expression(paste("Specificity",(pi[2]))),bty="n",cex=cex.axis) 
      } ## if sens
    }  #metric mean & !relativeBias

    if(metric == "variance" & relativeBias == FALSE) {
      sigdat <- x$dtrdat[,c("sens", "spec", "sigsq_Barm")]
      Blabel <- paste0("B",x$Barm[1], "B", x$Barm[2])        
      ylab <- ifelse(is.null(ylab), paste0("Variance of ", Blabel), ylab)
      if(is.null(ylim)) {
          ylim <- range(sigdat[,3], na.rm=TRUE)
      }
      if(xtype == "spec") {
        xlab <- ifelse(is.null(xlab), expression(paste("Specificity",(pi[2]))), xlab)
        par(mar=mar, xpd=TRUE)
        interaction.plot(sigdat$spec, sigdat$sens, sigdat$sigsq_Barm,
           xlab=xlab, ylab=ylab, ylim=ylim,
           col=legend.col, lty=legend.lty, lwd=lwd, xtick=FALSE,
           legend=FALSE,cex.lab=cex.lab, cex.axis=cex.axis,
           cex.main=1.1*cex.lab, ...)
        legend("topright", inset=legend.inset,
           legend = c(seq(seqst,seqend,by=seqtick)),
           lty=legend.lty,col=legend.col,
           lwd=lwd,title=expression(paste("Sensitivity",(pi[1]))),bty="n",cex=cex.axis) 
      }  ## if spec
      if(xtype == "sens") {
        xlab <- ifelse(is.null(xlab), expression(paste("Sensitivity",(pi[1]))), xlab)  
        par(mar=mar, xpd=TRUE)
        interaction.plot(sigdat$sens, sigdat$spec, sigdat$sigsq_Barm,
           xlab=xlab,   ylab=ylab, ylim=ylim,
           col=legend.col, lty=legend.lty, lwd=lwd,xtick=FALSE,
           legend=FALSE,cex.lab=cex.lab, cex.axis=cex.axis,
           cex.main=1.1*cex.lab, ...)
           
        legend("topright", inset=legend.inset,
           legend = c(seq(seqst,seqend,by=seqtick)),
           lty=legend.lty, col=legend.col,
           lwd=lwd, title=expression(paste("Specificity",(pi[2]))),
           bty="n",cex=cex.axis) 
      } ## if sens
    }  #metric variance
    
    if(metric == "mean" & relativeBias==TRUE) {
      mudat <- x$dtrdat[,c("sens", "spec", "mu_Barm")]

      mudat$bias_Barm <- 100 * (x$true_mumix - mudat$mu_Barm)/x$true_mumix
      Blabel <- paste0("B",x$Barm[1], "B", x$Barm[2])        
      colnames(mudat) <- gsub("Barm", Blabel, colnames(mudat))
      ylab <- ifelse(is.null(ylab), paste0("Bias (%) of mixed mean ", Blabel), ylab)
      if(is.null(ylim)) {
          ylim <- range(mudat[,ncol(mudat)], na.rm=TRUE)
      }
      if(xtype == "spec") {
        xlab <- ifelse(is.null(xlab), expression(paste("Specificity",(pi[2]))), xlab)  
        par(mar=mar, xpd=TRUE)
        interaction.plot(mudat$spec, mudat$sens,
                         mudat[,grep("bias", colnames(mudat))],
           xlab=xlab, ylab=ylab, ylim=ylim,
           col=legend.col, lty=legend.lty, lwd=lwd, xtick=FALSE,
           legend=FALSE,cex.lab=cex.lab, cex.axis=cex.axis,
           cex.main=1.1*cex.lab, ...)           
        legend("topright", inset=legend.inset,
           legend = c(seq(seqst,seqend,by=seqtick)),
           lty=legend.lty, col=legend.col,
           lwd=lwd, title=expression(paste("Sensitivity",(pi[1]))),
           bty="n",cex=cex.axis) 
      }  ## if spec
      if(xtype == "sens") {
        xlab <- ifelse(is.null(xlab), expression(paste("Sensitivity",(pi[1]))), xlab)  
        par(mar=mar, xpd=TRUE)
        interaction.plot(mudat$sens, mudat$spec, mudat[,grep("bias", colnames(mudat))],
           xlab=xlab, ylab=ylab, ylim=ylim,
           col=legend.col, lty=legend.lty, lwd=lwd, xtick=FALSE,
           legend=FALSE,cex.lab=cex.lab, cex.axis=cex.axis,
           cex.main=1.1*cex.lab, ...)
        legend("topright", inset=legend.inset,
           legend = c(seq(seqst,seqend,by=seqtick)),
           lty=legend.lty, col=legend.col,
           lwd=lwd, title=expression(paste("Specificity",(pi[2]))),bty="n",cex=cex.axis) 
      } ## if sens
    }  #metric bias of mixed mean

    ##Bias of mixed variance
    if(metric == "variance" & relativeBias==TRUE) {
      sigdat <- x$dtrdat[,c("sens", "spec", "sigsq_Barm")]

      sigdat$bias_Barm <- 100*(x$true_sigmix - sigdat$sigsq_Barm)/x$true_sigmix

      Blabel <- paste0("B",x$Barm[1], "B", x$Barm[2])        
      colnames(sigdat) <- gsub("Barm", Blabel, colnames(sigdat))
      ylab <- ifelse(is.null(ylab), paste0("Bias of mixed variance ", Blabel, "(%)"), ylab)
      if(is.null(ylim)) {
          ylim <- range(sigdat[,ncol(sigdat)], na.rm=TRUE)
      }
        
      if(xtype == "spec") {
        xlab <- ifelse(is.null(xlab), expression(paste("Specificity",(pi[2]))), xlab)  
        par(mar=mar, xpd=TRUE)
        interaction.plot(sigdat$spec, sigdat$sens, sigdat[,grep("bias", colnames(sigdat))],
           xlab=xlab, ylab=ylab, ylim=ylim,
           col=legend.col, lty=legend.lty, lwd=lwd, xtick=FALSE,
           legend=FALSE,cex.lab=cex.lab, cex.axis=cex.axis,
           cex.main=1.1*cex.lab, ...)
        legend("topright", inset=legend.inset,
           legend = c(seq(seqst,seqend,by=seqtick)),
           lty=legend.lty, col=legend.col,
           lwd=lwd, title=expression(paste("Sensitivity",(pi[1]))),
           bty="n",cex=cex.axis) 
      }  ## if spec
      if(xtype == "sens") {
        xlab <- ifelse(is.null(xlab), expression(paste("Sensitivity",(pi[1]))), xlab)  
        par(mar=mar, xpd=TRUE)
        interaction.plot(sigdat$sens, sigdat$spec, sigdat[,grep("bias", colnames(sigdat))],
           xlab=xlab, ylab=ylab, ylim=ylim,
           col=legend.col, lty=legend.lty, lwd=lwd, xtick=FALSE,
           legend=FALSE,cex.lab=cex.lab, cex.axis=cex.axis,
           cex.main=1.1*cex.lab, ...)
        legend("topright", inset=legend.inset,
           legend = c(seq(seqst,seqend,by=seqtick)),
           lty=legend.lty, col=legend.col, lwd=lwd,
           title=expression(paste("Specificity",(pi[2]))),bty="n",cex=cex.axis) 
      } ## if sens
    }  #metric bias of mixed variance

    invisible()
}
