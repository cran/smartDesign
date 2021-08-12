#' plot.smartSST
#'
#' Plot for \code{smartSST} object.
#'
#' @param x an object of class \code{\link{meanSST}}
#' @param Barm explain
#' @param xtype explain
#' @param metric explain
#' @param ylim explain
#' @param mar explain here
#' @param legend.inset explain here
#' @return Nothing is returned 
#' @seealso \code{\link{meanSST}}
#' @examples
#' @author AEC
#' @name plotSmart
NULL
#> NULL

#' @rdname plot.meanSST
#' @export
#'

plot.smartSST <- function(x, xtype = "spec", metric="mean", ylim=NULL, ylab=NULL, xlab=NULL,
                      relativeBias=FALSE, mar=c(4, 4, 4, 6), lwd=2,
                      legend.inset=c(-.2,0), cex.axis = .7, cex.lab=.9, ...) {
  if(!("smartSST" %in% class(x))) {
      stop("Object not a smartSST class")
  }
    
  ## verify inputs
  if(any(mar < 0) | length(mar) < 4) {
      stop("invalid vector for mar")
  }  
  if(xtype != "sens" & xtype != "spec"){
    stop("Please enter valid xtype (spec or sens): ", xtype)
  }
  if(metric != "mean" & metric != "variance") {
    stop("Please enter valid metric (mean, variance): ", metric)
  }
  ## ticks for legend, so if xtype is sens, then legend is spec
  if(xtype == "sens") {
      ## AEC: u-sens and u-spec to be used below in legend
      ## rather than nrow(dat)
      uspec <- unique(x$sstdat$spec)
      seqst <- min(uspec)
      seqend <- max(uspec)
      seqtick <- uspec[2] - uspec[1]
      legend.col <- legend.lty <- 1:length(uspec)
  } else {
      usens <- unique(x$sstdat$sens)
      seqst <- min(usens)
      seqend <- max(usens)
      seqtick <- usens[2] - usens[1]
      legend.col <- legend.lty <- 1:length(usens)
  }
  
  usrpar <- par()
  on.exit(par(usrpar))
  
  if(metric == "mean" & relativeBias == FALSE) {  
    title = "Expected mean of "
    BLabel = paste0("B", as.character(x$Barm))
    metricname <- "mu"
    
    dat <- x$sstdat[,c("sens","spec","mu_Barm")]
    colnames(dat)[ncol(dat)] <- paste0("mu_B", x$Barm)
    if(is.null(ylim)) {
       ylim <- range(dat[,3], na.rm=TRUE)
    }
    ylab <- ifelse(is.null(ylab), paste0(title, BLabel), ylab)
    par(mar=mar, xpd=TRUE)
    if(xtype == "spec"){
        xlab <- ifelse(is.null(xlab), expression(paste("Specificity",(pi[2]))), xlab)      
        interaction.plot(dat$spec, dat$sens, dat[,3],
                 xlab=xlab, xtick=TRUE, 
                 ylab=ylab, col=c(1:nrow(dat)),
                 lwd=lwd, ylim=ylim, 
                 legend=FALSE, cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.lab, ...)
        
        legend("topright", inset = legend.inset, 
               legend = c(seq(seqst, seqend, by= seqtick)), 
               lty=legend.lty, col=legend.col, ncol=1,
               lwd=lwd, title = expression(paste("Sensitivity", (pi[1]))), bty="n", cex=cex.axis)
        ## if(spec) 
    } else if(xtype == "sens") {
        xlab <- ifelse(is.null(xlab), expression(paste("Sensitivity",(pi[1]))), xlab)
        interaction.plot(dat$sens, dat$spec, dat[,3],
                 xlab=xlab, xtick=TRUE, ylim=ylim,
                 ylab=ylab, col=legend.col, lty=legend.lty, lwd=lwd,
                 legend=FALSE, cex.lab=cex.lab, cex.axis=cex.axis, cex.main=cex.lab, ...)
        
        legend("topright", inset=legend.inset, 
               legend = c(seq(seqst, seqend, by=seqtick)), 
               lty=legend.lty, col=legend.col, ncol=1,
               lwd=lwd, title = expression(paste("Specificity", (pi[2]))), bty="n", cex=cex.axis)
        
    } ## if(sens)

  } ## metric == "mean" & !relativeBias

  if(metric == "variance" & relativeBias == FALSE) {

      title = "Variance of "
      BLabel = paste0("B", as.character(x$Barm))
      dat <- x$sstdat[,c("sens","spec","sigsq_Barm")]
      colnames(dat)[ncol(dat)] <- paste0("sigsq_B", x$Barm)
      if(is.null(ylim)) {
          ylim <- range(dat[,3], na.rm=TRUE)
      }
      ylab <- ifelse(is.null(ylab), paste0(title, BLabel), ylab)
      par(mar=mar, xpd=TRUE)
      if(xtype == "spec") {
          xlab <- ifelse(is.null(xlab), expression(paste("Specificity",(pi[2]))), xlab)
          interaction.plot(dat$spec, dat$sens, dat[,3],
             xlab=xlab, xtick = TRUE,
             ylab=ylab, ylim=ylim, col=legend.col, lwd=lwd, lty=legend.lty,
             legend=FALSE,cex.lab=cex.lab,cex.axis=cex.axis, cex.main=cex.lab, ...)
          legend("topright", inset=legend.inset, 
                 legend = c(seq(seqst,seqend,by=seqtick)),
                 lty=legend.lty,col=legend.col,ncol=1,
                 lwd=lwd,title=expression(paste("Sensitivity",(pi[1]))),bty="n",cex=cex.axis)
      } else if(xtype == "sens") {
          xlab <- ifelse(is.null(xlab), expression(paste("Sensitivity",(pi[1]))), xlab)
          interaction.plot(dat$sens, dat$spec, dat[,3],
               xlab=xlab, xtick = TRUE,
               ylab=ylab, ylim=ylim, col=legend.col,lwd=lwd, lty=legend.lty,
               legend=FALSE,cex.lab=cex.lab,cex.axis=cex.axis, cex.main=cex.lab, ...)
          legend("topright", inset=legend.inset, 
               legend = c(seq(seqst,seqend,by=seqtick)),
               lty=legend.lty,col=legend.col, ncol=1,lwd=lwd,
             title=expression(paste("Specificity",(pi[2]))),bty="n",cex=cex.axis)
    }
  } ## variance metric

  if(metric == "mean" & relativeBias == TRUE) {
    title = "Relative bias of mean  "
    BLabel = paste0("B", as.character(x$Barm), " (%)")
  
    dat <- x$sstdat[,c("sens","spec","mu_Barm", "sigsq_Barm")]
    mu_true <- ifelse(x$Barm %in% c(1,2,5,6), x$mu_Barm["G1"], x$mu_Barm["G0"])
    #mu_spec1 <- dat$mu_Barm[dat$spec==1 & dat$sens == 1]
    dat$bias_Barm <- 100 * (mu_true - dat$mu_Barm)/mu_true
    colnames(dat) <- gsub("Barm", BLabel, colnames(dat))
    if(is.null(ylim)) {
        ylim <- range(dat[,ncol(dat)], na.rm=TRUE)
    }
    ylab <- ifelse(is.null(ylab), paste0(title, BLabel), ylab)
    
    par(mar=mar, xpd=TRUE)
    if(xtype == "spec") {
      xlab <- ifelse(is.null(xlab), expression(paste("Specificity",(pi[2]))), xlab)  
      interaction.plot(dat$spec,dat$sens,dat[,grep("bias", colnames(dat))],
                 xlab=xlab, xtick=TRUE, ylab=ylab, ylim=ylim, col=legend.col,
                 lty=legend.lty, lwd=lwd, legend=FALSE,cex.lab=cex.lab,
                 cex.axis=cex.axis, cex.main=1.1*cex.lab, ...)
         ##,trace.label="Sensitivity")
      legend("topright",inset=legend.inset, 
             legend = c(seq(seqst,seqend,by=seqtick)),
             lty=legend.lty, col=legend.col, ncol=1,bty="n",cex=cex.axis,
             lwd=lwd,title=expression(paste("Sensitivity",(pi[1]))))
    }
    if(xtype == "sens") {
        xlab <- ifelse(is.null(xlab), expression(paste("Sensitivity",(pi[1]))), xlab)
        interaction.plot(dat$sens, dat$spec,dat[,grep("bias", colnames(dat))],
                 xlab=xlab, xtick=TRUE, ylim=ylim,
                 ylab=ylab, col=legend.col, lwd=lwd, lty=legend.lty,
                 legend=FALSE,cex.lab=cex.lab, cex.axis=cex.axis, cex.main=1.1*cex.lab)
       
        legend("topright", inset=legend.inset, 
               legend = c(seq(seqst,seqend,by=seqtick)),
               lty=legend.lty,col=legend.col,
        ncol=1,lwd=2,title=expression(paste("Specificity",(pi[2]))),
        bty="n",cex=cex.axis)
    }
  }
  if(metric == "variance" & relativeBias == TRUE) {
      title = "Relative bias of variance "
      BLabel = paste0("B", as.character(x$Barm), " (%)")
      
      dat <- x$sstdat[,c("sens","spec","mu_Barm", "sigsq_Barm")]
      sigsq_true <- ifelse(x$Barm %in% c(1,2,5,6), x$sigsq_Barm["G1"], x$sigsq_Barm["G0"])
      dat$bias_Barm <- 100 * (sigsq_true - dat$sigsq_Barm)/sigsq_true
      colnames(dat) <- gsub("Barm", BLabel, colnames(dat))
      if(is.null(ylim)) {
          ylim <- range(dat[,ncol(dat)], na.rm=TRUE)
      }
      ylab <- ifelse(is.null(ylab), paste0(title, BLabel), ylab)
      
      par(mar=mar, xpd=TRUE)
      if(xtype == "spec") {
          xlab <- ifelse(is.null(xlab), expression(paste("Specificity",(pi[2]))), xlab)
          interaction.plot(dat$spec,dat$sens,dat[,grep("bias", colnames(dat))],
                           xlab=xlab, ylab=ylab, xtick=TRUE, ylim=ylim,
                           col=legend.col, lwd=lwd, legend=FALSE,
                           cex.lab=cex.lab, cex.axis=cex.axis, cex.main=1.1*cex.lab, ...)

          legend("topright",inset=legend.inset, 
                 legend = c(seq(seqst,seqend,by=seqtick)),
                 lty=legend.lty,col=legend.col,ncol=1,
                 lwd=2,title=expression(paste("Sensitivity",(pi[1]))),bty="n",cex=cex.axis)
      }
      if(xtype == "sens") {
          xlab <- ifelse(is.null(xlab), expression(paste("Sensitivity",(pi[1]))), xlab)
          interaction.plot(dat$sens, dat$spec,dat[,grep("bias", colnames(dat))],
                           xlab=xlab, xtick=TRUE, ylim=ylim,
                           ylab=ylab, col=legend.col, lwd=lwd,
                   legend=FALSE,cex.lab=cex.lab, cex.axis=cex.axis, cex.main=1.1*cex.lab, ...)
          
          legend("topright", inset=legend.inset, 
                 legend = c(seq(seqst,seqend,by=seqtick)),
                 lty=legend.lty,col=legend.col,
                 ncol=1,lwd=lwd, title=expression(paste("Specificity",(pi[2]))),
                 bty="n",cex=cex.axis)
      }
    } ## variance & relativeBias

  invisible()      
}

