

smartSST <- function(mu_Barm=c(G1=30, G0=20), sigsq_Barm=c(G1=100, G0=100),
                     nsubject=500,
                     Barm=1, type="continuous",
                     sens=seq(0.5,1, by=0.1), spec=seq(0.5, 1, by=0.1),
                     pG_A1 = 0.8, pG_A2=0.8, pran_A1 = 0.5, pran_Barm = 0.5) {
    
    ## validate input
    ## pG < 1
    if(any(c(pG_A1, pG_A2) >= 1) | any(c(pG_A1,pG_A2) >= 1)) {
        stop("pG values outside of range (0,1)\n")
    }
    ## Barm in 1:8
    if(!(Barm %in% c(1:8))) {
        stop("Please choose Barm between 1 and 8: ", Barm)
    }
    ## names of mu and sigma just have G1 and G0
    if(!all(c("G1","G0") %in% names(mu_Barm))) {
        stop("mu vector should have named elements 'G1', 'G0'\n")
    }
    if(!all(c("G1","G0") %in% names(sigsq_Barm))) {
        stop("sigsq vector should have named elements 'G1', 'G0'\n")
    }
    ## sens/spec in (0,1)
    if(length(spec) <2 | length(sens) < 2) {
       stop("Specifity and Sensitivity need to have length > 1.")
    }
    if(any(spec > 1 | spec <= 0)) {
       stop("Specificity vector contains values <=0 or >1")
    }
    if(any(sens > 1 | sens <= 0)) {
       stop("Sensitivity vector contains values <=0 or >1")
    }
    if(pran_Barm > 1 | pran_Barm <=0) {
        stop("pran_Barm should be in (0,1)\n")
    }
     if(pran_A1 > 1 | pran_A1 <=0) {
        stop("pran_A1 should be in (0,1)\n")
    }
    ## for SST and Barm,
    ## only need two means, G0 and G1, which are correct and mis-classified
    ## two variances, for G0 and G1
    ## N_Barm is calculated from N, and pran_A1 and pran_Barm
    ## responder to A1 for that b-level.
    if(length(mu_Barm) != 2 | length(sigsq_Barm) != 2 ) {
        stop("length of mu, sigsq, and n for Barm should be 2")
    }
    
    ## nB1 is N * pG_A1 * pran_B1
    ## 
    sens_mat <- as.matrix(sens, ncol=1) #pie1 is sensitivity in treatment A1
    spec_mat <- as.matrix(spec, ncol=1) #pie2 is specificity in treatment A1    
  
    if(Barm %in% 1:4) {
       ## note: when sens_mat=spec_mat=1, no misclassification;
       n_Bmat <- sigsq_B_matrix <- mu_B_matrix <- w1_A1_matrix <- w2_A1_matrix <-
         pR_A1B34_matrix <- pR_A1B12_matrix <- matrix(0,nrow=nrow(sens_mat), ncol=nrow(spec_mat))
       n_A1 <- nsubject * pran_A1
       
       for(i in 1:nrow(sens_mat)){
           for(j in 1:nrow(spec_mat)){
               ##w1 is wt related to responders, not going to show in the plot
               w1_A1_matrix[i,j] <- sens_mat[i,]*pG_A1/(sens_mat[i,]*pG_A1 + 
                                    (1-spec_mat[j,])*(1-pG_A1))    
          
               ##w2 is the weight related to non-responders
               w2_A1_matrix[i,j] <- spec_mat[j,]*(1-pG_A1)/(spec_mat[j,]*(1-pG_A1) +
                                           (1-sens_mat[i,])*pG_A1)

               pR_A1B12_matrix[i,j] <- sens_mat[i,]*pG_A1+(1-spec_mat[j,])*(1-pG_A1)
               pR_A1B34_matrix[i,j] <- (1-sens_mat[i,])*pG_A1+(spec_mat[j,])*(1-pG_A1)
               ##mean B1
               if(Barm == 1) {
                   n_Bmat[i,j] <- n_A1 * pR_A1B12_matrix[i,j] * pran_Barm
                  # n_B1 <- n_A1 * pR_A1_dat$pR_A1 * Pran_B1B2     # the number of responders who are randomly assigned into treatment B1
                   mu_B_matrix[i,j] <- w1_A1_matrix[i,j] *
                       mu_Barm['G1'] + (1-w1_A1_matrix[i,j]) * mu_Barm['G0']
                   sigsq_B_matrix[i,j] <- w1_A1_matrix[i,j] *
                         (sigsq_Barm["G1"] + mu_Barm["G1"]^2) + (1-w1_A1_matrix[i,j])*
                         (sigsq_Barm["G0"] + mu_Barm["G0"]^2) - mu_B_matrix[i,j]^2 
               }
               ## mean B2
               if(Barm == 2) {
                   n_Bmat[i,j] <- n_A1 * pR_A1B12_matrix[i,j] * (1-pran_Barm)
                   mu_B_matrix[i,j] <- w1_A1_matrix[i,j]*
                       mu_Barm['G1']+(1-w1_A1_matrix[i,j])*mu_Barm['G0']
                    sigsq_B_matrix[i,j] <- w1_A1_matrix[i,j] *
                         (sigsq_Barm["G1"] + mu_Barm["G1"]^2) + (1-w1_A1_matrix[i,j])*
                         (sigsq_Barm["G0"] + mu_Barm["G0"]^2) - mu_B_matrix[i,j]^2
               }
               ## mean B3
               if(Barm == 3) {
                   n_Bmat[i,j] <- n_A1*pR_A1B34_matrix[i,j]*pran_Barm
                   mu_B_matrix[i,j] <- w2_A1_matrix[i,j] * mu_Barm['G0'] +
                       (1-w2_A1_matrix[i,j]) * mu_Barm['G1']                 
                   sigsq_B_matrix[i,j]<- w2_A1_matrix[i,j] *
                       (sigsq_Barm["G0"] + mu_Barm["G0"]^2) +
                       (1-w2_A1_matrix[i,j])*(sigsq_Barm["G1"]+mu_Barm["G1"]^2) -
                       mu_B_matrix[i,j]^2 
               }
               ## mean B4
               if(Barm == 4) {
                   n_Bmat[i,j] <- n_A1*pR_A1B34_matrix[i,j] * (1-pran_Barm)
                  ## n_Bmat[i,j] <- n_A1*pR_A1B12_matrix[i,j] * (1-pran_Barm)
                   mu_B_matrix[i,j]<-w2_A1_matrix[i,j]*
                       mu_Barm['G0']+(1-w2_A1_matrix[i,j])*mu_Barm['G1']
                   sigsq_B_matrix[i,j] <- w2_A1_matrix[i,j] *
                       (sigsq_Barm["G0"] + mu_Barm["G0"]^2) +
                       (1-w2_A1_matrix[i,j]) * (sigsq_Barm["G1"] + mu_Barm["G1"]^2) -
                          mu_B_matrix[i,j]^2
               }                              
           } ## for(j)
       }  ## for(i)
        
       sens.dat <- data.frame(sens=as.vector(t(matrix(sens_mat,nrow(sens_mat), nrow(sens_mat))))) 
       spec.dat <- data.frame(spec=rep(spec_mat,nrow(spec_mat)))   
    }  # Barm in 1:4
    
    if(Barm %in% c(5:8)) {
       n_Bmat <- sigsq_B_matrix <- mu_B_matrix <- w1_A2_matrix <- w2_A2_matrix <-
         pR_A2B56_matrix <- pR_A2B78_matrix <- matrix(0,nrow=nrow(sens_mat), ncol=nrow(spec_mat))
       n_A2 <- nsubject * (1-pran_A1) 
       for(i in 1:nrow(sens_mat)){
           for(j in 1:nrow(spec_mat)){
               ## w1 is positive predict value
               w1_A2_matrix[i,j] <- sens_mat[i,]*pG_A2/(sens_mat[i,]*
                                                      pG_A2+(1-spec_mat[j,])*(1-pG_A2))
               ## w2 is negative predict value
               w2_A2_matrix[i,j] <- spec_mat[j,]*(1-pG_A2)/(spec_mat[j,]*
                                                         (1-pG_A2)+(1-sens_mat[i,])*pG_A2)
               
               ## pR_Amatricies, used in sample size calculations
               pR_A2B56_matrix[i,j] <- sens_mat[i,]*pG_A2+(1-spec_mat[j,])*(1-pG_A2)
               pR_A2B78_matrix[i,j] <- (1-sens_mat[i,])*pG_A2+(spec_mat[j,])*(1-pG_A2)
               ## mean B5
               if(Barm == 5) {
                   n_Bmat[i,j] <- n_A2*pR_A2B56_matrix[i,j] * pran_Barm
                   mu_B_matrix[i,j]<-w1_A2_matrix[i,j]*
                       mu_Barm["G1"]+(1-w1_A2_matrix[i,j])*mu_Barm["G0"]
                   sigsq_B_matrix[i,j]<-w1_A2_matrix[i,j] *
                       (sigsq_Barm["G1"] + mu_Barm["G1"]^2) + (1-w1_A2_matrix[i,j])*
                       (sigsq_Barm["G0"] + mu_Barm["G0"]^2) - mu_B_matrix[i,j]^2
               }
               ## mean B6
               if(Barm == 6) {
                   n_Bmat[i,j] <- n_A2*pR_A2B56_matrix[i,j] * 1-(pran_Barm)
                   mu_B_matrix[i,j]<-w1_A2_matrix[i,j]*
                       mu_Barm["G1"]+(1-w1_A2_matrix[i,j])*mu_Barm["G0"]
                   sigsq_B6_matrix[i,j]<-w1_A2_matrix[i,j] *
                       (sigsq_Barm["G1"] + mu_Barm["G1"]^2)+(1-w1_A2_matrix[i,j])*
                       (sigsq_Barm["G0"] + mu_Barm["G0"]^2)-mu_B_matrix[i,j]^2 
               }
     
               ## mean/variance B7
               if(Barm == 7) {
                   n_Bmat[i,j] <- n_A2*pR_A2B78_matrix[i,j] * (pran_Barm)
                   mu_B_matrix[i,j] <- w2_A2_matrix[i,j]*
                       mu_Barm["G0"] + (1-w2_A2_matrix[i,j]) * mu_Barm["G1"]
                   sigsq_B_matrix[i,j]<-w2_A2_matrix[i,j] *
                       (sigsq_Barm["G0"] + mu_Barm["G0"]^2) + (1-w2_A2_matrix[i,j])*
                       (sigsq_Barm["G1"] + mu_Barm["G1"]^2) - mu_B_matrix[i,j]^2  
               }
               ## mean/var B8
               if(Barm == 8) {
                   n_Bmat[i,j] <- n_A2*pR_A2B78_matrix[i,j] * 1-(pran_Barm)
                   mu_B_matrix[i,j]<-w2_A2_matrix[i,j]*
                       mu_Barm["G0"]+(1-w2_A2_matrix[i,j]) * mu_Barm["G1"]
                   sigsq_B8_matrix[i,j]<-w2_A2_matrix[i,j] *
                       (sigsq_Barm["G0"] + mu_Barm["G0"]^2) + (1-w2_A2_matrix[i,j])*
                       (sigsq_Barm["G1"] + mu_Barm["G1"]^2) - mu_B_matrix[i,j]^2  
               }         
       }  ## for(j)
     }  ## for(i)
#     pR_A_matrix <- pR_A2_matrix
#     n_Barm <- (nsubject * (1-pran_A1)) * as.vector(t(pR_A_matrix)) * pran_Barm  
     sens.dat <- data.frame(sens=as.vector(t(matrix(sens_mat,nrow(sens_mat),nrow(sens_mat)))))
     spec.dat <- data.frame(spec=rep(spec_mat,nrow(spec_mat)))
       
   }  #if: Barm 5:8                   

   ## JPS changing on 12/13, to make mu_B_dat and sigsq_B_dat,
   ## which is how the plot and DTR needs them anyway
   sstdat <- data.frame(sens=sens.dat, spec=spec.dat,
                        mu_Barm=as.vector(t(mu_B_matrix)),
                        sigsq_Barm=as.vector(t(sigsq_B_matrix)),
                        n_Barm=as.vector(t(n_Bmat)))
                          
   out <- list(sstdat=sstdat,
             mu_Barm = mu_Barm, sigsq_Barm=sigsq_Barm,
             nsubject=nsubject, Barm = Barm,
             # n_Barm = n_Barm,
             #mu_B_matrix = mu_B_matrix, sigsq_B_matrix=sigsq_B_matrix,
             #sens.dat = sens.dat, spec.dat = spec.dat,
             #pR_A_matrix = pR_A_matrix,
             pG_A1 = pG_A1, pG_A2=pG_A2, pran_A1=pran_A1,  pran_Barm=pran_Barm)
  
    ## to assign class, for S3 method dispatch (plot, print)    
   class(out) <- c("smartSST", "list")
   return(out)
}

print.smartSST <- function(x, ...) {
  cat("smart SST object with B-level:", x$Barm, "\n")
  print(x$sstdat, ...)
  invisible(x$sstdat)
}
