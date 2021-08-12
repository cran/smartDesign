
## combine mean and variance calculations
smartDTR <- function(mu_Barm=cbind(G1=c(30,25), G0=c(20,20)),
                     sigsq_Barm=cbind(G1=c(100,100), G0=c(100,100)),
                     nsubject=500,  Barm=c(1,3), type="continuous",  
                     sens=seq(0.5,1, by=0.1), spec=seq(0.5, 1, by=0.1),
                     pG_A1 = 0.8, pG_A2 = 0.8, pran_A1 = 0.5, 
                     pran_Barm = c(0.5, 0.5)) {
    
    ## to-do: validate input
    ## call smartSST for each blevel in the bgroup,
    ## use sstdat from each
    
    sst1 <- smartSST(mu_Barm[1,], sigsq_Barm[1,], nsubject=nsubject, type=type,
                     Barm=Barm[1], sens=sens, spec=spec,pG_A1=pG_A1,
                     pG_A2=pG_A2, pran_A1=pran_A1, pran_Barm=pran_Barm[1])
    sst2 <- smartSST(mu_Barm[2,], sigsq_Barm[2,], nsubject=nsubject, type=type,
                     Barm=Barm[2], sens=sens, spec=spec, pG_A1=pG_A1,
                     pG_A2=pG_A2, pran_A1=pran_A1, pran_Barm=pran_Barm[2])

    nA1 <- nsubject*pran_A1
    if(all(Barm %in%  c(1,3))) {
        n_B1B3 <- sst1$sstdat$n_Barm + sst2$sstdat$n_Barm
        
        n_B1G1 <- nA1 * sst1$sstdat$sens * sst1$pG_A1 * sst1$pran_Barm
        n_B1G0 <- nA1 * (1-sst1$sstdat$spec) * (1-sst1$pG_A1) * sst1$pran_Barm
        n_B3G1 <- nA1 * (1-sst2$sstdat$sens) * sst2$pG_A1 * sst2$pran_Barm
        n_B3G0 <- nA1 * (sst2$sstdat$spec) * (1-sst2$pG_A1) * sst2$pran_Barm
       
        ##mixed mean 
        mu_B1B3 <- (sst1$sstdat$n_Barm * sst1$sstdat$mu_Barm +
                    sst2$sstdat$n_Barm * sst2$sstdat$mu_Barm)/(n_B1B3)
        ## true for rel bias plots
        true_mumix <- (sst1$pG_A1*sst1$pran_Barm*mu_Barm[1,1]+
                      (1-sst1$pG_A1)*sst2$pran_Barm*mu_Barm[2,2])/
               (sst1$pG_A1*sst1$pran_Barm + (1-sst1$pG_A1)*sst2$pran_Barm)
        ##true_mu_B1B3 <- (pG_A1*Pran_B1B2*mu_A1B1G1R1+(1-pG_A1)*Pran_B3B4*mu_A1B3G0R0)/(pG_A1*Pran_B1B2+(1-pG_A1)*Pran_B3B4)
        ##mixed sigmasquare
        sigsq_B1B3 <- (n_B1G1*(sigsq_Barm[1,1] + mu_Barm[1,1]^2) +
                       n_B1G0*(sigsq_Barm[1,2] + mu_Barm[1,2]^2) +
                       n_B3G1*(sigsq_Barm[2,1] + mu_Barm[2,1]^2) +
                       n_B3G0*(sigsq_Barm[2,2] + mu_Barm[2,2]^2))/n_B1B3 -
                     mu_B1B3^2
  ## from Jessie's code:
     #   sigsq_B1B3=(n_B1G1*(sigsq_A1B1G1R1+mu_A1B1G1R1^2) +
     #       n_B1G0*(sigsq_A1B1G0R1+mu_A1B1G0R1^2) +
     #       n_B3G1*(sigsq_A1B3G1R0+mu_A1B3G1R0^2) + 
     #       n_B3G0*(sigsq_A1B3G0R0+mu_A1B3G0R0^2))/n_B1B3 - mu_B1B3^2  # DTRs variance 
        true_sigmix <- (sst1$pG_A1*(sigsq_Barm[1,1]+mu_Barm[1,1]^2) * sst1$pran_Barm +
                       (1-sst1$pG_A1)*(sigsq_Barm[2,2]+mu_Barm[2,2]^2)*(1-sst2$pran_Barm))/
             (sst1$pG_A1*sst1$pran_Barm+(1-sst1$pG_A1)*sst2$pran_Barm) - (true_mumix)^2 


        dtrdat <- data.frame(sens=sst1$sstdat$sens, spec=sst1$sstdat$spec,
                     n_Barm=n_B1B3, mu_Barm=mu_B1B3,
                     sigsq_Barm=sigsq_B1B3)
    }
    if(all(Barm %in% c(1,4))) {
        n_B1B4 <- sst1$sstdat$n_Barm + sst2$sstdat$n_Barm
        
        n_B1G1 <- nA1 * sst1$sstdat$sens * sst1$pG_A1 * sst1$pran_Barm
        n_B1G0 <- nA1 * (1-sst1$sstdat$spec) * (1-sst1$pG_A1) * sst1$pran_Barm
        n_B4G1 <- nA1 * (1-sst2$sstdat$sens) * sst2$pG_A1 * sst2$pran_Barm
        n_B4G0 <- nA1 * (sst2$sstdat$spec) * (1-sst2$pG_A1) * sst2$pran_Barm

        mu_B1B4 <- (sst1$sstdat$n_Barm * sst1$sstdat$mu_Barm +
                    sst2$sstdat$n_Barm * sst2$sstdat$mu_Barm) / (n_B1B4)
        ## true for rel bias plots
        true_mumix <- (sst1$pG_A1*sst1$pran_Barm*mu_Barm[1,1]+
                      (1-sst1$pG_A1)*(1-sst2$pran_Barm)*mu_Barm[2,2]) /
               (sst1$pG_A1*sst1$pran_Barm + (1-sst1$pG_A1)*(1-sst2$pran_Barm))

        ##mixed sigmasquare
        sigsq_B1B4 <- (n_B1G1*(sigsq_Barm[1,1] + mu_Barm[1,1]^2) +
                       n_B1G0*(sigsq_Barm[2,1] + mu_Barm[2,1]^2) +
                       n_B4G1*(sigsq_Barm[1,2] + mu_Barm[1,2]^2) +
                       n_B4G0*(sigsq_Barm[2,2] + mu_Barm[2,2]^2))/n_B1B4 -
                      mu_B1B4^2
        true_sigmix <- (sst1$pG_A1*(sigsq_Barm[1,1]+mu_Barm[1,1]^2)*sst1$pran_Barm +
                       (1-sst1$pG_A1)*(sigsq_Barm[2,2]+mu_Barm[2,2]^2)*(1-sst2$pran_Barm))/
           (sst1$pG_A1*sst1$pran_Barm+(1-sst1$pG_A1)*(1-sst2$pran_Barm)) - (true_mumix)^2 
        
        dtrdat <- data.frame(sens=sst1$sstdat$sens, spec=sst1$sstdat$spec,
                             n_Barm=n_B1B4, mu_Barm=mu_B1B4,
                             sigsq_Barm=sigsq_B1B4)
    }
    if(all(Barm %in% c(2,4))) {
        
        n_B2B4 <- sst1$sstdat$n_Barm + sst2$sstdat$n_Barm
        
        n_B2G1 <- nA1 * sst1$sstdat$sens * sst1$pG_A1 * (1-sst1$pran_Barm)
        n_B2G0 <- nA1 * (1-sst1$sstdat$spec) * (1-sst1$pG_A1) * (1-sst1$pran_Barm) 
        n_B4G1 <- nA1 * (1-sst2$sstdat$sens) * sst2$pG_A1 * sst2$pran_Barm
        n_B4G0 <- nA1 * (sst2$sstdat$spec) * (1-sst2$pG_A1) * sst2$pran_Barm

        mu_B2B4 <- (sst1$sstdat$n_Barm * sst1$sstdat$mu_Barm +
                    sst2$sstdat$n_Barm * sst2$sstdat$mu_Barm)/(n_B2B4)

        ## true for rel bias plot
        true_mumix <- (sst1$pG_A1*(1-sst1$pran_Barm)*mu_Barm[1,1]+
                      (1-sst1$pG_A1)*(1-sst2$pran_Barm)*mu_Barm[2,2]) /
               (sst1$pG_A1*(1-sst1$pran_Barm) + (1-sst1$pG_A1)*(1-sst2$pran_Barm))

        ##mixed sigmasquare
        sigsq_B2B4 <- (n_B2G1*(sigsq_Barm[1,1] + mu_Barm[1,1]^2) +
                       n_B2G0*(sigsq_Barm[1,2] + mu_Barm[1,2]^2) +
                       n_B4G1*(sigsq_Barm[2,1] + mu_Barm[2,1]^2) +
                       n_B4G0*(sigsq_Barm[2,2] + mu_Barm[2,2]^2))/n_B2B4 -
                      mu_B2B4^2
        true_sigmix <- (sst1$pG_A1*(sigsq_Barm[1,1]+mu_Barm[1,1]^2)*(1-sst1$pran_Barm) +
                       (1-sst1$pG_A1)*(sigsq_Barm[2,2]+mu_Barm[2,2]^2)*(1-sst2$pran_Barm))/
          (sst1$pG_A1*(1-sst1$pran_Barm) + (1-sst1$pG_A1)*(1-sst2$pran_Barm)) - (true_mumix)^2 
        
        dtrdat <- data.frame(sens=sst1$sstdat$sens, spec=sst1$sstdat$spec,
                             n_Barm=n_B2B4, mu_Barm=mu_B2B4,
                             sigsq_Barm=sigsq_B2B4)

    }
    if(all(Barm %in%  c(2,3))) {
        n_B2B3 <- sst1$sstdat$n_Barm + sst2$sstdat$n_Barm
        
        n_B2G1 <- nA1 * sst1$sstdat$sens * sst1$pG_A1 * (1-sst1$pran_Barm)
        n_B2G0 <- nA1 * (1-sst1$sstdat$spec) * (1-sst1$pG_A1) * (1-sst1$pran_Barm) 
        n_B3G1 <- nA1 * (1-sst2$sstdat$sens) * sst2$pG_A1 * sst2$pran_Barm
        n_B3G0 <- nA1 * (sst2$sstdat$spec) * (1-sst2$pG_A1) * sst2$pran_Barm
       
        mu_B2B3 <- (sst1$sstdat$n_Barm * sst1$sstdat$mu_Barm +
                    sst2$sstdat$n_Barm * sst2$sstdat$mu_Barm)/(n_B2B3)
        ## true for rel bias plot
        true_mumix <- (sst1$pG_A1*(1-sst1$pran_Barm)*mu_Barm[1,1]+
                      (1-sst1$pG_A1)*(sst2$pran_Barm)*mu_Barm[2,2]) /
               (sst1$pG_A1*(1-sst1$pran_Barm) + (1-sst1$pG_A1)*(sst2$pran_Barm))
        
        ##mixed sigmasquare
        sigsq_B2B3 <- (n_B2G1*(sigsq_Barm[1,1] + mu_Barm[1,1]^2) +
                       n_B2G0*(sigsq_Barm[1,2] + mu_Barm[1,2]^2) +
                       n_B3G1*(sigsq_Barm[2,1] + mu_Barm[2,2]^2) +
                       n_B3G0*(sigsq_Barm[2,2] + mu_Barm[2,2]^2))/n_B2B3 -
                      mu_B2B3^2
        true_sigmix <- (sst1$pG_A1*(sigsq_Barm[1,1]+mu_Barm[1,1]^2) *(1-sst1$pran_Barm) +
                       (1-sst1$pG_A1)*(sigsq_Barm[2,2]+mu_Barm[2,2]^2)*sst2$pran_Barm)/
             (sst1$pG_A1*(1-sst1$pran_Barm) + (1-sst1$pG_A1)*sst2$pran_Barm)-(true_mumix)^2 
        
        dtrdat <- data.frame(sens=sst1$sstdat$sens, spec=sst1$sstdat$spec,
                             n_Barm=n_B2B3, mu_Barm=mu_B2B3,
                             sigsq_Barm=sigsq_B2B3)

    }
    nA2 <- nsubject - nA1
    if(all(Barm %in% c(5,7))) {

        n_B5B7 <- sst1$sstdat$n_Barm + sst2$sstdat$n_Barm
        
        n_B5G1 <- nA2 * sst1$sstdat$sens * (1-sst1$pG_A1) * sst1$pran_Barm
        n_B5G0 <- nA2 * (1-sst1$sstdat$spec) * (sst1$pG_A1) * sst1$pran_Barm
        n_B7G1 <- nA2 * (1-sst2$sstdat$sens) * (1-sst2$pG_A1) * sst2$pran_Barm
        n_B7G0 <- nA2 * (sst2$sstdat$spec) * (sst2$pG_A1) * sst2$pran_Barm

        ##mixed mean  ## JPS use mu_Barm and n_Barm from smartSST call
        mu_B5B7 <- (sst1$sstdat$n_Barm * sst1$sstdat$mu_Barm +
                    sst2$sstdat$n_Barm * sst2$sstdat$mu_Barm)/(n_B5B7)
        true_mumix <- (sst1$pG_A2*sst1$pran_Barm*mu_Barm[1,1]+
                      (1-sst1$pG_A2)*sst2$pran_Barm*mu_Barm[2,2])/
               (sst1$pG_A1*sst1$pran_Barm + (1-sst1$pG_A2)*sst2$pran_Barm)
        
        ##mixed sigmasquare
        sigsq_B5B7 <- (n_B5G1*(sigsq_Barm[1,1] + mu_Barm[1,1]^2) +
                       n_B5G0*(sigsq_Barm[1,2] + mu_Barm[1,2]^2) +
                       n_B7G1*(sigsq_Barm[2,1] + mu_Barm[2,1]^2) +
                       n_B7G0*(sigsq_Barm[2,2] + mu_Barm[2,2]^2))/n_B5B7 -
                      mu_B5B7^2
        true_sigmix <- (sst1$pG_A2*(sigsq_Barm[1,1]+mu_Barm[1,1]^2) * sst1$pran_Barm +
                       (1-sst1$pG_A2)*(sigsq_Barm[2,2]+mu_Barm[2,2]^2)*(1-sst2$pran_Barm))/
             (sst1$pG_A2*sst1$pran_Barm+(1-sst1$pG_A2)*sst2$pran_Barm)-(true_mumix)^2
        
        dtrdat <- data.frame(sens=sst1$sstdat$sens, spec=sst1$sstdat$spec,
                     n_Barm=n_B5B7, mu_Barm=mu_B5B7,
                     sigsq_Barm=sigsq_B5B7)

    }
    if(all(Barm %in% c(5,8))) {
        n_B5B8 <- sst1$sstdat$n_Barm + sst2$sstdat$n_Barm
        
        n_B5G1 <- nA2 * sst1$sstdat$sens * (1-sst1$pG_A1) * sst1$pran_Barm
        n_B5G0 <- nA2 * (1-sst1$sstdat$spec) * (sst1$pG_A1) * sst1$pran_Barm

        n_B8G1 <- nA2 * (1-sst2$sstdat$sens) * (1-sst2$pG_A1) * (1-sst2$pran_Barm)        
        n_B8G0 <- nA2 * (sst2$sstdat$spec) * (sst2$pG_A1) * (1-sst2$pran_Barm)
        
        ##mixed mean  ## JPS use mu_Barm and n_Barm from smartSST call
        mu_B5B8 <- (sst1$sstdat$n_Barm * sst1$sstdat$mu_Barm +
                    sst2$sstdat$n_Barm * sst2$sstdat$mu_Barm)/(n_B5B8)
                ## true for rel bias plots
        true_mumix <- (sst1$pG_A2*sst1$pran_Barm*mu_Barm[1,1]+
                      (1-sst1$pG_A2)*(1-sst2$pran_Barm)*mu_Barm[2,2]) /
               (sst1$pG_A2*sst1$pran_Barm + (1-sst1$pG_A2)*(1-sst2$pran_Barm))

        ##mixed sigmasquare
        sigsq_B5B8 <- (n_B5G1*(sigsq_Barm[1,1] + mu_Barm[1,1]^2) +
                       n_B5G0*(sigsq_Barm[1,2] + mu_Barm[1,2]^2) +
                       n_B8G1*(sigsq_Barm[2,1] + mu_Barm[2,1]^2) +
                       n_B8G0*(sigsq_Barm[2,2] + mu_Barm[2,2]^2))/n_B5B8 -
                      mu_B5B8^2
        true_sigmix <- (sst1$pG_A2*(sigsq_Barm[1,1]+mu_Barm[1,1]^2)*sst1$pran_Barm +
                       (1-sst1$pG_A2)*(sigsq_Barm[2,2]+mu_Barm[2,2]^2)*(1-sst2$pran_Barm))/
           (sst1$pG_A2*sst1$pran_Barm+(1-sst1$pG_A2)*(1-sst2$pran_Barm)) - (true_mumix)^2 

        dtrdat <- data.frame(sens=sst1$sstdat$sens, spec=sst1$sstdat$spec,
                     n_Barm=n_B5B8, mu_Barm=mu_B5B8,
                     sigsq_Barm=sigsq_B5B8)
    }
    if(all(Barm %in% c(6,8))) {
        n_B6B8 <- sst1$sstdat$n_Barm + sst2$sstdat$n_Barm
        
        n_B6G1 <- nA2 * sst1$sstdat$sens * (1-sst1$pG_A1) * (1-sst1$pran_Barm)
        n_B6G0 <- nA2 * (1-sst1$sstdat$spec) * (sst1$pG_A1) * (1-sst1$pran_Barm)

        n_B8G1 <- nA2 * (1-sst2$sstdat$sens) * (1-sst2$pG_A1) * (1-sst2$pran_Barm)        
        n_B8G0 <- nA2 * (sst2$sstdat$spec) * (sst2$pG_A1) * (1-sst2$pran_Barm)
        
        ##mixed mean  ## JPS use mu_Barm and n_Barm from smartSST call
        mu_B6B8 <- (sst1$sstdat$n_Barm * sst1$sstdat$mu_Barm +
                    sst2$sstdat$n_Barm * sst2$sstdat$mu_Barm)/(n_B6B8)
        ## true for rel bias plot
        true_mumix <- (sst1$pG_A2*(1-sst1$pran_Barm)*mu_Barm[1,1]+
                      (1-sst1$pG_A2)*(1-sst2$pran_Barm)*mu_Barm[2,2]) /
               (sst1$pG_A2*(1-sst1$pran_Barm) + (1-sst1$pG_A2)*(1-sst2$pran_Barm))
        ##mixed sigmasquare
        sigsq_B6B8 <- (n_B6G1*(sigsq_Barm[1,1] + mu_Barm[1,1]^2) +
                       n_B6G0*(sigsq_Barm[1,2] + mu_Barm[1,2]^2) +
                       n_B8G1*(sigsq_Barm[2,1] + mu_Barm[2,1]^2) +
                       n_B8G0*(sigsq_Barm[2,2] + mu_Barm[2,2]^2))/n_B6B8 -
                      mu_B6B8^2
        true_sigmix <- (sst1$pG_A2*(sigsq_Barm[1,1]+mu_Barm[1,1]^2)*(1-sst1$pran_Barm) +
                       (1-sst1$pG_A2)*(sigsq_Barm[2,2]+mu_Barm[2,2]^2)*(1-sst2$pran_Barm))/
          (sst1$pG_A2*(1-sst1$pran_Barm) + (1-sst1$pG_A2)*(1-sst2$pran_Barm)) - (true_mumix)^2 

        dtrdat <- data.frame(sens=sst1$sstdat$sens, spec=sst1$sstdat$spec,
                     n_Barm=n_B6B8, mu_Barm=mu_B6B8,
                     sigsq_Barm=sigsq_B6B8)

    }
     if(all(Barm %in% c(6,8))) {
        n_B6B8 <- sst1$sstdat$n_Barm + sst2$sstdat$n_Barm
        
        n_B6G1 <- nA2 * sst1$sstdat$sens * (1-sst1$pG_A1) * (1-sst1$pran_Barm)
        n_B6G0 <- nA2 * (1-sst1$sstdat$spec) * (sst1$pG_A1) * (1-sst1$pran_Barm)

        n_B8G1 <- nA2 * (1-sst2$sstdat$sens) * (1-sst2$pG_A1) * (1-sst2$pran_Barm)        
        n_B8G0 <- nA2 * (sst2$sstdat$spec) * (sst2$pG_A1) * (1-sst2$pran_Barm)
        
        ##mixed mean  ## JPS use mu_Barm and n_Barm from smartSST call
        mu_B6B8 <- (sst1$sstdat$n_Barm * sst1$sstdat$mu_Barm +
                    sst2$sstdat$n_Barm * sst2$sstdat$mu_Barm)/(n_B6B8)

        ##mixed sigmasquare
        sigsq_B6B8 <- (n_B6G1*(sigsq_Barm[1,1] + mu_Barm[1,1]^2) +
                       n_B6G0*(sigsq_Barm[1,2] + mu_Barm[1,2]^2) +
                       n_B8G1*(sigsq_Barm[2,1] + mu_Barm[2,1]^2) +
                       n_B8G0*(sigsq_Barm[2,2] + mu_Barm[2,2]^2))/n_B6B8 -
                      mu_B6B8^2
        dtrdat <- data.frame(sens=sst1$sstdat$sens, spec=sst1$sstdat$spec,
                     n_Barm=n_B6B8, mu_Barm=mu_B6B8,
                     sigsq_Barm=sigsq_B6B8)
    }
    if(all(Barm %in% c(6,7))) {
        n_B6B7 <- sst1$sstdat$n_Barm + sst2$sstdat$n_Barm
        
        n_B6G1 <- nA2 * sst1$sstdat$sens * (1-sst1$pG_A1) * (1-sst1$pran_Barm)
        n_B6G0 <- nA2 * (1-sst1$sstdat$spec) * (sst1$pG_A1) * (1-sst1$pran_Barm)

        n_B7G1 <- nA2 * (1-sst2$sstdat$sens) * (1-sst2$pG_A1) * (sst2$pran_Barm)        
        n_B7G0 <- nA2 * (sst2$sstdat$spec) * (sst2$pG_A1) * (sst2$pran_Barm)
        
        ##mixed mean  ## JPS use mu_Barm and n_Barm from smartSST call
        mu_B6B7 <- (sst1$sstdat$n_Barm * sst1$sstdat$mu_Barm +
                    sst2$sstdat$n_Barm * sst2$sstdat$mu_Barm)/(n_B6B7)
        ## true for rel bias plot
        true_mumix <- (sst1$pG_A2*(1-sst1$pran_Barm)*mu_Barm[1,1]+
                      (1-sst1$pG_A2)*(sst2$pran_Barm)*mu_Barm[2,2]) /
               (sst1$pG_A2*(1-sst1$pran_Barm) + (1-sst1$pG_A2)*(sst2$pran_Barm))

        ##mixed sigmasquare
        sigsq_B6B7 <- (n_B6G1*(sigsq_Barm[1,1] + mu_Barm[1,1]^2) +
                       n_B6G0*(sigsq_Barm[1,2] + mu_Barm[1,2]^2) +
                       n_B7G1*(sigsq_Barm[2,1] + mu_Barm[2,1]^2) +
                       n_B7G0*(sigsq_Barm[2,2] + mu_Barm[2,2]^2))/n_B6B7 -
                      mu_B6B7^2
        true_sigmix <- (sst1$pG_A2*(sigsq_Barm[1,1]+mu_Barm[1,1]^2) *(1-sst1$pran_Barm) +
                       (1-sst1$pG_A2)*(sigsq_Barm[2,2]+mu_Barm[2,2]^2)*sst2$pran_Barm)/
             (sst1$pG_A2*(1-sst1$pran_Barm) + (1-sst1$pG_A2)*sst2$pran_Barm)-(true_mumix)^2 

        dtrdat <- data.frame(sens=sst1$sstdat$sens, spec=sst1$sstdat$spec,
                     n_Barm=n_B6B7, mu_Barm=mu_B6B7,
                     sigsq_Barm=sigsq_B6B7)
    }
    
    out <- list(dtrdat=dtrdat, sst1=sst1, sst2=sst2,
                true_mumix=true_mumix, true_sigmix=true_sigmix, 
                mu_Barm=mu_Barm, sigsq_Barm=sigsq_Barm,Barm = Barm)
#                pG_A1 = pG_A1, pG_A2=pG_A2, pran_A1=pran_A1,  pran_Barm=pran_Barm)   
    ## pG_A1, pran, etc are all within sst objects
    
    ## to assign class, for S3 method dispatch (plot, print)    
    class(out) <- c("smartDTR", "list")
    return(out)
}

print.smartDTR <- function(x, ...) {
  cat("smart DTR object with B-level:", x$Barm, "\n")
  print(x$dtrdat, ...)
  invisible()
}

