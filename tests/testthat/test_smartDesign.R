## Tests for smartDesign


context("Testing smartSST")


mumat13 <- cbind(G1=c(30,35), G0=c(20,28))
varmat13 <- cbind(G1=c(100,100),G0=c(100,100))
sst1new <- smartSST(mu_Barm=c(G1=30, G0=20), sigsq_Barm=c(G1=16,G0=16),
         Barm=1, sens=seq(.6, 1, by=.1),  spec=seq(.6, 1, by=.1),
         nsubject=252)
sst2new <- smartSST(mu_Barm=c(G1=20, G0=30), sigsq_Barm=c(G1=16,G0=16),
         Barm=2, sens=seq(.6, 1, by=.1),  spec=seq(.6, 1, by=.1),
         nsubject=252)

psst12new <- powerSST(sst1new, sst2new)
#psst12$powerdat
load("SSTsave.RData")
#sst1$sstdat[1:10,]
#sst1new$sstdat[1:10,]
#save(sst1, sst2, psst12, file="SSTsave.RData")
context("Testing smartSST")

test_that("smartSST", {
  expect_equal(sst1new$sstdat, expected=sst1$sstdat)
  expect_equal(sst2new$sstdat, expected=sst2$sstdat)
  expect_equal(psst12new$powerdat, expected=psst12$powerdat)
  })

context("testing smartDTR")
mumat13 <- cbind(G1=c(30,35), G0=c(20,28))
varmat13 <- cbind(G1=c(100,100),G0=c(100,100))

dtr13new <- smartDTR(mu_Barm=mumat13, sigsq_Barm=varmat13,
                 Barm=c(1,3), nsubject=252, pG_A1=0.8)

mumat24 <- cbind(G1=c(25,32), G0=c(18,23))
varmat24 <- cbind(G1=c(100,100),G0=c(100,100))

dtr24new <- smartDTR(mu_Barm=mumat24, sigsq_Barm=varmat24,
                 Barm=c(2,4), nsubject=252, pG_A1=0.8, pG_A2=0.8)


pdtr13vs24new <- powerDTR(dtr13new, dtr24new)
#save(dtr13, dtr24, pdtr13vs24, file="DTRsave.RData")
load("DTRsave.RData")

test_that("smartDTR", {
  expect_equal(dtr13new$dtrdat, expected=dtr13$dtrdat)
  expect_equal(dtr24new$dtrdat, expected=dtr24$dtrdat)
  expect_equal(pdtr13vs24new$powerdat, expected=pdtr13vs24$powerdat)
  })
