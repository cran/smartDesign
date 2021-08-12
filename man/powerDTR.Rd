\name{powerDTR}
\alias{powerDTR}
\alias{print.powerDTR}
\alias{plot.powerDTR}
\title{Power Dynamic Treatment Regimen (DTR) Trial design clinical trial calculations}
\description{
Power Calculations Comparing two Dynamic Treatment Regimen (DTR) Trial design clinical trial calculations
}
\usage{
powerDTR(dtr1, dtr2, pG_A1 = 0.8, pG_A2 = 0.8, alpha=0.05)
}
\arguments{
  \item{dtr1}{an object of smartDTR class, created by function of the
    same name}
  \item{dtr2}{an object of smartDTR class, created by function of the
    same name}
  \item{pG_A1}{probability of response to therapy given assignment to A1}
  \item{pG_A2}{probability of response to therapy given assignment to A2}
  \item{alpha}{accepted type-I error rate for power calculations}
}
\details{
 more details on power DTR
}
\value{
  An object of the powerDTR S3 class, with the following elements:
  \item{powerdat: }{data.frame with sens, spec, mu, sigsq and sample
    size, power
  }
}
\author{
Jun (Jessie) He, Aberaham Eyman-Casey, Jason P. Sinnwell, Mayo Clinic
}
\examples{
  mumat13 <- cbind(G1=c(30,35), G0=c(20,28))
  varmat13 <- cbind(G1=c(100,100),G0=c(100,100))

  dtr13 <- smartDTR(mu_Barm=mumat13, sigsq_Barm=varmat13,
                   Barm=c(1,3), nsubject=252, pG_A1=0.8)

  mumat24 <- cbind(G1=c(25,32), G0=c(18,23))
  varmat24 <- cbind(G1=c(100,100),G0=c(100,100))

  dtr24 <- smartDTR(mu_Barm=mumat24, sigsq_Barm=varmat24,
                   Barm=c(2,4), nsubject=252, pG_A1=0.8, pG_A2=0.8)

  pdtr13vs24 <- powerDTR(dtr13, dtr24)
  print(pdtr13vs24)  ## plot(pdtr13vs24)
}

