\name{smartDTR}
\alias{smartDTR}
\alias{print.smartDTR}
\alias{plot.smartDTR}
\title{Dynamic Treatment Regimen (DTR) Trial design clinical trial calculations}
\description{
Dynamic Treatment Regimen (DTR) Trial design clinical trial calculations
}
\usage{
smartDTR(mu_Barm=cbind(G1=c(30,25), G0=c(20,20)),
                     sigsq_Barm=cbind(G1=c(100,100), G0=c(100,100)),
                     nsubject=500,  Barm=c(1,3), type="continuous",  
                     sens=seq(0.5,1, by=0.1), spec=seq(0.5, 1, by=0.1),
                     pG_A1 = 0.8, pG_A2 = 0.8, pran_A1 = 0.5, 
                     pran_Barm = c(0.5, 0.5))
}
\arguments{
  \item{mu_Barm}{matrix of two named vectors of the means for the two
    B arms (columns) for the smart DTR trial, with rows as 'G1' and 'G0'}
  \item{sigsq_Barm}{matrix of two named vectors of the variances (sigma-squared) for the two
    Blevels (columns) for the smart DTR trial, with rows as 'G1' and 'G0'}
  \item{nsubject}{total sample size for the trial}
  \item{Barm}{for the second phase of the trial, the 'B' levels for
    which the DTR means/variances apply}
  \item{type}{trial response variable type; only continuous is
    implemented currently}
  \item{sens}{range of sensitivity for smart SST calculations; (0,1]}
  \item{spec}{range of specificity for smart SST calculations; (0,1]}
  \item{pG_A1}{probability of response to therapy given assignment to A1}
  \item{pG_A2}{probability of response to therapy given assignment to A2}
  \item{pran_A1}{probability of random assignment to A1}
  \item{pran_Barm}{probability of assignment to Barms}
}
\details{
 see details in the reference
}
\references{
  Jun He, Donna K. McClish & Roy T. Sabo (2021) Evaluating Misclassification Effects on Single Sequential Treatment in Sequential Multiple Assignment Randomized Trial (SMART) Designs, Statistics in Biopharmaceutical Research, DOI: 10.1080/19466315.2021.1883472
}
\value{
  An object of the smartDTR S3 class, with the following elements:
  \item{dtrdat: }{data.frame with sens, spec, mu, sigsq and sample size (n)}
  \item{sst1: }{ smartSST object from the first Barm }
  \item{sst2: }{ smartSST object from the second Barm }
  \item{true_mumix: }{ true mu mixture }
  \item{true_sigmix: }{ true sigma mixture }
  \item{mu_Barm, sigsq_Barm, Barm: }{ input B-arm, mu, and sigsq for DTR }
}
\author{
Jun (Jessie) He, Aberaham Eyman-Casey, Jason P. Sinnwell, Mayo Clinic
}

\examples{
mumat13 <- cbind(G1=c(30,35), G0=c(20,28))
varmat13 <- cbind(G1=c(100,100),G0=c(100,100))

dtr13 <- smartDTR(mu_Barm=mumat13, sigsq_Barm=varmat13,
                 Barm=c(1,3), nsubject=252, pG_A1=0.8)

print(dtr13)
}

