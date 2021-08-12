\name{smartSST}
\alias{smartSST}
\alias{print.smartSST}
\alias{plot.smartSST}
\title{Single Sequential Trial design clinical trial calculations}
\description{
  Single Sequential Trial design clinical trial calculations
}
\usage{
smartSST(mu_Barm=c(G1=30, G0=20), sigsq_Barm=c(G1=100, G0=100),
                     nsubject=500,
                     Barm=1, type="continuous",
                     sens=seq(0.5,1, by=0.1), spec=seq(0.5, 1, by=0.1),
                     pG_A1 = 0.8, pG_A2=0.8, pran_A1 = 0.5, pran_Barm = 0.5)
}
\arguments{
  \item{mu_Barm}{named vector of the means for the Barm for the
    smart SST trial, with names 'G1' and 'G0'}
  \item{sigsq_Barm}{named vector of the variances (sigma-squared) for the Barm for the
    smart SST trial, with names 'G1' and 'G0'}
  \item{nsubject}{total sample size for the trial}
  \item{Barm}{for the second phase of the trial, the 'B' level for
    which the means/variances apply}
  \item{type}{trial response variable type; only continuous is
    implemented currently}
  \item{sens}{range of sensitivity for smart SST calculations; (0,1]}
  \item{spec}{range of specificity for smart SST calculations; (0,1]}
  \item{pG_A1}{probability of response to therapy given assignment to A1}
  \item{pG_A2}{probability of response to therapy given assignment to A2}
  \item{pran_A1}{probability of random assignment to A1}
  \item{pran_Barm}{probability of assignment to Barm}
}
\details{
 more details on smart SST
}
\value{
  An object of the smartSST S3 class, with the following elements:
  \item{sstdat: }{data.frame with sens, spec, mu, sigsq and sample size (n)
  }
  \item{mu_Barm: }{The value of mu_Barm passed to the function}
  \item{sigsq_Barm: }{The value of sigsq_Barm passed to the function}
}
\references{
  Jun He, Donna K. McClish & Roy T. Sabo (2021) Evaluating Misclassification Effects on Single Sequential Treatment in Sequential Multiple Assignment Randomized Trial (SMART) Designs, Statistics in Biopharmaceutical Research, DOI: 10.1080/19466315.2021.1883472
}
\author{
Jun (Jessie) He, Aberaham Eyman-Casey, Jason P. Sinnwell, Mayo Clinic
}

\examples{
sst1 <- smartSST(mu_Barm=c(G1=30, G0=20), sigsq_Barm=c(G1=16,G0=16),
        Barm=1, sens=seq(.6, 1, by=.1),  spec=seq(.6, 1, by=.1),
        nsubject=252)
print(sst1$sstdat, digits=2)
}

