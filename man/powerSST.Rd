\name{powerSST}
\alias{powerSST}
\alias{print.powerSST}
\alias{plot.powerSST}
\title{Power for Single Sequential Treatment (SST) Trial design clinical trial calculations}
\description{
Power Calculations Comparing two Single Sequential Treatment Treatment (SST) Trial design clinical trial calculations
}
\usage{
powerSST(sst1, sst2, pG_A1 = 0.8, pG_A2 = 0.8, alpha=0.05)
}
\arguments{
  \item{sst1}{an object of smartSST class, created by function of the
    same name}
  \item{sst2}{an object of smartSST class, created by function of the
    same name}
  \item{pG_A1}{probability of response to therapy given assignment to A1}
  \item{pG_A2}{probability of response to therapy given assignment to A2}
  \item{alpha}{accepted type-I error rate for power calculations}
}
\details{
 more details to come
}
\value{
  An object of the powerSST S3 class, with the following elements:
  \item{powerdat: }{data.frame with sens, spec, mu, sigsq and sample size, power
  }
}
\author{
Jun (Jessie) He, Aberaham Eyman-Casey, Jason P. Sinnwell, Mayo Clinic
}
\examples{
  sst1 <- smartSST(mu_Barm=c(G1=30, G0=20), sigsq_Barm=c(G1=16,G0=16),
         Barm=1, sens=seq(.6, 1, by=.1),  spec=seq(.6, 1, by=.1),
         nsubject=252)
  sst2 <- smartSST(mu_Barm=c(G1=20, G0=30), sigsq_Barm=c(G1=16,G0=16),
         Barm=2, sens=seq(.6, 1, by=.1),  spec=seq(.6, 1, by=.1),
         nsubject=252)

  psst12 <- powerSST(sst1, sst2)
  print(psst12) ## plot(psst12)
  }
