\name{MCMCtGARCH}
\alias{mctgarch}
\alias{plotpara}
\alias{aus2005}
\title{ MCMC estimation of GARCH (1.1) model with Student-t Distribution }
\description{ 
\code{MCMCtGARCH} uses the Markov chain Monte Carlo sampling method and Metropolis-Hastings algorithm to estimate the parameters of GARCH(1,1) model with Student-t Distribution. Simulation inefficiency factor (SIF) is calculated and shown in result which is estimated as the variance of the sample mean divided by the variance of the sample mean from a hypothetical sampler that draws independent random observations from the posterior distribution.The ineration times for burn-in period is set to m/10 in our code.
}

\usage{ 
mctgarch(y, m)
plotpara(para.matrix, m)
}

\arguments{
  \item{y}{financial return data }
  \item{m}{ MCMC sampling iteration times }
  \item{para.matrix}{returned sampled sequences of parameters}
}

\value{
   a list containing:

  \item{y}{financial return data }
  \item{m}{ MCMC sampling iteration times. m=10000 is recommended to achieve sound performance}
  \item{para.matrix}{returned sampled sequences of parameters}
  \item{aus2005}{Percentage Return of AOI stock price data from 2005 to 2008 i.e.[LN(Pt)-LN(Pt-1)]*100}
}

\author{ Ken Ren \email{ken.ren@abs.gov.au}}

\examples{

##-- 
##--   
##--  
#-- Purpose: use the Markov chain Monte Carlo sampling method 
#-- and Metropolis-Hastings algorithm to estimate the parameters
#-- of GARCH(1,1) model with Student-t Distribution. 
# 
# arguments:
#xt             Financial return data 
#m              MCMC sampling iteration times 
#para.matrix    Sampled sequences of parameters for M iterations 
#aus2005        Test data in example: Return of AOI stock price data from 2005 to 2008 
# output:
# matrix of estimated GARCH(1 ,1) parameters a0, a1, b1 and nv. 
# More information please refer to Pdf Help document

#xt<-matrix(scan(file="AUS2005.txt"),ncol=1,byrow=T)

data(aus2005)
xt<-aus2005 
m<-1000
mc=mctgarch(xt, m)
plotpara(mc, m)

}
\keyword{ distribution }
\keyword{ models }


