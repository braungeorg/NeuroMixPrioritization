% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_hep_fu.R
\name{calc_hep_fu}
\alias{calc_hep_fu}
\title{Calculate the free chemical in the hepaitic clearance assay}
\usage{
calc_hep_fu(
  chem.cas = NULL,
  chem.name = NULL,
  dtxsid = NULL,
  parameters = NULL,
  Vr = 0.005,
  pH = 7.4
)
}
\arguments{
\item{chem.cas}{Chemical Abstract Services Registry Number (CAS-RN) -- if
parameters is not specified then the chemical must be identified by either
CAS, name, or DTXISD}

\item{chem.name}{Chemical name (spaces and capitalization ignored) --  if
parameters is not specified then the chemical must be identified by either
CAS, name, or DTXISD}

\item{dtxsid}{EPA's 'DSSTox Structure ID (\url{https://comptox.epa.gov/dashboard})
 -- if parameters is not specified then the chemical must be identified by 
either CAS, name, or DTXSIDs}

\item{parameters}{Parameters from the appropriate parameterization function
for the model indicated by argument model}

\item{Vr}{Rratio of cell volume to incubation volume. Default is taken from}

\item{pH}{pH of the incupation medium.}
}
\value{
A numeric fraction between zero and one
}
\description{
Method from Kilford et al. (2008) for fraction of unbound chemical in the 
 hepatocyte intrinsic clearance assay
}
\references{
Kilford, Peter J., et al. "Hepatocellular binding of drugs: 
correction for unbound fraction in hepatocyte incubations using microsomal 
binding or drug lipophilicity data." Drug Metabolism and Disposition 36.7 
(2008): 1194-1197.

Wetmore, Barbara A., et al. "Incorporating high-throughput exposure 
predictions with dosimetry-adjusted in vitro bioactivity to inform chemical 
toxicity testing." Toxicological Sciences 148.1 (2015): 121-136.
}
\author{
John Wambaugh and Robert Pearce
}
\keyword{in-vitro}
