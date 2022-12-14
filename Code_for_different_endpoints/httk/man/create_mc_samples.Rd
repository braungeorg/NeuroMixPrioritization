% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_mc_samples.R
\name{create_mc_samples}
\alias{create_mc_samples}
\title{Create a data table of draws of parameter values for Monte Carlo}
\usage{
create_mc_samples(
  chem.cas = NULL,
  chem.name = NULL,
  dtxsid = NULL,
  parameters = NULL,
  samples = 1000,
  species = "Human",
  suppress.messages = FALSE,
  model = "3compartmentss",
  httkpop = TRUE,
  invitrouv = TRUE,
  calcrb2p = TRUE,
  censored.params = list(),
  vary.params = list(),
  return.samples = FALSE,
  tissue = NULL,
  httkpop.dt = NULL,
  invitro.mc.arg.list = list(adjusted.Funbound.plasma = TRUE, poormetab = TRUE,
    fup.censored.dist = FALSE, fup.lod = 0.01, fup.meas.cv = 0.4, clint.meas.cv = 0.3,
    fup.pop.cv = 0.3, clint.pop.cv = 0.3),
  httkpop.generate.arg.list = list(method = "direct resampling", gendernum = NULL,
    agelim_years = NULL, agelim_months = NULL, weight_category = c("Underweight",
    "Normal", "Overweight", "Obese"), gfr_category = c("Normal", "Kidney Disease",
    "Kidney Failure"), reths = c("Mexican American", "Other Hispanic",
    "Non-Hispanic White", "Non-Hispanic Black", "Other")),
  convert.httkpop.arg.list = list(),
  propagate.invitrouv.arg.list = list(),
  parameterize.arg.list = list(restrictive.clearance = T, default.to.human = FALSE,
    clint.pvalue.threshold = 0.05, regression = TRUE)
)
}
\arguments{
\item{chem.cas}{Chemical Abstract Services Registry Number (CAS-RN) -- if
parameters is not specified then the chemical must be identified by either
CAS, name, or DTXISD}

\item{chem.name}{Chemical name (spaces and capitalization ignored) --  if
parameters is not specified then the chemical must be identified by either
CAS, name, or DTXISD}

\item{dtxsid}{EPA's DSSTox Structure ID (\url{https://comptox.epa.gov/dashboard})
 -- if parameters is not specified then the chemical must be identified by 
either CAS, name, or DTXSIDs}

\item{parameters}{Parameters from the appropriate parameterization function
for the model indicated by argument model}

\item{samples}{Number of samples generated in calculating quantiles.}

\item{species}{Species desired (either "Rat", "Rabbit", "Dog", "Mouse", or
default "Human"). Species must be set to "Human" to run httkpop model.}

\item{suppress.messages}{Whether or not to suppress output message.}

\item{model}{Model used in calculation: 'pbtk' for the multiple compartment
model,'3compartment' for the three compartment model, '3compartmentss' for
the three compartment steady state model, and '1compartment' for one
compartment model.  This only applies when httkpop=TRUE and species="Human",
otherwise '3compartmentss' is used.}

\item{httkpop}{Whether or not to use the Ring et al. (2017) "httkpop"
population generator. Species must be 'Human'.}

\item{invitrouv}{Logical to indicate whether to include in vitro parameters
such as intrinsic hepatic clearance rate and fraction unbound in plasma
in uncertainty and variability analysis}

\item{calcrb2p}{Logical determining whether or not to recalculate the 
chemical ratio of blood to plasma}

\item{censored.params}{The parameters listed in censored.params are sampled
from a normal distribution that is censored for values less than the limit
of detection (specified separately for each paramter). This argument should
be a list of sub-lists. Each sublist is named for a parameter in
"parameters" and contains two elements: "CV" (coefficient of variation) and
"LOD" (limit of detection, below which parameter values are censored. New
values are sampled with mean equal to the value in "parameters" and standard
deviation equal to the mean times the CV.  Censored values are sampled on a
uniform distribution between 0 and the limit of detection. Not used with
httkpop model.}

\item{vary.params}{The parameters listed in vary.params are sampled from a
normal distribution that is truncated at zero. This argument should be a
list of coefficients of variation (CV) for the normal distribution. Each
entry in the list is named for a parameter in "parameters". New values are
sampled with mean equal to the value in "parameters" and standard deviation
equal to the mean times the CV. Not used with httkpop model.}

\item{return.samples}{Whether or not to return the vector containing the
samples from the simulation instead of the selected quantile.}

\item{tissue}{Desired steady state tissue conentration.}

\item{httkpop.dt}{A data table generated by \code{\link{httkpop_generate}}.
This defaults to NULL, in which case \code{\link{httkpop_generate}} is 
called to generate this table.}

\item{invitro.mc.arg.list}{Additional parameters passed to 
\code{\link{invitro_mc}}.}

\item{httkpop.generate.arg.list}{Additional parameters passed to 
\code{\link{httkpop_generate}}.}

\item{convert.httkpop.arg.list}{Additional parameters passed to the 
convert_httkpop_* function for the model.}

\item{propagate.invitrouv.arg.list}{Additional parameters passed to model's
associated in vitro uncertainty and variability propagation function}

\item{parameterize.arg.list}{Additional parameters passed to the 
parameterize_* function for the model.}
}
\value{
A data table where each column corresponds to parameters needed for the 
specified model and each row represents a different Monte Carlo sample of
parameter values.
}
\description{
This function creates a data table of draws of parameter values for use with 
Monte Carlo methods
}
\examples{

\donttest{
sample_set = create_mc_samples(chem.name = 'bisphenol a')
}

}
\references{
Wambaugh, John F., et al. "Toxicokinetic triage for 
environmental chemicals." Toxicological Sciences 147.1 (2015): 55-67.

Ring, Caroline L., et al. "Identifying populations sensitive to
environmental chemicals by simulating toxicokinetic variability."
Environment international 106 (2017): 105-118.
}
\author{
Caroline Ring, Robert Pearce, and John Wambaugh
}
\keyword{Monte-Carlo}
