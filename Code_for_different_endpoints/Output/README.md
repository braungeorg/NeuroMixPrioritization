
# summarized_data.csv

An overview of compounds included in mixture simulation with information about physicochemical and kinetic properties, effect data benchmark values.

## Columns: 
- Compound = Name of the chemical
- CAS = CAS number of the chemical
- MW = molecular weight in g/mol of the chemical
- Intake_rate = intake_rate / exposure in mg/kg/day 
- Fub = fraction unbound (plasma)
- Clint = intrinsic clearance in µL/min/10^6 cells
- pKa_Acceptor/_Donor = pKa values of the chemical
- logKow = logKow of the chemical
- Neutral_From = fraction of uncharged species at pH 7.4
- ..._Reference = the source for the respective value
- Plasma_Conc_0 = plasma concentration in µM at 0 % quantile from MonteCarlo prediction
- Plasma_Conc_10 = plasma concentration in µM at 10 % quantile from MonteCarlo prediction
- Plasma_Conc_50 = plasma concentration in µM at 50 % quantile from MonteCarlo prediction
- Plasma_Conc_90 = plasma concentration in µM at 90 % quantile from MonteCarlo prediction
- logKlipw = liposome-water partitioning coefficient of the chemical
- logDlipw = liposome-water distribution coefficient of the chemical
- IC10_Baseline = calculated logarithmic concentration in µM needed to result in 10 % effect based on the baseline model
- AC10 = logarithmic activity concentration in µM resulting in 10 % effect, merged if multiple endpoints
- AC10_Top = % of top of the merged curve for specific activity concentrations
- AC10_Hillslope = hill slope of the merged curve for specific activity concentrations
- AC50 = logarithmic activity concentration in µM resulting in 50 % effect, merged if multiple endpoints
- Endpoint_effect = endpoint or merged endpoints used from endpoint-specific assays, i.e. receptor binding
- IC10 = logarithmic inhibitory concentration in µM resulting in 10 % effect, merged if multiple endpoints
- IC10_Top = % of top of the merged curve for unspecific inhibitory concentrations
- IC10_Hillslope = hill slope of the merged curve for unspecific inhibitory concentrations
- IC50 = logarithmic inhibitory concentration in µM resulting in 50 % effect, merged if multiple endpoints
- Endpoint_cytotoxicity = endpoint or merged endpoints used from non-endpoint-specifc assays, i.e. viability assays
- SR = specificity ratio
- TR = toxicity ratio

If you sort by highest SR, you get specifically acting chemicals for your endpoint of interest. Note, that if baseline was used as form of cytotoxicity, the SR can be very large and should be used with caution. 

-----------------------------------------------------------------------------------

# Mix_median.csv

Result of a mixture where all compounds of interest are present at median predicted plasma concentrations up to an effect of 10 % total.

## Columns: 
- Name = Name of the chemical
- CAS = CAS of the chemical
- Conc = predicted plasma concentration in µM
- Effecti = individual effect in % of the respective chemical
- AC10µM = the µM activity concentration resulting in 10 % of the individual chemical
- pi = concentration fraction of the respective chemical within this mixture
- Endpoint(s) = which endpoints were considered for potency of this chemical

If you filter by highest Effecti, you can select chemicals which contribute highly to the mixture effect. 

-------------------------------------------------------------------------------------

# Summary_Mixtures_Mix_random.csv

A summary of all mixtures regarding number and identify of dominant chemicals. 

## Columns: 
- Mixture = Name/Nr of the mixture
- CumNr90Effect = Number of chemicals needed to explain 90 % of mixture effect
- TopContributors = Names of chemicals needed to explain 90 % of mixture effect

-------------------------------------------------------------------------------------

# Summary_Contributors_Mix_random.csv

A summary of all chemicals which at least once were contributing to 90 % mixture effect. 

## Columns: 
- Compound = Name of the chemical
- CAS = CAS of the chemical
- Frequency_Factor = frequency factor of the chemical over all simulated mixtures






 


