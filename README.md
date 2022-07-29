# NeuroMixPrioritization
Code for "Prioritization of mixtures of neurotoxic chemicals for biomonitoring using high-throughput toxicokinetics and mixture toxicity modeling"

# Code
**1_httk_configuration.R**  
Used for setting up the httk R package 

**2_CRC_fit.R**  
Uses available CRC data and applies a fit. Further, checks for linearity up to 10 % effect.

**3_Mixture_simulation.R**  
Here, plasma concentration quantiles are calculated, effect data is merged, baseline toxicity is predicted and mixtures are simulated. Mixmedian and Mixrandom are calculated with and without blood-brain-barrier permeability restriction. 

# httk
httk with version 2.0.4 from https://github.com/USEPA/CompTox-ExpoCast-httk/releases/tag/2.0.4  
Needed for pre-configuration

# SI_Tables
**Table_B1:**  
Chemical identifiers from CompTox Dashboard for all compounds of interest.

**Table_B2:**  
Relevant properties and values for chemicals of interest with their sources.
MW = molecular weight in g/mol
Expocast = ExpoCast value in mg/kg/day
Fub = fraction unbound (plasma)
Clint = intrinsic clearance in µL/min/10^6 cells 
CNS = predicted Blood-Brain-Barrier permeability
Neutral_Form = portion of uncharged species at pH = 7.4

**Table_B3:**  
Experimental data before own fit
Activity_level = effect level relative to reference in %
log_Benchmark_concentration = log10 of concentration at activity level
Concentrations = used concentrations in µM
Responses = effects relative to reference in %
Endpoint = used assay

**Table_B4:**  
Concatenated experimental data per compound after fit
log_AC10 = log10 of the activity concentration at 10 % maximum effect in µM
log_AC50 = log10 of the activity concentration at 50 % maximum effect in µM
Slope = slope of the fitted curve
Top_curve = top effect level of fitted curve
Endpoints = all merged endpoints
Viability = if endpoints are based on cell viability (TRUE) or specific effects (FALSE)

**Table_B5:**  
Summary of many variables of interest per compound
MW = molecular weight in g/mol
Expocast = ExpoCast value in mg/kg/day
Fub = fraction unbound (plasma)
Clint = intrinsic clearance in µL/min/10^6 cells 
CNS = predicted Blood-Brain-Barrier permeability
Neutral_Form = portion of uncharged species at pH = 7.4
Plasma_Conc_X = predicted plasma concentration in µM at quantile X
logKlipw = liposome-water partition constant
logDlipw = liposome-water distribution ratio
IC10_Baseline = concentration causing 10 % effect in µM; calculated by using the baseline toxicity model
AC10 = concentration causing 10 % (neurotoxic) effect in µM; calculated from merged assay data fit
AC10_Top = top effect level of fitted curve used for AC10 in %
AC10_HillSlope = slope of fitted curve used for AC10
EC50 = concentration causing 50 % (neurotoxic) effect in µM; calculated from merged assay data fit
Endpoint_effect = merged endpoints from the fitted (neurotoxic) assays
IC10 = concentration causing 10 % (cytotoxic) effect in µM; calculated from merged assay data fit
IC10_Top = top effect level of fitted curve used for IC10 in %
IC10_HillSlope = slope of fitted curve used for IC10
IC50 = concentration causing 50 % (cytotoxic) effect in µM; calculated from merged assay data fit
Endpoint_cytotoxicity = merged endpoints from the fitted (cytotoxic) assays
SR = specificity ratio (=IC10/AC10)

**Table_B6:**  
Summary of Mixmedian without CNS restriction
Conc. = predicted plasma concentration in µM at 50 % quantile
Effecti = effect per compound in %
AC10µM = concentration causing 10 % effect in µM
pi = concentration fraction
Endpoint(s) = merged endpoints used

**Table_B7:**  
Summary of Mixmedian with CNS restriction
Conc. = predicted plasma concentration in µM at 50 % quantile
Effecti = effect per compound in %
AC10µM = concentration causing 10 % effect in µM
pi = concentration fraction
Endpoint(s) = merged endpoints used

**Table_B8:**  
Overview of mixtures in Mixrandom without CNS restriction
Mixture = Name of mixture
Effectm = mixture effect
CumNr90Effect = how many compounds were needed to result in 90 % total mixture effect
TopContributors = respective compounds neededto result in 90 % total mixture effect

**Table_B9:**  
Summary of top contributors per mixture in Mixrandom without CNS restriction
Frequency_Factor = factor of how often the chemical contributed to causing 90 % of total mixture effect relative to its overall ocurrence 

**Table_B10:**  
Overview of mixtures in Mixrandom with CNS restriction
Mixture = Name of mixture
Effectm = mixture effect
CumNr90Effect = how many compounds were needed to result in 90 % total mixture effect

**Table_B11:**  
Summary of top contributors per mixture in Mixrandom with CNS restriction
Frequency_Factor = factor of how often the chemical contributed to causing 90 % of total mixture effect relative to its overall ocurrence 
TopContributors = respective compounds neededto result in 90 % total mixture effect

**Table_B12:**  
Summary of Mixmedian without CNS restriction (only experimental values)
Conc. = predicted plasma concentration in µM at 50 % quantile
Effecti = effect per compound in %
AC10µM = concentration causing 10 % effect in µM
pi = concentration fraction
Endpoint(s) = merged endpoints used

**Table_B13:**  
Summary of Mixmedian with CNS restriction (only experimental values)
Conc. = predicted plasma concentration in µM at 50 % quantile
Effecti = effect per compound in %
AC10µM = concentration causing 10 % effect in µM
pi = concentration fraction
Endpoint(s) = merged endpoints used

**Table_B14:**  
Overview of mixtures in Mixrandom without CNS restriction (only experimental values)
Mixture = Name of mixture
Effectm = mixture effect
CumNr90Effect = how many compounds were needed to result in 90 % total mixture effect
TopContributors = respective compounds neededto result in 90 % total mixture effect

**Table_B15:**  
Summary of top contributors per mixture in Mixrandom without CNS restriction (only experimental values)
Frequency_Factor = factor of how often the chemical contributed to causing 90 % of total mixture effect relative to its overall ocurrence 

**Table_B16:**  
Overview of mixtures in Mixrandom with CNS restriction (only experimental values)
Mixture = Name of mixture
Effectm = mixture effect
CumNr90Effect = how many compounds were needed to result in 90 % total mixture effect

**Table_B17:**  
Summary of top contributors per mixture in Mixrandom with CNS restriction (only experimental values)
Frequency_Factor = factor of how often the chemical contributed to causing 90 % of total mixture effect relative to its overall ocurrence 
TopContributors = respective compounds neededto result in 90 % total mixture effect

**Table_B18:**  
Summary of Mixmedian without CNS restriction (only baseline predictions)
Conc. = predicted plasma concentration in µM at 50 % quantile
Effecti = effect per compound in %
AC10µM = concentration causing 10 % effect in µM
pi = concentration fraction
Endpoint(s) = merged endpoints used

**Table_B19:**  
Summary of Mixmedian with CNS restriction (only baseline predictions)
Conc. = predicted plasma concentration in µM at 50 % quantile
Effecti = effect per compound in %
AC10µM = concentration causing 10 % effect in µM
pi = concentration fraction
Endpoint(s) = merged endpoints used

**Table_B20:**  
Overview of mixtures in Mixrandom without CNS restriction (only baseline predictions)
Mixture = Name of mixture
Effectm = mixture effect
CumNr90Effect = how many compounds were needed to result in 90 % total mixture effect
TopContributors = respective compounds neededto result in 90 % total mixture effect

**Table_B21:**  
Summary of top contributors per mixture in Mixrandom without CNS restriction
Frequency_Factor = factor of how often the chemical contributed to causing 90 % of total mixture effect relative to its overall ocurrence 

**Table_B22:**  
Overview of mixtures in Mixrandom with CNS restriction (only baseline predictions)
Mixture = Name of mixture
Effectm = mixture effect
CumNr90Effect = how many compounds were needed to result in 90 % total mixture effect

**Table_B23:**  
Summary of top contributors per mixture in Mixrandom with CNS restriction (only baseline predictions)
Frequency_Factor = factor of how often the chemical contributed to causing 90 % of total mixture effect relative to its overall ocurrence 
TopContributors = respective compounds neededto result in 90 % total mixture effect
