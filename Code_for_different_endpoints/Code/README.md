# 1_httk_configuration.R

This script is used for setting up httk and install a local version including all information that was added to "physicochemical_and_kinetic_parameters.csv" in the *Input* folder. 

## Requirements
- R>=4.3.1
- RTools==4.0.0
Packages
- stringr==1.4.0
- devtools==2.4.3

----------------------
# 2_CRC_fit.R

This script is used for fitting raw assay data and merging different endpoints to consensus benchmark concentrations per compound. Data for processing needs to be added to "raw_assay_data.csv" in the *Input* folder. Results in "fitted_assay_data.csv" in the *Output* folder. 

## Requirements
- R>=4.3.1
Packages
- dplyr==1.0.9
- drc==3.0-1
- tcpl==2.1.0
- EnvStats==2.7.0
- data.table==1.14.2
- stringr==1.4.0
- httk==2.0.4 (configured with own data)

----------------------
# 3_Mixture_simulation.R

This script is used for simulating mixtures considering parameters included in "physicochemical_and_kinetic_paramters.csv" and "fitted_assay_data.csv" from *Input* to calculate different prioritization approaches. Results in "Mix_median.csv", "summarized_data.csv", "Summary_Contributors_Mix_random.csv", and "Summary_Mixtures_Mix_random.csv".

## Requirements
- R>=4.3.1
Packages
- dplyr==1.0.9
- data.table==1.14.2
- stringr==1.4.0
- httk==2.0.4 (configured with own data)
- doParallel==1.0.17
- foreach==1.5.2
- deSolve==1.33
- doSNOW==1.0.20

----------------------
# 3_Mixture_simulation_without_parallelization.R

This script is used for simulating mixtures considering parameters included in "physicochemical_and_kinetic_paramters.csv" and "fitted_assay_data.csv" from *Input* to calculate different prioritization approaches. Results in "Mix_median.csv", "summarized_data.csv", "Summary_Contributors_Mix_random.csv", and "Summary_Mixtures_Mix_random.csv".

Adjusted so that no parallel processing is necessary (for small simulations). 

## Requirements
- R>=4.3.1
Packages
- dplyr==1.0.9
- data.table==1.14.2
- stringr==1.4.0
- httk==2.0.4 (configured with own data)
- deSolve==1.33

 


