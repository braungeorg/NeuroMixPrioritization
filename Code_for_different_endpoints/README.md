# A Quick-Guide on how to use the compound prioritization approach using your own data and endpoints

## 1 - Select your compounds of interest
First, identify a list of compounds you would be interested in. This can be lists of compounds which were already grouped for your endpoint of interest like the lists we used in our neurotoxic approach from U.S. EPA, or every other source. Once you have identified your candidates, you need to collect the necessary data and add them to the file "physicochemical_and_kinetic_parameters.csv". More information is available in the README included in the Input folder. 

## 2 - Configure httk for your compounds of interest
Within the Code folder, you can find the script "1_httk_configuration.R". Make sure that the "httk" folder is in your working directory and run the script. The configured httk package including all data you have added to "physicochemical_and_kinetic_parameters.csv" is installed to your local R library. 

## 3 - Identify (bio)assays of interest and add raw data
Assemble assays regarding your toxicity endpoint of interest. To allow for fitting the respective concentration response curves and check for linearity, it would be optimal if you have raw data (i.e. concentrations and their responses). Add your information to the file "raw_assay_data.csv" in the Input folder. Information regarding the column titles can be found in the README within the Input folder. 

## 4 - Fit your assay data and merge endpoints per compound
Once you are finished with setting up "raw_assay_data.csv" you can start the "2_CRC_fit.R" script from the Code folder. This will both fit data with concentration-response data as well as merging specific and non-specific effects per compound. The result is the file "fitted_assay_data.csv" which is a list of compounds about which assays were merged and benchmark values like the logarithmic AC10 or AC50 with merged curve characteristics like slope and top of the curve.

## 5 - Simulate mixtures and create results
Once the steps 1 - 4 are finished and files are set up, you can run the script "3_Mixture_simulation.R". Here, you can add filters and configure httk. Configurations could be if you want to use restrictive or non-restrictive clearance, how many samples to draw from MonteCarlo simulations or the concentration quantiles used for calculating Mix_random. You can also set the number of mixtures to simulate. You have three major outputs which are written to the output folder.
 - summarized data: An overview of your final compounds considered for mixture simulation. Allows to filter for example for high SR. Note, if you select SR, to also check if the SR is based on experimental data or against Baseline toxicity. If baseline was used, you may not consider the resulting high SR.
 - Mix_median. Is a simulated mixture where all compounds of interest are considered present at the median plasma concentration quantile. The most relevant parameter is Effecti, the individual effects of the participating chemicals. The higher Effecti, the more contributing is this chemical to the overall mixture effect. 
 - Mix_random has two outputs. Summary_Mixtures gives an overview of all simulated mixtures, how many and which chemicals were needed to result in 90 % total mixture effect. In Summary_Contributors you have a list of all chemicals which at least once dominated a mixture by contributing to 90 % mixture effect. Here, you have the Frequency Factors and high FF (near to 1) represent relevant and frequent contributors whilst low values represent unlikely dominant chemicals. 

-------------------------------------------------------------------------------------------------------------------

For all files you can find a README files in the respective folders, explaining values and their meaning. 
