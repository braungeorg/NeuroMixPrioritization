##########################################
#3_Mixture_simulation
##########################################
#This script is used to fit your raw CRC data, check for outliers and linearity and merge assay data per compound
##########################################
#Necessary R packages
list.of.packages <- c("dplyr","doParallel","foreach","data.table","deSolve","doSNOW","stringr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
#call installed packages from library

library(httk)
library(dplyr)
library(doParallel)
library(foreach)
library(data.table)
library(deSolve)
library(doSNOW)
library(stringr)

#This will set your working directory to the folder where the "Code" folder with this script is contained
#You can set it freely, but need to consider to change paths accordingly so all files needed are found
#These are: "physicochemical_and_kinetic_parameters.csv" and "chem.physical_and_invitro.data_for_httk_configuration.csv"
setwd(str_remove(dirname(rstudioapi::getActiveDocumentContext()$path),"/Code"))

#########################################
#########################################

#Detect and use available cores for multiple processing ...
#cores <- detectCores()-1
#...or set manually
cores <- 60

Data_new = read.csv("Input/physicochemical_and_kinetic_parameters.csv")
Effect_data = read.csv("Input/fitted_assay_data.csv")
Effects = Effect_data[Effect_data$Viability!=T,]
Viability = Effect_data[Effect_data$Viability==T,]

################################################################################
#calculate plasma concentrations in µM
################################################################################

#Three quantiles - 10 % (low), 50 % (median), 90 % (high)
#MonteCarlo simulation using the httk-included httk-pop, which was a survey from NHANES
#Calculated with physiology data from 20,000 people

pb <- txtProgressBar(max=length(Data_new$CAS), style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

Quantiles = c(0,0.1,0.5,0.9)

cl <- makeCluster(cores)
doSNOW::registerDoSNOW(cl)
Concentrations = foreach (i= 1:length(Data_new$CAS), .combine = "rbind",.packages = "httk",.options.snow=opts) %dopar% {
  set.seed(18)
  y = tryCatch(as.numeric(calc_mc_css(samples = 20000,which.quantile = Quantiles,chem.cas = Data_new$CAS[i], suppress.messages = T,model = "pbtk",output.units = "uM", calc.analytic.css.arg.list = list(well.stirred.correction = TRUE,adjusted.Funbound.plasma = TRUE, regression = TRUE, IVIVE = NULL,restrictive.clearance = T, bioactive.free.invivo = T,daily.dose = Data_new$Intake_rate[i]))),
               error=function(e){y=c(NA,NA,NA,NA)})
  y
}
close(pb)
stopCluster(cl)

Concentrations1 <- as.data.frame(Concentrations)
Filter1 = is.na(Concentrations1$V1)
Concentrations2 = filter(Concentrations1,Filter1==F)
Data_new = filter(Data_new,Filter1==F)
Filter2 = Filter1[1:length(Concentrations2$V1)]
Concentrations2[,1] = 0
for(i in 1:length(Concentrations2$V1)){
  if(max(Concentrations2[i,])==min(Concentrations2[i,])){
    Filter2[i] = T
  }
}
Concentrations3 = filter(Concentrations2,Filter2==F)

colnames(Concentrations3) = c("Plasma_Conc_0",
                              "Plasma_Conc_10",
                              "Plasma_Conc_50",
                              "Plasma_Conc_90")

Data_new = filter(Data_new,Filter2==F)

Final_Data <- data.frame(Data_new,Concentrations3)

################################################################################
#Calculate baseline toxicity based on Lee et al. 2021
################################################################################

exps <- function(x){
  10^x
}

#Lee et al. 2021 - equation 11
calc_logKlipw <- function(logKow){
  logKow*1.01+0.12
}

Final_Data$logKlipw <- sapply(Final_Data$logKow,calc_logKlipw)

#Lee et al. 2021 - equation 14
calc_logDlipw <- function(logKlipw,Neutral_Form){
  log10(10^(logKlipw)*(Neutral_Form+(10^(-1))*(1-Neutral_Form)))
}

Final_Data$logDlipw = mapply(calc_logDlipw,Final_Data$logKlipw,Final_Data$Neutral_Form)

#Lee et al. 2021 - equation 19 with a = 1.23, b =  4.97, and c = -0.236
calc_IC10_Baseline = function(logDlipw){
  (1/(10^(1.23+4.97*(1-exp(-0.236*logDlipw)))))*1000000
}

Final_Data$IC10_Baseline = sapply(Final_Data$logDlipw,calc_IC10_Baseline)

#for calculations the IC10 has to be temporarily set to log value
Final_Data$IC10_Baseline = log10(Final_Data$IC10_Baseline)

#Check if concatenated assay data is available for each compound, otherwise use baseline toxicity
#1st - Check if specific effect data is available and use this
#2nd - Check if viability data is available and use this
#3d - Use Baseline predictions and extrapolate IC50 assuming hill equation with slope of 1
Final_Data$AC10 <- 1:length(Final_Data$CAS)
Final_Data$AC10 = as.numeric(Final_Data$AC10)
Final_Data$AC10_Top <- Final_Data$AC10
Final_Data$AC10_HillSlope <- Final_Data$AC10
Final_Data$AC50 <- Final_Data$AC10
Final_Data$Endpoint_effect <- c(1:length(Final_Data$CAS))
Final_Data$IC10 <- Final_Data$AC10
Final_Data$IC10_Top <- Final_Data$IC10
Final_Data$IC10_HillSlope <- Final_Data$IC10
Final_Data$IC50 <- Final_Data$IC10
Final_Data$Endpoint_cytotox <- Final_Data$Endpoint_effect

#Derive AC10 values
for (i in 1:length(Final_Data$CAS)){
  if(Final_Data$CAS[i] %in% Effects$CAS){
    z <- match(Final_Data$CAS[i],Effects$CAS)
    Final_Data$AC10[i] =Effects$log_AC10[z]
    Final_Data$AC10_Top[i] =Effects$Top_curve[z]
    Final_Data$AC10_HillSlope[i] =Effects$Slope[z]
    Final_Data$AC50[i] =Effects$log_AC50[z]
    Final_Data$Endpoint_effect[i] = Effects$Endpoints[z]
  } else if (Final_Data$CAS[i] %in% Viability$CAS){
    z <- match(Final_Data$CAS[i],Viability$CAS)
    Final_Data$AC10[i] = Viability$log_AC10[z]
    Final_Data$AC10_Top[i] = Viability$Top_curve[z]
    Final_Data$AC10_HillSlope[i] = Viability$Slope[z]
    Final_Data$AC50[i] = Viability$log_AC50[z]
    Final_Data$Endpoint_effect[i] = Viability$Endpoints[z]
  } else {
    Final_Data$AC10[i] = Final_Data$IC10_Baseline[i]
    Final_Data$Endpoint_effect[i] = "Baseline"
    Final_Data$AC10_Top[i] = 100
    Final_Data$AC10_HillSlope[i] = 1
    Final_Data$AC50[i] = Final_Data$IC10_Baseline[i]-log10(10/(100-10))
  }
}

#Derive IC10 values (needed for calculation of SR)
for (i in 1:length(Final_Data$CAS)){
  if (Final_Data$CAS[i] %in% Viability$CAS){
    z <- match(Final_Data$CAS[i],Viability$CAS)
    Final_Data$IC10[i] = Viability$log_AC10[z]
    Final_Data$IC10_Top[i] = Viability$Top_curve[z]
    Final_Data$IC10_HillSlope[i] = Viability$Slope[z]
    Final_Data$IC50[i] = Viability$log_AC50[z]
    Final_Data$Endpoint_cytotox[i] = Viability$Endpoints[z]
  } else {
    Final_Data$IC10[i] = Final_Data$IC10_Baseline[i]
    Final_Data$Endpoint_cytotox[i] = "Baseline"
    Final_Data$IC10_Top[i] = 100
    Final_Data$IC10_HillSlope[i] = 1
    Final_Data$IC50[i] = Final_Data$IC10_Baseline[i]-log10(10/(100-10))
  }
}

#Effects are masked if IC10 is below AC10, hence, IC10 is set as new measure for potency
for(i in 1:length(Final_Data$AC10)){
  if(Final_Data$IC10[i]<Final_Data$AC10[i]){
    Final_Data$AC10[i] = Final_Data$IC10[i]
    Final_Data$Endpoint_effect[i] = Final_Data$Endpoint_cytotox[i]
    Final_Data$AC10_Top[i] = Final_Data$IC10_Top[i]
    Final_Data$AC10_HillSlope[i] = Final_Data$IC10_HillSlope[i]
    Final_Data$AC50[i] = Final_Data$IC50[i]
  }
}

#Specificity Ratio SR - high = effects are driven by specific endpoints, not cytotoxicity
#Baseline is included, but if you filter your results for high SR you should not mix SR based on experimental and baseline values
Final_Data$SR <- exps(Final_Data$IC10)/exps(Final_Data$AC10)

#Toxic ratio (TR) is used as a measure in how more toxic experimental cytotoxicity values are compared to baseline predictions
Final_Data$TR <- exps(Final_Data$IC10_Baseline)/exps(Final_Data$IC10)

#Apply here general filter, in example exclude chemicals > 80 g/mol as we did in our approach
#Final_Data = Final_Data[Final_Data$MW>=80,]

#To include the additional information, this summarized data is written to the output
fwrite(Final_Data,"Output/summarized_data.csv",row.names = F)

#For mixture effect assessment the linear portion of the CRC is used -> logConc transformed to Conc
Final_Data$AC10 = exps(Final_Data$AC10)
Final_Data$IC10 = exps(Final_Data$IC10)
Final_Data$IC10_Baseline = exps(Final_Data$IC10_Baseline)

################################################################################
#Mixture analysis
################################################################################
#Get as separate data.frame of the calculated concentrations for later on use
Concentrations_final = data.frame(Final_Data$Compound,
                                  Final_Data$CAS,
                                  rep(0,length(Final_Data$Compound)),
                                  Final_Data$Plasma_Conc_10,
                                  Final_Data$Plasma_Conc_50,
                                  Final_Data$Plasma_Conc_90)

colnames(Concentrations_final) = c("Compound",
                                   "CAS",
                                   "0%",
                                   "10%",
                                   "50%",
                                   "90%")

################################################################################
#Calculate Mix_median assuming all chemicals present with their median concentrations
Conc <- Final_Data$Plasma_Conc_50
Ctot <- sum(Conc)
pi <- as.numeric(Conc/Ctot)
slopei <- 10/Final_Data$AC10
Sum_pi_slopei = sum(pi*slopei/10)
AC10.Mix = (1/Sum_pi_slopei)
Conc_new = pi*AC10.Mix
slopei = 10/Final_Data$AC10
pi_slopei = pi*slopei
slopem = sum(pi_slopei)
Effecti = Conc_new*slopei
Effectm = slopem*sum(Conc_new)

BigMixture <- data.frame(Final_Data$Compound,
                         Final_Data$CAS,
                         Conc_new,
                         Effecti,
                         Final_Data$AC10,
                         pi,
                         Final_Data$Endpoint_effect)

colnames(BigMixture) = c("Name","CAS","Conc","Effecti","AC10µM","pi","Endpoint(s)")
data.table::fwrite(BigMixture,"Output/Mix_median.csv",row.names = F)

################################################################################
#To address variability, mixtures are calculated where
#a) not all compounds are present and 
#b) if present, it is randomly selected if 10%, 50%, or 90% quantile
#This is calculated for "Nr_of_Mixtures" mixtures
#A frequency factor checks the overall frequency and the frequency of a compound in the top contributing compounds
#Top contributing = is within average number of compounds which are needed to cumulatively get 90 % of total effect

#All calculations are based on Escher et al. 2020

Nr_of_Mixtures <- 50
Compounds_per_Mixture <- length(Final_Data$Compound)

Mix_CAS = Final_Data$CAS

Mix_Exp = Final_Data$Expocast

Mix_Nam = Final_Data$Compound

Mix_AC10 = Final_Data$AC10

Mix_slopei = 0.1/Mix_AC10

Names = paste0("Mixture_",rep(1:Nr_of_Mixtures)) #used as header/name of each mixture

pb <- txtProgressBar(max=Nr_of_Mixtures, style=3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress=progress)

cl <- makeCluster(cores)
doSNOW::registerDoSNOW(cl)

#set.seed(18) #you can choose any number or delete this line; needed to get identical results per run (otherwise slightly different due to random functions)
#Here, you can select the probability of the quantiles of the concentrations (0%,10%,50%,90%)
#Currently, each quantile is equally likely (25 % for four possibilities)
Mix_Conc = foreach (i = 1:Nr_of_Mixtures,.combine = "cbind",.options.snow=opts)%dopar%{
  set.seed(i)
  y = 1:Compounds_per_Mixture
  for(j in 1:Compounds_per_Mixture){
    y[j] = Concentrations_final[j,sample(c(3,4,5,6),size=1,prob = c(0.25,0.25,0.25,0.25))]
  }
  y
}
stopCluster(cl)

for(c in cores:1){
  slices = Nr_of_Mixtures/c
  if(slices %% 2 == 0){
    cores = c
    break
  }
}

#To allow for fast parallel processing (depending on how many mixtures you want to simulate) a folder will be generated where your data is temporarily stored in a sliced format
slices = Nr_of_Mixtures/cores
dir.create("Output/sliced")

Mix_Conc = as.data.frame(Mix_Conc)
colnames(Mix_Conc) = Names

Ctot = 1:Nr_of_Mixtures
Ctot = colSums(Mix_Conc)

Ctot_frame = matrix(ncol=length(Ctot),nrow=Compounds_per_Mixture)
for(ct in 1:length(Ctot)){
  Ctot_frame[,ct] = Ctot[ct]
}

dir.create("Output/sliced2")

Start = 1
for(l in 1:cores){
  fwrite(Ctot_frame[,(Start:(Start+slices-1))],paste0("Output/sliced2/",l,".csv"),row.names = F)
  Start=Start+slices
}

Start = 1

for(l in 1:cores){
  fwrite(Mix_Conc[,(Start:(Start+slices-1))],paste0("Output/sliced/",l,".csv"),row.names = F)
  Start=Start+slices
}

cl <- makeCluster(cores)
doSNOW::registerDoSNOW(cl)

Mix_pi = foreach::foreach(i=1:cores,.combine = "cbind")%dopar%{
  data = as.data.frame(data.table::fread(paste0("Output/sliced/",i,".csv")))
  data2 = as.data.frame(data.table::fread(paste0("Output/sliced2/",i,".csv")))
  
  result = data
  for(j in 1:slices){
    result[,j] = data[,j]/data2[,j]
  }
  result
}
stopCluster(cl)

Start=1
for(l in 1:cores){
  fwrite(Mix_pi[,(Start:(Start+slices-1))],paste0("Output/sliced/",l,".csv"),row.names = F)
  Start=Start+slices
}

cl <- makeCluster(cores)
doSNOW::registerDoSNOW(cl)
Mix_pi_slopei = foreach::foreach(i=1:cores,.combine = "cbind")%dopar%{
  data = as.data.frame(data.table::fread(paste0("Output/sliced/",i,".csv")))
  result = data*Mix_slopei
  result
}
stopCluster(cl)

Mix_Slopem = as.numeric(colSums(Mix_pi_slopei))
Mix_AC10.Mix = as.numeric((0.1/Mix_Slopem))

Mix_AC10.Mix_frame = matrix(nrow=Compounds_per_Mixture,ncol=Nr_of_Mixtures)
for(ac in 1:ncol(Mix_AC10.Mix_frame)){
  Mix_AC10.Mix_frame[,ac] = Mix_AC10.Mix[ac]
}

Start = 1
for(l in 1:cores){
  fwrite(Mix_AC10.Mix_frame[,(Start:(Start+slices-1))],paste0("Output/sliced2/",l,".csv"),row.names = F)
  Start=Start+slices
}

cl <- makeCluster(cores)
doSNOW::registerDoSNOW(cl)

Mix_Conc = foreach::foreach(i=1:cores,.combine = "cbind")%dopar%{
  data = as.data.frame(data.table::fread(paste0("Output/sliced/",i,".csv")))
  data2 = as.data.frame(data.table::fread(paste0("Output/sliced2/",i,".csv")))
  result = data
  for(j in 1:slices){
    result[,j] = data[,j]*data2[,j]
  }
  result
}
stopCluster(cl)

Ctot = as.numeric(colSums(Mix_Conc))

Mix_Effecti = Mix_Conc
Mix_Effectm = 1:Nr_of_Mixtures

Start=1
for(l in 1:cores){
  fwrite(Mix_Conc[,(Start:(Start+slices-1))],paste0("Output/sliced/",l,".csv"),row.names = F)
  Start=Start+slices
}

cl <- makeCluster(cores)
doSNOW::registerDoSNOW(cl)

Mix_Effecti = foreach::foreach(i=1:cores,.combine = "cbind")%dopar%{
  data = as.data.frame(data.table::fread(paste0("Output/sliced/",i,".csv")))
  result = data*Mix_slopei
  result
}
stopCluster(cl)

Mix_Effectm = as.numeric(Ctot*colSums(Mix_pi_slopei,na.rm = T))

Mix_Effecti=Mix_Effecti*100
Mix_Effectm=Mix_Effectm*100

Start = 1
for(l in 1:cores){
  fwrite(Mix_Effecti[,(Start:(Start+slices-1))],paste0("Output/sliced/",l,".csv"),row.names = F)
  Start=Start+slices
}

cl <- makeCluster(cores)
doSNOW::registerDoSNOW(cl)

Ordered_Mixtures_Index = foreach::foreach(i=1:cores,.combine = "cbind")%dopar%{
  data = as.data.frame(data.table::fread(paste0("Output/sliced/",i,".csv")))
  result = data
  for(j in 1:slices){
    result[,j] = order(data[,j],decreasing = T)
  }
  result
}
stopCluster(cl)

Start = 1
for(l in 1:cores){
  fwrite(Ordered_Mixtures_Index[,(Start:(Start+slices-1))],paste0("Output/sliced2/",l,".csv"),row.names = F)
  Start=Start+slices
}

cl <- makeCluster(cores)
doSNOW::registerDoSNOW(cl)

Ordered_Mixtures_Effects = foreach::foreach(i=1:cores,.combine = "cbind")%dopar%{
  data = as.data.frame(data.table::fread(paste0("Output/sliced/",i,".csv")))
  data2 = as.data.frame(data.table::fread(paste0("Output/sliced2/",i,".csv")))
  result = data
  for(j in 1:slices){
    result[,j] = data[data2[,j],j]
  }
  result
}
stopCluster(cl)

Cumulative_frame = apply(as.data.frame(Ordered_Mixtures_Effects),2,cumsum)

Start = 1
for(l in 1:cores){
  fwrite(Cumulative_frame[,(Start:(Start+slices-1))],paste0("Output/sliced/",l,".csv"),row.names = F)
  Start=Start+slices
}

cl <- makeCluster(cores)
doSNOW::registerDoSNOW(cl)

Cumulative = foreach::foreach(i=1:cores,.combine = "c")%dopar%{
  data = as.data.frame(data.table::fread(paste0("Output/sliced/",i,".csv")))
  result = 1:ncol(data)
  for(j in 1:slices){
    result[j] = match(T,data[,j]>9)
  }
  result
}
stopCluster(cl)
unlink("Output/sliced",recursive = T)
unlink("Output/sliced2",recursive = T)

Top_Contributors_per_Mixture = matrix(ncol = Nr_of_Mixtures,nrow=max(Cumulative))

for(i in 1:Nr_of_Mixtures){
  Top_Contributors_per_Mixture[1:Cumulative[i],i] = Mix_CAS[as.numeric(c(Ordered_Mixtures_Index[1:Cumulative[i],i]))]
}

Top_Contributors_per_Mixture_Name = matrix(ncol = Nr_of_Mixtures,nrow=max(Cumulative))

for(i in 1:Nr_of_Mixtures){
  Top_Contributors_per_Mixture_Name[1:Cumulative[i],i] = Mix_Nam[as.numeric(c(Ordered_Mixtures_Index[1:Cumulative[i],i]))]
}

Main_Contributors = rep("1",Nr_of_Mixtures)
for(i in 1:Nr_of_Mixtures){
  Main_contri = Top_Contributors_per_Mixture_Name[!is.na(Top_Contributors_per_Mixture_Name[,i]),i]
  Main_Contributors[i]=paste0(Main_contri,collapse = ";")
}

Frequ_Top_Contrib_Mixtures = as.data.table(table(Top_Contributors_per_Mixture))

row_not_0 = function(x){
  return(length(x[x!=0]))
}

Frequency = apply(Mix_Conc,1,row_not_0)

Frequ_All = data.table("CAS"=Mix_CAS,"Frequ"=Frequency)

Factor_Frequency = 1:length(Frequ_Top_Contrib_Mixtures$Top_Contributors_per_Mixture)
for(i in 1:length(Frequ_Top_Contrib_Mixtures$Top_Contributors_per_Mixture)){
  k = match(Frequ_Top_Contrib_Mixtures$Top_Contributors_per_Mixture[i],Frequ_All$CAS)
  Factor_Frequency[i] = Frequ_Top_Contrib_Mixtures[i,2]/Frequ_All$Frequ[k]
}

Summary_Mixtures = data.table(Names,Cumulative,Main_Contributors)
colnames(Summary_Mixtures) = c("Mixture","CumNr90Effect","TopContributors")

Summary_Contributors = as.data.frame(matrix(ncol = 3,nrow = length(Frequ_Top_Contrib_Mixtures$Top_Contributors_per_Mixture)))

Summary_Contributors[,2] = Frequ_Top_Contrib_Mixtures$Top_Contributors_per_Mixture

for(i in 1:length(Summary_Contributors$V1)){
  k = match(Summary_Contributors$V2[i],Final_Data$CAS)
  Summary_Contributors$V1[i] = Final_Data$Compound[k]
}

Summary_Contributors[,3] = as.numeric(unlist(Factor_Frequency))

colnames(Summary_Contributors) = c("Compound","CAS","Frequency_Factor")

fwrite(Summary_Mixtures,"Output/Summary_Mixtures_Mix_random.csv",row.names = F)
fwrite(Summary_Contributors,"Output/Summary_Contributors_Mix_random.csv",row.names = F)
