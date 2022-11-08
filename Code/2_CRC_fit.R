library(httk)
library(dplyr)
library(drc)
library(tcpl)
library(EnvStats)
library(data.table)
library(stringr)

setwd(str_remove(dirname(rstudioapi::getActiveDocumentContext()$path),"/Code"))


#Detect and use available cores for multiple processing or set manually
#cores <- detectCores()-1
cores <- 60

exps <- function(x){
  10^x
}

#Import in-vitro effect and cytotoxicity data form literature, ICE, NTP
Data <- fread("SI_Tables/Table_B3.csv")

Data = Data[!is.null(Data$log_Benchmark_concentration),]

#Check if each concentration value has a correspondant response value
for(i in 1:length(Data$Compound)){
  a = length(as.numeric(unlist(strsplit(Data$Concentrations[i],split=","))))
  b = length(as.numeric(unlist(strsplit(Data$Responses[i],split=","))))
  if(a != b){
    Data$log_Benchmark_concentration[i] = NA
  }
}

Data = Data[!is.na(Data$log_Benchmark_concentration),]

#To allow for a fit, the concentration and response values are saved as data frames in a list
Data_list = as.list(1:length(Data$Compound))

for(i in 1:length(Data$Compound)){
  if(Data$Concentrations[i]=="none"){
    Data_list[[i]] = "noCRC"
    next
  }
  Data_list[[i]] = data.frame("Dose"=as.numeric(unlist(strsplit(Data$Concentrations[i],split=","))),
                         "Response"=as.numeric(unlist(strsplit(Data$Responses[i],split=","))))
}

Negative_CRC = rep(F,length(Data$Compound))

#Warnings occur since the data.frames have length >1
for(i in 1:length(Data_list)){
  if(Data_list[[i]]=="noCRC"){
    next
  }
  frame = as.data.frame(Data_list[[i]])
  if(max(frame$Response)<=0){
    Negative_CRC[i] = T
  }
}

for(i in 1:length(Negative_CRC)){
  if(Negative_CRC[i]==T){
    Data_list[[i]]$Response = Data_list[[i]]$Response*-1
  }
}

#Columns are pre-defined
Data$AC50_fit = NA
Data$Slope_fit = NA
Data$Top_fit = NA
Data$Bottom_fit = NA
#List which contains the respective fits
Data_fit = Data_list

#Fits data which has CRC data:
#-Rosner test used to check two most extreme values as outlieres
#-Curve needs at least two data points at different concentrations above and below 50 % relative effect
#-If top exceeds 115, a new fit with a fixed top of 100 is tried
#-If no CRC is available: Bottom = 0, Top = 100, Slope = 1, AC50 = value extrapolated from given benchmark value using hill equation

Data_safe = Data
Data_list_safe = Data_list
Data_fit_safe = Data_fit

Data = Data_safe
Data_list = Data_list_safe
Data_fit = Data_fit_safe

for (i in 1:length(Data$Compound)){
  if(Data_list[[i]]!="noCRC"){
    z = as.data.frame(Data_list[[i]])
    z$Dose = log10(z$Dose) #since we are using a log-logistic CRC model, concentrations need to be logarithmic
    z_original = z
    if(var(z$Response) == 0){
      x = NA} else {
        #Combine replicates
        replicates_table = as.data.frame(table(z$Dose))
        replicates = replicates_table$Freq[1]
        if(replicates>1){
          u = replicates-1
          resp_new = 1:(length(z$Dose)/replicates)
          dose_new = z$Dose[seq(1,length(z$Dose),by=replicates)]
          for(d in 1:length(dose_new)){
            a = z[z$Dose==dose_new[d],]
            resp_new[d] = mean(a$Response)
          }
          z = z[1:length(dose_new),]
          z$Dose = dose_new
          z$Response = resp_new
        }
        z$rel = z$Response/max(z$Response)*100
        if(max(z$Response)==0){
          x = NA
          Data$AC50_fit[i] = NA
          Data_fit[[i]] = x
          Data_list[[i]] = z
          next
        }
        #Check if the form of the equation is promising
        Filter = z$rel == 100
        if(is.nan(z$rel[1])){
          x = NA} else {
            coord_max = max(which(Filter))
            for(j in coord_max:length(z$rel)){
              if(z$rel[j]<80){
                z$rel[j] = NaN
              }
            }
            Filter = is.nan(z$rel)
            z = filter(z,Filter==F)
            above = z$rel[z$rel>50]
            below = z$rel[z$rel<50]
          }
        if(length(above)<2 | length(below)<2){
          x = NA
          Data$AC50_fit[i] = NA
          Data_fit[[i]] = x
          Data_list[[i]] = z
          next
        }
        #For fit and outliers - use original data but within accepted margins
        z = z_original[z_original$Dose%in%z$Dose,]
        #Outlier removal
        rosner = rosnerTest(z$Response,k=2)
        resulta = as.data.frame(rosner$all.stats)
        Outliers = resulta$Outlier == T
        Outliers = filter(resulta,Outliers==T)
        for (k in 1:length(Outliers$Outlier)){
          r = Outliers$Obs.Num[k]
          z$Dose[r] = NaN
        }
        Filter = is.nan(z$Dose)
        z = filter(z,Filter==F)
        z$rel = z$Response/max(z$Response)*100
        if(max(z$Response)==0){
          x = NA
          Data$AC50_fit[i] = NA
          Data_fit[[i]] = x
          Data_list[[i]] = z
          next
        }
        x <- tryCatch(drm(Response~Dose, data=z, fct=LL.4(fixed = c(NA,0,NA,NA),names = c("Slope", "Lower Limit", "Upper Limit", "AC50")), logDose = 10),
                            error=function(e){x = NA})
          }
    if(is.na(x)){
      Data$AC50_fit[i] = NA
    } else {
      y = as.data.frame(x$coefficients)
      if(y$`x$coefficients`[1]<0){
        Data$Slope_fit[i] = (-1)*y$`x$coefficients`[1]
        x$coefficients[1] = (-1)*y$`x$coefficients`[1]
      } else {Data$Slope_fit[i] = y$`x$coefficients`[1]}
      Data$Bottom_fit[i] = 0
      Data$Top_fit[i] = y$`x$coefficients`[2]
      Data$AC50_fit[i] = log10(y$`x$coefficients`[3])
      if(Data$Top_fit[i]>115){
        x <- tryCatch(drm(Response~Dose, data=z, fct=LL.4(fixed = c(NA,0,100,NA),names = c("Slope", "Lower Limit", "Upper Limit", "AC50")), logDose = 10),
                      error=function(e){x=NA})
        if(is.na(x)){
          Data$AC50_fit[i] = NA
        } else {
          y = as.data.frame(x$coefficients)
          if(y$`x$coefficients`[1]<0){
            Data$Slope_fit[i] = (-1)*y$`x$coefficients`[1]
            x$coefficients[1] = (-1)*y$`x$coefficients`[1]
          } else {Data$Slope_fit[i] = y$`x$coefficients`[1]}
          Data$Bottom_fit[i] = 0
          Data$Top_fit[i] = 100
          Data$AC50_fit[i] = log10(y$`x$coefficients`[2])
        }}
    }
    Data_fit[[i]] = x
    Data_list[[i]] = z
    }
}

for(i in 1:length(Data$Compound)){
  if(Data_list[[i]]=="noCRC"){
    Data$Top_fit[i] = 100
    Data$Slope_fit[i] = 1
    Data$Bottom_fit[i] = 0
    Data$AC50_fit[i] = Data$log_Benchmark_concentration[i] - log(Data$Activity_level[i]/(100-Data$Activity_level[i]))
  }
}

#Filter extreme values
Filter = Data$Slope_fit > 10 | Data$Slope_fit < 0.5 | Data$Bottom_fit == Data$Top_fit | is.na(Data$AC50_fit)

Data = dplyr::filter(Data,Filter==F)

for(i in 1:length(Filter)){
  if(Filter[i]==T){
    Data_fit[[i]] = NA
    Data_list[[i]] = NA
  }
}

Data_fit = Data_fit[!is.na(Data_fit)]
Data_list = Data_list[!is.na(Data_list)]

Data$AC10 = NA

for(i in 1:length(Data$Compound)){
  Data$AC10[i] = tcplHillACXX(XX=10,
                              tp=Data$Top_fit[i],
                              ga=Data$AC50_fit[i],
                              gw=Data$Slope_fit[i],
                              bt=Data$Bottom_fit[i])
}

Filter = is.nan(Data$AC10) | Data$AC10 == Inf | Data$AC10 == -Inf | is.na(Data$AC10)
Data = filter(Data,Filter==F)
for(i in 1:length(Filter)){
  if(Filter[i]==T){
    Data_fit[[i]] = NA
    Data_list[[i]] = NA
  }
}

Data_fit = Data_fit[!is.na(Data_fit)]
Data_list = Data_list[!is.na(Data_list)]

Effects = Data[Data$Is_viability==F,]
Effects_CAS = unique(Effects$CAS)
Viability = Data[Data$Is_viability==T,]
Viability_CAS = unique(Viability$CAS)

AC10然 = rep(NA,length(Effects_CAS))
Endpoint_Effect = rep(NA,length(Effects_CAS))
AC10然_Top = rep(NA,length(Effects_CAS))
AC10然_slope = rep(NA,length(Effects_CAS))
AC10然_AC50 = rep(NA,length(Effects_CAS))

for (i in 1:length(Effects_CAS)){
  Filter = Effects$CAS == Effects_CAS[i]
  A = filter(Effects,Filter==T)
  AC10然[i] = median(A$AC10)
  if(length(A$Endpoint)>1){
    Endpoint_Effect[i] = paste(A$Endpoint,collapse = "|")
  } else {Endpoint_Effect[i] = A$Endpoint[1]}
  AC10然_Top[i] = median(A$Top)
  AC10然_slope[i] = median(A$Slope)
  AC10然_AC50[i] = median(A$AC50)
}

IC10然 = rep(NA,length(Viability_CAS))
Endpoint_Cytotox = rep(NA,length(Viability_CAS))
IC10然_Top = rep(NA,length(Viability_CAS))
IC10然_slope = rep(NA,length(Viability_CAS))
IC10然_AC50 = rep(NA,length(Viability_CAS))

for (i in 1:length(Viability_CAS)){
  Filter = Viability$CAS == Viability_CAS[i]
  A = filter(Viability,Filter==T)
  IC10然[i] = median(A$AC10)
  if(length(A$Endpoint)>1){
    Endpoint_Cytotox[i] = paste(A$Endpoint,collapse = "|")
  } else {Endpoint_Cytotox[i] = A$Endpoint[1]}
  IC10然_Top[i] = median(A$Top)
  IC10然_slope[i] = median(A$Slope)
  IC10然_AC50[i] = median(A$AC50)
}

Effects_summarized <- data.frame(Effects_CAS,AC10然,AC10然_AC50,AC10然_slope,AC10然_Top,Endpoint_Effect)
Effects_summarized$Viability = F
Viability_summarized <- data.frame(Viability_CAS,IC10然,IC10然_AC50,IC10然_slope,IC10然_Top,Endpoint_Cytotox)
Viability_summarized$Viability = T

colnames(Effects_summarized) = c("CAS","AC10","AC50","Slope","Top_curve","Endpoints","Viability")
colnames(Viability_summarized) = colnames(Effects_summarized)
Activit_fit = rbind(Effects_summarized,Viability_summarized)
Activit_fit$Name = "unknown"
for(i in 1:length(Activit_fit$CAS)){
  k = match(Activit_fit$CAS[i],Data$CAS)
  Activit_fit$Name[i] = Data$Compound[k]
}

Activitydata_fit = data.frame(Activit_fit[,8],Activit_fit[,1:7])
colnames(Activitydata_fit) = c("Compound","CAS","log_AC10","log_AC50","Slope","Top_curve","Endpoints","Viability")

data.table::fwrite(Activitydata_fit,"SI_Tables/Table_B4.csv")

############################################################################################
#Dose response curves - does our linear approach below 10 % make sense?
Data_with_CRC <- Data[Data$Concentrations!="none",]
Data_with_CRC_list = Data_list[Data_list!="noCRC"]
Data_with_CRC_fit = Data_fit[Data_fit!="noCRC"]

Data_with_CRC$Top_fit <- as.numeric(Data_with_CRC$Top_fit)
Data_with_CRC$Bottom_fit <- as.numeric(Data_with_CRC$Bottom_fit)
Data_with_CRC$Slope_fit <- as.numeric(Data_with_CRC$Slope_fit)
Data_with_CRC$AC50_fit <- as.numeric(Data_with_CRC$AC50_fit)
Curves <- as.list(1:length(Data_with_CRC$Name))
SS_good_pred <- as.numeric(1:length(Data_with_CRC$Name))
SS_good_exp = SS_good_pred
SS_curve_pred = SS_good_pred
SS_curve_exp = SS_curve_pred
slopes = SS_good_pred

for (i in 1:length(Data_with_CRC$Compound)){
  Concentration <- seq(0, exps(Data_with_CRC$AC10[i]), by=(exps(Data_with_CRC$AC10[i])/100))
  Response <- Concentration
  for (j in 1:length(Response)){
    Response[j] = tcplHillVal(logc = log10(Concentration[j]), tp=Data_with_CRC$Top_fit[i],bt=Data_with_CRC$Bottom_fit[i],ga=Data_with_CRC$AC50_fit[i],gw=Data_with_CRC$Slope_fit[i])
  }
  x = data.frame(Concentration,Response)
  colnames(x) = c("Concentration","Response")
  x$Response = (x$Response/Data_with_CRC$Top_fit[i])*100
  Filter = is.nan(x$Concentration)
  x = filter(x,Filter==F)
  x$CIu = 1.1*x$Response
  x$CIb = 0.9*x$Response
  slope_linear = 10/tail(x$Concentration,1)
  slopes[i] = slope_linear
  x$linear = slope_linear*x$Concentration
  SS_good_pred[i] = sum((x$Response-x$CIu)^2)
  SS_curve_pred[i] = sum((x$Response-x$linear)^2)
  #Experimental values
  y = as.data.frame(Data_with_CRC_list[[i]])
  y$Dose = 10^y$Dose
  Filter = is.nan(y$Dose) | is.nan(y$Response) | y$rel >30 | y$rel <0 
  y = filter(y,Filter==F)
  if(length(y$Dose)<4){
    SS_good_exp[i] = Inf
    SS_curve_exp[i] = Inf
    next
  }
  model = tryCatch(lm(y$rel~0+y$Dose),error=function(e){NA})
  if(is.na(model)){
    Curves[[i]] = as.data.frame(x)
    SS_good_exp[i] = Inf
    SS_curve_exp[i] = Inf
    next
  }
  slope_exp = as.numeric(model$coefficients)
  x$exp = x$Concentration*slope_exp
  x$CIexpu = 1.1*x$exp
  x$CIexpb = 0.9*x$exp
  SS_good_exp[i] = sum((x$exp-x$CIexpu)^2)
  SS_curve_exp[i] = sum((x$exp-x$linear)^2)
  Curves[[i]] = as.data.frame(x)
}

SS_curve_exp = SS_curve_exp[SS_curve_exp!=Inf]
SS_good_exp = SS_curve_exp[SS_curve_exp!=Inf]
SS_curve_pred = SS_curve_exp[SS_curve_exp!=Inf]
SS_good_pred = SS_curve_exp[SS_curve_exp!=Inf]

Filter = SS_curve_pred <= SS_good_pred | SS_curve_exp <= SS_good_exp

Good_curves = Filter[Filter==T]
print(paste0("Sufficiently linear CRC: ",round(length(Good_curves)/length(Filter)*100,digits=0)," %"))
###################################################################################
#Look at DRC
f = 1 #define f to whatever curve you want to investigate

Min = -3
for(i in 1:length(Data_with_CRC_list)){
  min = min(Data_with_CRC_list[[i]]$Dose)
  if(min<Min){
    Min = min
  }
}

#The minimal dose used was 10^-5 然 -> no AC10 lower than that is possible
  Curve = Curves[[f]]
  plot(Curve$Concentration,Curve$Response,type="l",col="blue",ylim = c(0,10), xlim = c(0,max(Curve$Concentration*1.1)))
  lines(Curve$Concentration,Curve$CIu,col="red")
  lines(Curve$Concentration,Curve$CIb,col="red")
  lines(Curve$Concentration,Curve$linear,col="purple")
  y = as.data.frame(Data_with_CRC_list[[f]])
  y$CI = 1.3*y$rel
  y$CI2 = 0.7*y$rel
  #max = (10^(Data_with_CRC$AC10[f]))
  Filter = is.nan(y$Dose) | is.nan(y$Response) | y$rel > 30 | y$rel <0 
  y = filter(y,Filter==F)
  model = lm(y$rel~0+exps(y$Dose))
  model2 = lm(y$CI~0+exps(y$Dose))
  model3 = lm(y$CI2~0+exps(y$Dose))
  abline(model,col="green")
  abline(model2,col="orange")
  abline(model3,col="orange")

