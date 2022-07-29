library(stringr)

setwd(str_remove(dirname(rstudioapi::getActiveDocumentContext()$path),"/Code"))

#This is the original data file which is included in httk from Wambaugh et al. 
chem.physical_and_invitro.data <- read.csv("chem.physical_and_invitro_data.csv")

#Those are the respective values of structural identifiers, molecular weight, plasma fraction unbound, intrinsic clearance, logKoa, pKa
Import_data <- read.csv("SI_Tables/Table_B2.csv")

#Import data into httk - necessary since we use parallel processing and add_chemtable() function in httk is not changing the data table in the package
#A: Append chem.physical_and_invitro.data for missing chemicals
to_be_added = Import_data[!(Import_data$CAS %in% chem.physical_and_invitro.data$CAS),]

if(length(to_be_added$CAS)>0){
  add_matrix <- matrix(ncol = length(chem.physical_and_invitro.data),nrow = length(to_be_added$CAS))
  colnames(add_matrix) = colnames(chem.physical_and_invitro.data)
  k <- length(chem.physical_and_invitro.data$CAS)+1
  chem.physical_and_invitro.data <- rbind(chem.physical_and_invitro.data, add_matrix)
  k1 <- k
  for (i in 1:length(to_be_added$CAS)){
    chem.physical_and_invitro.data$CAS[k1] = to_be_added$CAS[i]
    k1 <- k1+1
  }
}
chem.physical_and_invitro.data$Human.Funbound.plasma <- as.numeric(chem.physical_and_invitro.data$Human.Funbound.plasma)
chem.physical_and_invitro.data$Human.Clint <- as.numeric(chem.physical_and_invitro.data$Human.Clint)

#B: Assemble data
for(i in 1:length(Import_data$CAS)){
  z = match(Import_data$CAS[i],chem.physical_and_invitro.data$CAS)
  chem.physical_and_invitro.data$Human.Clint[z] = as.numeric(Import_data$Clint[i])
  chem.physical_and_invitro.data$Human.Funbound.plasma[z] = as.numeric(Import_data$Fub[i])
  chem.physical_and_invitro.data$logP[z] = as.numeric(Import_data$logKoa[i])
  chem.physical_and_invitro.data$MW[z] = as.numeric(Import_data$MW[i])
  chem.physical_and_invitro.data$Human.Clint.Reference[z] = Import_data$Clint_Reference[i]
  chem.physical_and_invitro.data$Human.Funbound.plasma.Reference[z] = Import_data$Fub_Reference[i]
  chem.physical_and_invitro.data$Compound[z] = Import_data$Compound[i]
  if(is.na(chem.physical_and_invitro.data$pKa_Accept.Reference[z])){
    chem.physical_and_invitro.data$pKa_Accept[z] = as.numeric(Import_data$pKa_Acceptor[i])
    chem.physical_and_invitro.data$pKa_Accept.Reference[z] = Import_data$pKa_Acceptor_Reference [i]
  }
  if(is.na(chem.physical_and_invitro.data$pKa_Donor.Reference[z])){
    chem.physical_and_invitro.data$pKa_Donor[z] = as.numeric(Import_data$pKa_Donor[i])
    chem.physical_and_invitro.data$pKa_Donor.Reference[z] = Import_data$pKa_Donor_Reference[i]
  }
}

write.csv(chem.physical_and_invitro.data,"chem.physical_and_invitro.data_for_httk_configuration.csv")
rm(list=ls())
load("httk/data/Tables.RData")
chem.physical_and_invitro.data = read.csv("chem.physical_and_invitro.data_for_httk_configuration.csv")
save.image("httk/data/Tables.RData")
devtools::install("httk")
