##################################################################
#               JRP data - ITC
#               Hieu, Masoomeh and Jaap August 2018
#               Estimation of TTC (target trajectory of CFI)
#               Fit a quadratic function for CFI data
#
##################################################################

  rm(list=ls(all=TRUE)) #To remove the hisory
  # dev.off() #To close all graphs
################################################################

# Set working directory
  setwd("C:/Users/Kevin Le/PycharmProjects/Pig Data Black Box")

# Data Processing process
  source("Package/Step0_DataTreat.R")
  # source("Package/Step0_Graphs.R")
  source("Package/Step1.1_MissData_Exceptions.R")
  source("Package/Step1.2_MissData_Estimation.R")
  source("Package/Step2.0_TTC_Filtration.R")
  source("Package/Step2.1_TTC_Q-LM.R")
  source("Package/Step2.2_TTC_QDR.R")
  source("Package/Step2.3_TTC_LM.R")
  source("Package/Step2.4_TTC_Graphs.R")
  source("Package/Step3.1_Detec_Pertub.R")