##################################################################              
#               JRP data - ITC
#               Hieu, Masoomeh and Jaap August 2018
#               Estimation of TTC (target trajectory of CFI)
#               Draw graph for TTC curves
#
##################################################################

rm(list=ls(all=TRUE)) #To remove the hisory 
# dev.off() #To close all graphs

#-------------------------------------------------------------------------------
# Packages needed for estimaton of Ideal trajectory - nonlinear regression
#-------------------------------------------------------------------------------
library("minpack.lm")
library("nlstools")
library("nlsMicrobio")
library("stats") 
library("tseries") #runs test for auto correlation
#Use NLS2
library(proto)
library(nls2)

################################################################

# Set working directory
setwd("C:/Users/Kevin Le/PycharmProjects/Pig Data Black Box")

#load dataset
source("Package/abcd.R")
load("Data/JRPData_TTC.Rdata")
load("Data/JRP.DFI.neg.Rdata")
load("Data/JRP.DFI.pos1.Rdata")
load("Data/JRP.DFI.pos2.Rdata")

#Check number of animals in each function type
#LM function for CFI
length(unique(ITC.param.neg$ANIMAL_ID)); length(unique(Age.remain.neg$ANIMAL_ID))
#Quadratic function for CFI
length(unique(ITC.param.pos1$ANIMAL_ID)); length(unique(Age.remain.pos1$ANIMAL_ID))
#Quadratic-linear function for CFI
length(unique(ITC.param.pos2$ANIMAL_ID)); length(unique(Age.remain.pos2$ANIMAL_ID))

#Order number of Animal_ID
ID <- unique(as.factor(No.NA.Data.1$ANIMAL_ID))

#===============================================================
# For loop for automatically estimating ITC of all pigs
#===============================================================
IDC <- seq_along(ID)
for (idc in IDC){
  # idc = 3
  i <- ID[idc]
  Data <- No.NA.Data.1[No.NA.Data.1$ANIMAL_ID == i,]
  idc1 <- unique(as.numeric(Data$idc.1))
  AGE <- Data$Age.plot # initial Ages of Animal i
  CFI.exp <- Data$CFI.plot # initial CFI of animal i
  DFI.exp <- Data$DFI.plot # initial DFI of animal i
  
  
  #Choose the right function for the right animal
  
  if(i %in% unique(ITC.param.neg$ANIMAL_ID)){
    param.i <- ITC.param.neg[ITC.param.neg$ANIMAL_ID == i,]
    FuncType <- param.i[dim(param.i)[1],]$FuncType
    Slope <- param.i[dim(param.i)[1],]$Slope
    param.i <- as.numeric(param.i[dim(param.i)[1], 4:5])
    ITC <- pred.abcd.0(param.i, AGE)[[1]]
    ITD <- rep(param.i[2], length(AGE))
    
    Data.remain <- Age.remain.neg[Age.remain.neg$ANIMAL_ID == i,]
    Age.remain <- Data.remain$Age
    DFI.remain <- Data.remain$Observed_DFI
    CFI.remain <- Data.remain$Observed_CFI
    
      
  } else if(i %in% unique(ITC.param.pos1$ANIMAL_ID)){
    param.i <- ITC.param.pos1[ITC.param.pos1$ANIMAL_ID == i,]
    FuncType <- param.i[dim(param.i)[1],]$FuncType
    Slope <- param.i[dim(param.i)[1],]$Slope
    param.i <- as.numeric(param.i[dim(param.i)[1], 5:7])
    ITC <- pred.abcd.1(param.i, AGE)[[1]]
    ITD <- pred.abcd.1(param.i, AGE)[[2]]
    
    Data.remain <- Age.remain.pos1[Age.remain.pos1$ANIMAL_ID == i,]
    Age.remain <- Data.remain$Age
    DFI.remain <- Data.remain$Observed_DFI
    CFI.remain <- Data.remain$Observed_CFI
    
  } else{
    param.i <- ITC.param.pos2[ITC.param.pos2$ANIMAL_ID == i,]
    FuncType <- param.i[dim(param.i)[1],]$FuncType
    Slope <- param.i[dim(param.i)[1],]$Slope
    Xs <- param.i[dim(param.i)[1],]$Xs
    param.i <- as.numeric(param.i[dim(param.i)[1], 6:8])
    ITC <- pred.abcd.2(param.i, AGE)[[1]]
    ITD <- pred.abcd.2(param.i, AGE)[[2]]
    
    Data.remain <- Age.remain.pos2[Age.remain.pos2$ANIMAL_ID == i,]
    Age.remain <- Data.remain$Age
    DFI.remain <- Data.remain$Observed_DFI
    CFI.remain <- Data.remain$Observed_CFI
    
  }
  
  #------------------------------------------------------
  #     Plot ITC of CFI 
  #------------------------------------------------------
  png(filename = paste0("Graphs/Step2_graphs/Final_graphs/", idc, ".", ID[idc], ".png"),
      height = 720, width = 1200, units = 'px', type="cairo-png")
  par(mfrow = c(1,2))
  par(mar=c(4.5,4.5,4.5,1.5))
  #CFI
  plot(AGE,
       ITC,
       main = "Cumulative Feed Intake",
       ylab = "Cumulative Feed Intake (kg)",
       xlab = "Age (days)",
       type="l", col ="blue", cex.main = 1.5,
       cex.axis = 1.5, cex.lab = 1.5, lwd = 2.7) # ITC of CFI
  points(AGE,
         CFI.exp,
         col="black", cex = 1.5) # Observation data of CFI
  points(Age.remain,
         CFI.remain,
         col="blue", pch = 16, cex = 1.3) # Datapoints left of CFI belong to ITC
  legend("topleft", c("Real data", "data points of ITC", "ITC"), 
         col=c("black", "blue", "blue"), 
         pch=c(1, 16, NA), lty = c(NA, NA,1), bty = "n")
  # DFI
  plot(AGE,
       DFI.exp, # Observation data of DFI
       main = "Daily Feed Intake",
       ylab = "Daily Feed Intake (kg)",
       xlab = "Age (days)",
       # ylim = c(min(ITC_DFI, Data$DFI.plot), max(ITC_DFI, Data$DFI.plot)),
       type="p", col ="black", cex.main = 1.5, cex = 1.5, 
       cex.axis = 1.5, cex.lab = 1.5) 
  lines(AGE,
        ITD, # ITC of DFI
        type = "l",col="blue", lwd = 2.7) 
  points(Age.remain,
         DFI.remain,
         col="blue", pch = 16, cex = 1.5) # Datapoints left of DFI belong to ITC
  mtext(paste("Pig ID:", ID[idc], ", FuncType:", FuncType, ", Slope =", Slope), side = 3, line = -1.8, outer = TRUE, cex = 2)
  par(mfrow = c(1,1))
  dev.off()
  
} # end FOR loop
