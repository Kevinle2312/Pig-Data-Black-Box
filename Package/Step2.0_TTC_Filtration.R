##################################################################              
#               JRP data - ITC
#               Hieu, Masoomeh and Jaap August 2018
#               Estimation of TTC (target trajectory of CFI)
#               Classify the functions used for TTC of CFI
#               The main function is quadratic-linear for CFI, linear-plateau for DFI
#               If the slope of DFI curve decreases -> linear for CFI
#               If the inflection point is too closed to last data point -> quadratic for CFI 
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
  library(ggplot2)
################################################################

# Set working directory
  setwd("C:/Users/Kevin Le/PycharmProjects/Pig Data Black Box")
  
#load dataset
  load("Data/JRPData.Rdata") #load dataset created in MissingData.step
  
  source("Package/abcd.R")
  
  #===============================================================
  # DATA PREPARATION
  #===============================================================
  
  #Create a new dataframe which will store Data after ITC estimation
  
  #Dataframe contains ITC parameters
  No.NA.Data.1 <- NULL
  
  #Order number of Animal_ID
  ID <- unique(as.factor(No.NA.Data$ANIMAL_ID))

  #===============================================================
  # For loop for automatically estimating ITC of all pigs
  #===============================================================
  IDC <- seq_along(ID)
  for (idc in IDC){
    # idc =100
    i <- ID[idc]
    Data <- No.NA.Data[No.NA.Data$ANIMAL_ID == i,]
    
    ####### Create data frame of x (Age) and y (CFI) ########
    x <- as.numeric(Data$Age.plot)
    Y <- as.numeric(Data$CFI.plot)
    Z <- as.numeric(Data$DFI.plot)
    Data.xy <- as.data.frame(cbind(x,Y))
    
    #Initial parameteres for parameter estimation
    X0.0 <- x[1]
    Xlast <- x[length(x)]
    
    #Provide set of initial parameters
    Xs.1 <- round(seq(X0.0 + 1, Xlast - 1, len = 30), digits = 0)
    X0.1 <- rep(X0.0, length(Xs.1))
    DFIs.1 <- NULL
    CFIs.1 <- NULL
    for(A in seq_along(Xs.1)){
      DFIs2 <- Data[Data$Age.plot == Xs.1[A],]$DFI.plot
      CFIs2 <- Data[Data$Age.plot == Xs.1[A],]$CFI.plot
      DFIs.1 <- c(DFIs.1, DFIs2)
      CFIs.1 <- c(CFIs.1, CFIs2)
    }
    
    st1 <- data.frame(cbind(X0.1,
                            Xs.1,
                            DFIs.1,
                            CFIs.1))
    names(st1) <- c("X0","Xs", "DFIs","CFIs")

    #RUN NLS2 with upper and lower bound
    # weight = 1/Y^2
    st2 <- nls2(Y ~ nls.func.2(X0, Xs, DFIs, CFIs),
                Data.xy,
                start = st1,
                # weights = weight,
                # trace = T,
               algorithm = "brute-force")
    par_init <- coef(st2)
    par_init
    
    c <- abcd.2(par_init)[3];c
    Xs.0 <- par_init[2]; Xs.0
    
    if(c < 0){
      FuncType <- "LM"
    } else if(c > 0 &&
              Xs.0 >= ((Xlast - X0.0)*0.2 + X0.0)
              # && Xs.0 <= ((Xlast - X0.0)*1 + X0.0)
              ){ 
      FuncType <- "QLM"
    } else {
      FuncType <- "QDR"
    } 
    #
    #Plot Function with initial parameters to data
    ITC_CFI.0 <- as.numeric(unlist(pred.func.2(par_init, x)[1]))
    ITC_DFI.0 <- as.numeric(unlist(pred.func.2(par_init, x)[2]))

    #Plot
    # png(filename = paste("Graphs/Step2_graphs/Initial_para/",idc,".",ID[idc],".png",sep=""),
    #     height = 720, width = 1200, units = 'px', type="cairo-png")
    par(mfrow = c(1,2))
    par(mar=c(4.5,4.5,4.5,1.5))
    #CFI
    plot(x, ITC_CFI.0,
         main = "Cumulative Feed Intake",
         ylab = "Cumulative Feed Intake (kg)",
         xlab = "Age (days)",
         ylim = c(0, max(ITC_CFI.0, Data$CFI.plot)),
         type="l", col ="blue", cex.main = 1.5,
         cex.axis = 1.5, cex.lab = 1.5, lwd = 2.7) # ITC of CFI
    points(x, Data$CFI.plot,
           col="black", cex = 1.5) # Observation data of CFI
    # DFI
    plot(x, Data$DFI.plot, # Observation data of DFI
         main = "Daily Feed Intake",
         ylab = "Daily Feed Intake (kg)",
         xlab = "Age (days)",
         type="p", col ="black", cex.main = 1.5, cex = 1.5,
         cex.axis = 1.5, cex.lab = 1.5)
    lines(x, ITC_DFI.0, # ITC of DFI
          type = "l",col="blue", lwd = 2.7)
    mtext(paste("Pig:", ID[idc], ", idc=", idc, ", FuncType:", FuncType), side = 3, line = -1.8, outer = TRUE, cex = 2)
    # par(mfrow = c(1,1))
    dev.off()
    
    #--------------------------------------------
    # Choose function for the TTC
    #--------------------------------------------
    
    idc.1 <- rep(idc, dim(Data)[1])
    FuncType <- rep(FuncType, dim(Data)[1])
    Data$FuncType <- FuncType
    
    Data$idc.1 <- idc.1
    
    
    No.NA.Data.1 <- rbind(No.NA.Data.1 , Data)
    
  } # end FOR loop
  
  # No.NA.Data <- No.NA.Data.1 
  
  DFI.neg <- subset(No.NA.Data.1, FuncType == "LM")
  length(unique(DFI.neg$idc.1))
  DFI.pos1 <- subset(No.NA.Data.1, FuncType == "QDR")
  length(unique(DFI.pos1$idc.1))
  DFI.pos2 <- subset(No.NA.Data.1, FuncType == "QLM")
  length(unique(DFI.pos2$idc.1))
  
  DFI.neg <- unique(DFI.neg$ANIMAL_ID); DFI.neg
  DFI.pos1 <- unique(DFI.pos1$ANIMAL_ID); DFI.pos1
  DFI.pos2 <- unique(DFI.pos2$ANIMAL_ID); DFI.pos2
  DFI.pos1 <- as.numeric(as.character(DFI.pos1))
  DFI.pos2 <- as.numeric(as.character(DFI.pos2))
  DFI.neg <- as.numeric(as.character(DFI.neg))
  #===============================================================
  #Save results to Rdata file
  #===============================================================

  save(No.NA.Data.1, DFI.neg, DFI.pos1, DFI.pos2, 
       file = "Data/JRPData_TTC.RData")

  # write.csv(No.NA.Data.1, DFI.neg, DFI.pos1, DFI.pos2, file = "JRP.param.csv")
