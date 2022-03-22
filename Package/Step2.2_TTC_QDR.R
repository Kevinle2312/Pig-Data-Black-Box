##################################################################              
#               JRP data - ITC
#               Hieu, Masoomeh and Jaap August 2018
#               Estimation of TTC (target trajectory of CFI)
#               Fit a quadratic function for CFI data
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
  
################################################################

# Set working directory
  setwd("C:/Users/Kevin Le/PycharmProjects/Pig Data Black Box")
  
#load dataset
  load("Data/JRPData_TTC.Rdata") #load dataset created in MissingData.step
  
  source("Package/Functions.R")
  
  #===============================================================
  # DATA PREPARATION
  #===============================================================
  ID <- DFI.pos1; ID
  length(ID)
  #Create a new dataframe which will store Data after ITC estimation
  
  #Dataframe contains ITC parameters
  ITC.param.pos1 <- data.frame(ANIMAL_ID=factor(),
                               X0=double(),
                               Y2=double(),
                               Ylast=double(),
                               a=double(),
                               b=double(),
                               c=double(),
                               stringsAsFactors=FALSE)
  
  #Dataframe contains data points on the ITC
  Data.remain <- data.frame(ANIMAL_ID=character(),
                            Age=double(),
                            obs.CFI=double(),
                            tt=double(),
                            ttt=double(),
                            stringsAsFactors=FALSE)

  #===============================================================
  # For loop for automatically estimating ITC of all pigs
  #===============================================================
  IDC <- seq_along(ID) # 67 JRP
  for (idc in IDC){
    # idc = 1
    i <- ID[idc]
    Data <- No.NA.Data.1[No.NA.Data.1$ANIMAL_ID == i,]
    idc1 <- unique(as.numeric(Data$idc.1))
    ####### Create data frame of x (Age) and y (CFI) ########
    x <- as.numeric(Data$Age.plot)
    Y <- as.numeric(Data$CFI.plot)
    Z <- as.numeric(Data$DFI.plot)
    Data.xy <- as.data.frame(cbind(x,Y))
    
    #Initial parameteres for parameter estimation
    xlast <- x[length(x)]
    x0.0 <- x[1]
    y2.0 <- Y[floor(length(x)*2/3)]
    ylast.0 <- Y[length(x)]
    #Vector contains 4 initial parameters 
    par_init <- c(x0.0, y2.0, ylast.0)  
    
    x2 <- 2*(xlast - x[1])/3+x[1]
    ##-------test the prediction by initial parametere values -------------

    ##################################################################
    # 1. reparametrization CFI at X0 = 0
    #function used for reparametrization in MAPLE 
    # solve({0=a+b*X0+c*X0**2+d*X0**3,
    # y1 = a + b*(X0+1/3*(xlast-X0)) +
    #   c*(X0+1/3*(xlast-X0))**2 + d*(X0+1/3*(xlast-X0))**3,
    # y2 = a + b*(X0+2/3*(xlast-X0)) + c*(X0+2/3*(xlast-X0))**2 +
    #   d*(X0+2/3*(xlast-X0))**3,
    # ylast = a+b*xlast+c*xlast**2+d*xlast**3},{a,b,c,d});
    #  2. with the source of the function abcd and pred
    ##################################################################    
   
    #--------------------------------------------
    # Create empty lists to store data after loop
    #--------------------------------------------
    
    par <- list()
    AC.res <- list()
    AC.pvalue <- NULL
    data2 <- list()
    data3 <- list()
    param <- data.frame(rbind(par_init))
    par.abcd <- data.frame(rbind(abcd.1(as.vector(par_init))))
    param.2 <- data.frame(X0=double(), 
                          Y2=double(),
                          Ylast=double(),
                          a=double(),
                          b=double(),
                          c=double(),
                          stringsAsFactors=FALSE) 
    j <- 2
    AC_pvalue <- 0
    AC.pvalue[1] <- AC_pvalue
    datapointsleft <- as.numeric(dim(Data)[1])
    dpl <- datapointsleft #vector of all dataponitsleft at each step
 
    #-------------------------------------------------------------------------------
    # Start the procedure of Non Linear Regression
    #-------------------------------------------------------------------------------
    # 
            while (AC_pvalue<=0.05 && datapointsleft >= 20){
              weight <- 1/Y^2
              #---------------- NON linear reg applied to log(Y) ---------------------------------
              
              nls.CFI <- nlsLM(Y ~ nls.func.1(X0, y2, ylast),
                             Data.xy, 
                             control = list(tol = 1e-2, printEval = TRUE, maxiter = 50),
                             start = list(X0 = par_init[1], y2 = par_init[2],
                                          ylast = par_init[3]), 
                             trace = F,
                             weights = weight)
              
              
              #--------RESULTS analysis GOODNESS of fit
              #estimate params
              par[[j]] <- coef(nls.CFI)
              par.abcd[j,] <- abcd.1(as.vector(coef(nls.CFI) )) #calculation of a, b, c and d
              param[j,] <- par[[j]]
              param.2[j-1,] <- cbind(param[j,], par.abcd[j,])
              
              #Calculation of days associated with Y1 and Y2
              x2[j] <- (xlast - coef(nls.CFI)[1])*2/3+coef(nls.CFI)[1]
              # #summary
              # summ = overview((nls.CFI))  #summary
              #residuals
              res1 <- nlsResiduals(nls.CFI) #residuals
              res2 <- nlsResiduals(nls.CFI)$resi1
              res <- res2[, 2]
              AC.res <- test.nlsResiduals(res1)
              AC.pvalue[j] <- AC.res$p.value
              
              #---------Check for negative residuals----------
              
              #Add filtration step order to data
              Step <- rep(j - 1, length(x)) 
              #create a new dataset with predicted CFI included
              Data.new <- data.frame(cbind(x, Z, Y, pred.func.1(par[[j]],x)[[1]], res, Step))
              names(Data.new) <- c("Age", "Observed_DFI","Observed_CFI", "Predicted_CFI", "Residual", "Step")
              # plot(Data.new$Age, Data.new$Predicted_CFI, type = "l", col = "black",lwd = 2,
              #      ylim = c(0, max(Data.new$Predicted_CFI, Data.new$Observed_CFI)))
              # lines(Data.new$Age, Data.new$Observed_CFI, type = "p", cex = 1.5)
              #
              #remove negative res
              Data.pos <- Data.new[!Data.new$Residual<0,]
              # lines(Data.pos$Age, Data.pos$Predicted_CFI, type = "l", col = j-1, lwd = 2)
              # lines(Data.pos$Age, Data.pos$Observed_CFI, type = "p", col = j, cex = 1.5)
              
              #restart 
              datapointsleft <- as.numeric(dim(Data.pos)[1])
              par_init <- par[[j]]
              AC_pvalue <- AC.pvalue[j]
              j <- j+1
              x <- Data.pos$Age
              Y <- Data.pos$Observed_CFI
              Z <- Data.pos$Observed_DFI
              Data.xy <- as.data.frame(cbind(x,Y))
              dpl <- c(dpl, datapointsleft); dpl
              }
    
    ANIMAL_ID <- rep(i, dim(param.2)[1])
    XLAST <- rep(xlast, dim(param.2)[1])
    FuncType <- unique(Data$FuncType) 
    param.2 <- cbind(ANIMAL_ID,
                     param.2,
                     x2[2:length(x2)],
                     XLAST,
                     AC.pvalue[2:length(AC.pvalue)],
                     dpl[seq_along(dpl) -1],
                     rep(FuncType, dim(param.2)[1]),
                     rep(idc1, dim(param.2)[1])
                    )
    colnames(param.2) <- c("ANIMAL_ID",
                           "X0",
                           "Y2",
                           "Ylast",
                           "a",
                           "b",
                           "c",
                           "X2",
                           "Xlast",
                           "P.runs.trst",
                           "DPL",
                           "FuncType",
                           "idc1"
                          )
    
    param.2  
    #Add one column about the slope of DFI to the function
    Slope <- NULL
    for(ii in 1:dim(param.2)[1]){
      
      if(param.2[ii,]$c < 0){
        Slope.1 <- -1
      } else {
        Slope.1 <- 1
      }
        
        Slope <- c(Slope, Slope.1)
      }
  
    param.2$Slope <- Slope  

 #-------------------------------------------------------------------------------
 # Give Animal ID for each animal and save them to dataframes in the loop
 #-------------------------------------------------------------------------------    
   ANIMAL_ID <- rep(i, dim(Data.new)[1])
   idc1 <- rep(idc1, dim(Data.new)[1])
   Data.new <- cbind(ANIMAL_ID , Data.new, idc1)
    
   Data.remain <- rbind(Data.remain, Data.new)
   
   Age.remain.pos1 <- Data.remain[, c(1:4, 8)]
   
   ITC.param.pos1 <- rbind(ITC.param.pos1 , param.2)
   
  } # end FOR loop
  
  #Check number of animals in this function type
  #Quadratic function for CFI
  length(unique(ITC.param.pos1$ANIMAL_ID)); length(unique(Age.remain.pos1$ANIMAL_ID))

  #==============================================================================
  # Regrouping the animals have negative slope in QDR function to Linear function
  #==============================================================================
  #Select animals have negative slope for DFI in Function type 2
  ID <- unique(ITC.param.pos1$ANIMAL_ID)
  ITC.param.p1 <- data.frame()

  for(AA in seq_along(ID)){
    i <- ID[AA]
    param.p1.neg <- ITC.param.pos1[ITC.param.pos1$ANIMAL_ID == i,]
    param.p1.neg <- param.p1.neg[dim(param.p1.neg)[1],]


    ITC.param.p1 <- rbind(ITC.param.p1, param.p1.neg)
  }
  ITC.param.p1.neg <- subset(ITC.param.p1, Slope == -1)

  #Remove animals have negative slope of DFI out of dataset for Function type 3
  ITC.param.pos1 <- subset(ITC.param.pos1, !(ANIMAL_ID %in% ITC.param.p1.neg$ANIMAL_ID))
  Age.remain.pos1 <- subset(Age.remain.pos1, !(ANIMAL_ID %in% ITC.param.p1.neg$ANIMAL_ID))

  #Add these animals to the group of linear function for CFI
  DFI.neg1 <- c( DFI.neg, unique(ITC.param.p1.neg$ANIMAL_ID))
  DFI.neg <- unique(DFI.neg1)
  # DFI.neg <- as.numeric(as.character(DFI.neg))
  length(DFI.neg)

  #Reset the animals in group of Quadratic function for CFI
  DFI.pos1 <- unique(Age.remain.pos1$ANIMAL_ID)

  # Reset all FuncType to "QDR"
  ITC.param.pos1$FuncType <- "QDR"

  #===============================================================
  #Save results to Rdata file
  #===============================================================

  save(ITC.param.pos1, Age.remain.pos1, 
       file = "Data/JRP.DFI.pos1.RData")

  save(No.NA.Data.1, DFI.neg, DFI.pos1, DFI.pos2,
       file = "Data/JRPData_TTC.RData")

#DFI.neg <- factor(DFI.neg1, levels = 1:nlevels(DFI.pos2), labels = levels(DFI.pos2))