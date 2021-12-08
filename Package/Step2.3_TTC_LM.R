##################################################################              
#               JRP data - ITC
#               Hieu, Masoomeh and Jaap August 2018
#               Estimation of TTC (target trajectory of CFI)
#               Fit a linear function for CFI data
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
  
  source("Package/abcd.R")
  
  #===============================================================
  # DATA PREPARATION
  #===============================================================
  
  #Order number of Animal_ID
  ID <- DFI.neg; ID
  length(ID)
  #Create a new dataframe which will store Data after ITC estimation
  
  #Dataframe contains ITC parameters
  ITC.param.neg <- data.frame(ANIMAL_ID=factor(),
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
  IDC <- seq_along(ID)
  for (idc in IDC){
    # idc = 3
    i <- ID[idc]
    Data <- No.NA.Data.1[No.NA.Data.1$ANIMAL_ID == i,]
    idc1 <- unique(as.numeric(Data$idc.1))
    ####### Create data frame of x (Age) and y (CFI) ########
    x <- as.numeric(Data$Age.plot)
    Y <- as.numeric(Data$CFI.plot)
    Z <- as.numeric(Data$DFI.plot)
    Data.xy <- as.data.frame(cbind(x,Y))
    
    #Initial parameteres for parameter estimation
    x0.0 <- x[1]
    xlast <- x[length(x)]
    ylast.0 <- Y[length(x)]
    
    #Vector contains 4 initial parameters 
    par_init <- c(x0.0, ylast.0)  
    
    ##-------test the prediction by initial parametere values -------------

    ##################################################################
    # 1. reparametrization CFI at X0 = 0
    #function used for reparametrization in MAPLE 
    # solve({0 = X0*b+a, Ylast = Xlast*b+a}, {a, b});
    # a = (Ylast*X0)/(X0 - Xlast), b = -(Ylast)/(X0-Xlast)
    #  2. with the source of the function abcd and pred
    ##################################################################    
   
    #
    #-------------------------------------------------------
    # Function to reparameter a, b, c and of cubic function:
    #-------------------------------------------------------
    
    #--------------------------------------------
    # Create empty lists to store data after loop
    #--------------------------------------------
    
    par <- list()
    AC.res <- list()
    AC.pvalue <- NULL
    data2 <- list()
    param <- data.frame(rbind(par_init))
    par.abcd <- data.frame(rbind(abcd.0(as.vector(par_init))))
    param.2 <- data.frame(X0=double(), 
                          Ylast=double(),
                          a=double(),
                          b=double(),
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
              
              nls.CFI <- nlsLM(Y ~ nls.func.0(X0, ylast),
                             Data.xy, 
                             control = list(tol = 1e-2, printEval = TRUE, maxiter = 50),
                             start = list(X0 = par_init[1],
                                          ylast = par_init[2]), 
                             trace = F,
                             weights = weight)
              
              
              #--------RESULTS analysis GOODNESS of fit
              #estimate params
              par[[j]] <- coef(nls.CFI)
              par.abcd[j,] <- abcd.0(as.vector(coef(nls.CFI) )) #calculation of a, b, c and d
              param[j,] <- par[[j]]
              param.2[j-1,] <- cbind(param[j,], par.abcd[j,])
              
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
              Data.new <- data.frame(cbind(x, Z, Y, pred.func.0(par[[j]],x)[[1]], res, Step))
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
                     XLAST,
                     AC.pvalue[2:length(AC.pvalue)],
                     dpl[seq_along(dpl) -1],
                     rep(FuncType, dim(param.2)[1]),
                     rep(idc1, dim(param.2)[1])
                    )
    colnames(param.2) <- c("ANIMAL_ID",
                           "X0",
                           "Ylast",
                           "a",
                           "b",
                           "Xlast",
                           "P.runs.trst",
                           "DPL",
                           "FuncType",
                           "idc1"
                         )
    
    param.2  

    #Add one column about the slope of DFI to the function
    Slope <- rep(0, dim(param.2)[1])
    
    param.2$Slope <- Slope  
 
    #-------------------------------------------------------------------------------
    # Give Animal ID for each animal and save them to dataframes in the loop
    #-------------------------------------------------------------------------------    
    ANIMAL_ID <- rep(i, dim(Data.new)[1])
    idc1 <- rep(idc1, dim(Data.new)[1])
    Data.new <- cbind(ANIMAL_ID , Data.new, idc1)
    
    Data.remain <- rbind(Data.remain, Data.new)
    Age.remain.neg <- Data.remain[, c(1:4, 8)]
    ITC.param.neg <- rbind(ITC.param.neg , param.2 )
    
  } # end FOR loop
  
#Check number of animals in this function type
#LM function for CFI
length(unique(ITC.param.neg$ANIMAL_ID)); length(unique(Age.remain.neg$ANIMAL_ID))
  
  
  #===============================================================
  #Save results to Rdata file
  #===============================================================

  save(ITC.param.neg, Age.remain.neg,         #Step 2: ITC
       file = "Data/JRP.DFI.neg.RData")
  
