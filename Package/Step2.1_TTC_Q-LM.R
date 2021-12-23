##################################################################              
#               JRP data - ITC
#               Hieu, Masoomeh and Jaap August 2018
#               Estimation of TTC (target trajectory of CFI)
#               Fit a quadratic-linear function for CFI data
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
  load("Data/JRPData_TTC.Rdata") #load dataset created in MissingData.step
  
  source("Package/Functions.R")
  
  #===============================================================
  # DATA PREPARATION
  #===============================================================
  # DFI.pos2 <- subset(No.NA.Data.1, FuncType == "QLM")
  # DFI.pos2 <- unique(DFI.pos2$ANIMAL_ID); length(DFI.pos2)
  ID <- DFI.pos2; ID; length(ID)
  #Create a new dataframe which will store Data after ITC estimation
  
  #Dataframe contains ITC parameters
  ITC.param.pos2 <- data.frame(ANIMAL_ID=factor(),
                               X0=double(),
                               Y1=double(),
                               Y2=double(),
                               Ylast=double(),
                               a=double(),
                               b=double(),
                               c=double(),
                               d=double(),
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
  IDC <- c(1:31, 32:length(ID)) # 17, 23, 52, 57, 116
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
    X0.0 <- x[1]
    Xlast <- x[length(x)]
    
    ##################################################################
    # 1. reparametrization CFI at X0 = 0
    #function used for reparametrization in MAPLE 
    # solve({
    # 0=a+b*X_0+c*X_0**2,
    # DFIs=b+2*c*Xs,CFIs=a+b*Xs+c*Xs**2},
    # {a,b,c});
    # a = -X0*(2*CFIs*Xs-CFIs*X0-Xs^2*DFIs+Xs*DFIs*X0)/(Xs^2-2*X0*Xs+X0^2)
    # b = (-Xs^2*DFIs+DFIs*X0^2+2*CFIs*Xs)/(Xs^2-2*X0*Xs+X0^2)
    # c = -(CFIs-Xs*DFIs+X0*DFIs)/(Xs^2-2*X0*Xs+X0^2)
    
    #  2. with the source of the function abcd and pred
    ##################################################################    
   
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
    
    st1 <- data.frame(cbind(X0.1, Xs.1, DFIs.1, CFIs.1))
    names(st1) <- c("X0","Xs", "DFIs","CFIs")

    #RUN NLS2 to find optimal initial parameters
    st2 <- nls2(Y ~ nls.func.2(X0, Xs, DFIs, CFIs),
                Data.xy,
                start = st1,
                # weights = weight,
                # trace = T,
               algorithm = "brute-force")
    
    par_init <- coef(st2); par_init
    
    #--------------------------------------------
    # Create empty lists to store data after loop
    #--------------------------------------------
    
    par <- list()
    AC.res <- list()
    AC.pvalue <- NULL
    data2 <- list()
    data3 <- list()
    param <- data.frame(rbind(par_init))
    par.abcd <- data.frame(rbind(abcd.2(as.vector(par_init))))
    param.2 <- data.frame(X0=double(), 
                          Xs=double(),
                          DFIs=double(),
                          CFIs=double(),
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
    
            while ((AC_pvalue<=0.05) && datapointsleft >= 20){
              weight <- 1/Y^2
              # ---------------- NON linear reg applied to log(Y) ---------------------------------
              st2 <- nls2(Y ~ nls.func.2(X0, Xs, DFIs, CFIs),
                          Data.xy,
                          start = st1,
                          weights = weight,
                          trace = F,
                          algorithm = "brute-force")
              par_init <- coef(st2)
              par_init
              st1 <- st1[!(st1$Xs == par_init[2]),]
              # nls.CFI <- nlsLM(Y ~ nls.func.2(X0, Xs, DFIs, CFIs),
              #                Data.xy,
              #                control = list(tol = 1e-2, printEval = TRUE, maxiter = 1024),
              #                start = list(X0 = par_init[1], Xs = par_init[2],
              #                             DFIs = par_init[3], CFIs = par_init[4]),
              #                weights = weight,
              #                algorithm = "port",
              #                lower = c(-100000,X0.0+1, -100000, -100000),
              #                upper = c(100000, Xlast-1, 100000, 100000),
              #                trace = F)

              nls.CFI <- nls2(Y ~ nls.func.2(X0, Xs, DFIs, CFIs),
                                              Data.xy,
                                              start = st1,
                                              weights = weight,
                                              control = nls.control(warnOnly = TRUE),
                                              trace = T,
                                              algorithm = "port",
                                              lower = c(-100000,X0.0+1, -100000, -100000),
                                              upper = c(100000, Xlast-1, 100000, 100000))
              
              #--------RESULTS analysis GOODNESS of fit
              #estimate params
              par[[j]] <- coef(nls.CFI)
              par.abcd[j,] <- abcd.2(as.vector(coef(nls.CFI) )) #calculation of a, b, c and d
              param[j,] <- par[[j]]
              param.2[j-1,] <- cbind(param[j,], par.abcd[j,])
              
              #summary
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
              Data.new <- data.frame(cbind(x, Z, Y, pred.func.2(par[[j]],x)[[1]], res, Step))
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
              
              #Criteria to stop the loop when the estimated parameters are equal to initial parameters
              # Crite <- sum(param.2[dim(param.2)[1],c(1:4)] == par_init) 
              
              datapointsleft <- as.numeric(dim(Data.pos)[1])
              par_init <- par[[j]]
              AC_pvalue <- AC.pvalue[j]
              j <- j+1
              x <- Data.pos$Age
              Y <- Data.pos$Observed_CFI
              Z <- Data.pos$Observed_DFI
              Data.xy <- as.data.frame(cbind(x,Y))
              dpl <- c(dpl, datapointsleft)
              dpl
              #Create again the grid
              X0.0 <- x[1]
              Xlast <- x[length(x)]
              #Xs
              if(par_init[2] -15 <= X0.0){
                Xs.1 <- round(seq(X0.0 + 5, Xlast - 5, len = 30), digits = 0)
              } else if(par_init[2] + 5 >= Xlast){
                Xs.1 <- round(seq(par_init[2]-10, par_init[2]-1, len = 6), digits = 0)
              } else{
                Xs.1 <- round(seq(par_init[2]-5, par_init[2] + 5, len = 6), digits = 0)
              }
              #
              X0.1 <- rep(X0.0, length(Xs.1))
              DFIs.1 <- NULL
              CFIs.1 <- NULL
              for(A in seq_along(Xs.1)){
                DFIs2 <- Data[Data$Age.plot == Xs.1[A],]$DFI.plot
                CFIs2 <- Data[Data$Age.plot == Xs.1[A],]$CFI.plot
                DFIs.1 <- c(DFIs.1, DFIs2)
                CFIs.1 <- c(CFIs.1, CFIs2)
              }
              st1 <- data.frame(cbind(X0.1, Xs.1, DFIs.1, CFIs.1))

              if(X0.0 <= par_init[2] && Xlast >=par_init[2]){
              st1 <- rbind(st1, par_init)
              }
              names(st1) <- c("X0","Xs", "DFIs","CFIs")
              }
    
      ANIMAL_ID <- rep(i, dim(param.2)[1])
      FuncType <- unique(Data$FuncType) 
      
      param.2 <- cbind(ANIMAL_ID,
                       param.2,
                       XLAST = rep( Data$Age.plot[length(Data$Age.plot)], dim(param.2)[1]),
                       AC.pvalue[2:length(AC.pvalue)],
                       dpl[seq_along(dpl) -1],
                       rep(FuncType, dim(param.2)[1]),
                       rep(idc1, dim(param.2)[1])
                      )
      colnames(param.2) <- c("ANIMAL_ID",
                             "X0",
                             "Xs",
                             "DFIs",
                             "CFIs",
                             "a",
                             "b",
                             "c",
                             "Xlast",
                             "P.runs.trst",
                             "DPL",
                             "FuncType",
                             "idc1"
                           )
      
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
      
      param.2$Slope <- Slope ; param.2 
  
 #-------------------------------------------------------------------------------
 # Give Animal ID for each animal and save them to dataframes in the loop
 #-------------------------------------------------------------------------------    
      ANIMAL_ID <- rep(i, dim(Data.new)[1])
      idc1 <- rep(idc1, dim(Data.new)[1])
      Data.new <- cbind(ANIMAL_ID , Data.new, idc1)
      
      Data.remain <- rbind(Data.remain, Data.new)
      
      Age.remain.pos2 <- Data.remain[, c(1:4, 8)]
      
      ITC.param.pos2 <- rbind(ITC.param.pos2 , param.2 )
  
  } # end FOR loop
  
#Check number of animals in this function type  
#Quadratic-linear function for CFI
length(unique(ITC.param.pos2$ANIMAL_ID)); length(unique(Age.remain.pos2$ANIMAL_ID))

#==============================================================================
# Regrouping the animals have negative slope in Q-L function to Linear function
#============================================================================== 
#Select animals have negative slope for DFI in Q-L Function for CFI 
ID <- unique(ITC.param.pos2$ANIMAL_ID)
ITC.param.p2 <- data.frame()

for(AA in seq_along(ID)){
  i <- ID[AA]
  param.p2.neg <- ITC.param.pos2[ITC.param.pos2$ANIMAL_ID == i,]
  param.p2.neg <- param.p2.neg[dim(param.p2.neg)[1],]
  
  
  ITC.param.p2 <- rbind(ITC.param.p2, param.p2.neg)
}
ITC.param.p2.neg <- subset(ITC.param.p2, Slope == -1)
length(unique(ITC.param.p2.neg$ANIMAL_ID))

#Remove animals have negative slope of DFI out of dataset for Function of Linear-Plateau
ITC.param.pos2 <- subset(ITC.param.pos2, !(ANIMAL_ID %in% ITC.param.p2.neg$ANIMAL_ID))
Age.remain.pos2 <- subset(Age.remain.pos2, !(ANIMAL_ID %in% ITC.param.p2.neg$ANIMAL_ID))

#Add these animals to the group of linear function for CFI
DFI.neg1 <- c( DFI.neg, unique(ITC.param.p2.neg$ANIMAL_ID))
DFI.neg <- unique(DFI.neg1)
# DFI.neg <- as.numeric(as.character(DFI.neg))
length(DFI.neg)

#==============================================================================
# Regrouping the animals have Xs too closed to Xlast in Q-L function to QDR function
#==============================================================================

#Select animals have XS too closed to Xlast Function type 3
ID <- unique(ITC.param.pos2$ANIMAL_ID)

ITC.param.p2 <- data.frame()

for(AA in seq_along(ID)){
  i <- ID[AA]
  param.p2.neg <- ITC.param.pos2[ITC.param.pos2$ANIMAL_ID == i,]
  param.p2.neg <- param.p2.neg[dim(param.p2.neg)[1],]
  
  
  ITC.param.p2 <- rbind(ITC.param.p2, param.p2.neg)
}

ITC.param.p2.FT2 <- subset(ITC.param.p2, Xs >= ((Xlast - X0)*0.9 + X0))
length(unique(ITC.param.p2.FT2$ANIMAL_ID))

#Remove animals have negative slope of DFI out of dataset for Function type 3
ITC.param.pos2 <- subset(ITC.param.pos2, !(ANIMAL_ID %in% ITC.param.p2.FT2$ANIMAL_ID))
Age.remain.pos2 <- subset(Age.remain.pos2, !(ANIMAL_ID %in% ITC.param.p2.FT2$ANIMAL_ID))

#Add these animals to the group of QDR function for CFI
DFI.pos11 <- c( DFI.pos1, unique(ITC.param.p2.FT2$ANIMAL_ID))
DFI.pos1 <- unique(DFI.pos11)
# DFI.pos1 <- as.numeric(as.character(DFI.pos1))

length(DFI.pos1)

#==============================================================================
# Regrouping the animals have error in Q-L function to QDR function
#==============================================================================

#The error pigs in Quadratic-linear function is identified
Last.pigs <- subset(No.NA.Data.1, !(ANIMAL_ID %in% DFI.neg)); length(unique(Last.pigs$ANIMAL_ID))   
Last.pigs <- subset(Last.pigs, !(ANIMAL_ID %in% DFI.pos1)); length(unique(Last.pigs$ANIMAL_ID)) 
Last.pigs <- subset(Last.pigs, !(ANIMAL_ID %in% ITC.param.pos2$ANIMAL_ID)); length(unique(Last.pigs$ANIMAL_ID))  

Last.pigs <- unique(Last.pigs$ANIMAL_ID)

#Add these animals to the group of QDR function for CFI
DFI.pos11 <- c( DFI.pos1, Last.pigs)
DFI.pos1 <- unique(DFI.pos11)
# DFI.pos1 <- as.numeric(as.character(DFI.pos1))


#Reset the animals in group of Quadratic-linear function for CFI
DFI.pos2 <- unique(Age.remain.pos2$ANIMAL_ID)

length(DFI.pos1); length(DFI.neg); length(DFI.pos2)
  
#===============================================================
#Save results to Rdata file
#===============================================================

  save(ITC.param.pos2, Age.remain.pos2,
       file = "Data/JRP.DFI.pos2.RData")

  save(No.NA.Data.1, DFI.neg, DFI.pos1, DFI.pos2,
        file = "Data/JRPData_TTC.RData")
# DFI.neg <- factor(DFI.neg1, levels = 1:nlevels(DFI.pos2), labels = levels(DFI.pos2))
# DFI.pos1 <- factor(DFI.pos11, levels = 1:nlevels(DFI.pos2), labels = levels(DFI.pos2))