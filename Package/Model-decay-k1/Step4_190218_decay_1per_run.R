##################################################################################              
#               Mycotoxin data - STEP 4: Quantify response to perturbation
#               Hieu 12/02/2019
#               Example of QDR function for target CFI
#               Fit perturbation model with 01, 02 and 03 perturbations to data
#               Compare these models to see which model is the best 
#               This code is for a model with 01 perturbation
#               Replace t_stop1 by decaying rate
##################################################################################

  rm(list=ls(all=TRUE)) #To remove the hisory 
  dev.off() #To close all graphs

#-------------------------------------------------------------------------------
#    Needed Packages
#-------------------------------------------------------------------------------
    library(deSolve)
    library(nlstools)
    library(proto)
    library(nls2)
    library(ggplot2)
    library(dplyr)
        
################################################################

  # Set working directory  
    setwd("C:/Users/hnguyenba/Dropbox/INRA/Modelling_Perturbation/Hieu_PHD_shared/Data/Finished-R-scripts/All steps/Mycotoxin/Step4_2_per/Test_multi_per")
  
#===============================================================
# DATA PREPARATION
#===============================================================

  # Load needed data and packages   
  load("DONData.QLM1.Rdata")
  source("Step4_190218_decay_1per_functions.R")
  source("Step4_190212_functions_target-CFI.R")
  
  # Create empty data frame to store results   
  Per.para <- data.frame()
  
  # Animal_ID
  ID <- No.NA.Data.TC %>% group_by(FuncType) %>% filter(FuncType == "QDR") %>% distinct(ANIMAL_ID)
  ID <- ID$ANIMAL_ID

 #-------------------------------------------------------------------------------
 # Extract data for animal i
 #-------------------------------------------------------------------------------
# for(idc in seq_along(ID)){
  idc <- 9
  i <- ID[idc] # order number of one animal
  Data <- No.NA.Data.TC[No.NA.Data.TC$ANIMAL_ID == i,]
  pertub.table <- fn.pertub.table[fn.pertub.table$ANIMAL_ID == i, ] 
  res.data <- fn.res[fn.res$ANIMAL_ID == i, ]
  
  Age <- Data$Age.plot
  DFI.obs <- Data$DFI.plot
  CFI.obs <- Data$CFI.plot
  
  # Difference between actual and target CFI (kg)    
  res <- res.data$res
  
 #-------------------------------------------------------------------------------
 # Calulate target CFI using suitable function
 #-------------------------------------------------------------------------------

     param.i <- ITC.param.pos1[ITC.param.pos1$ANIMAL_ID == i,]
     TTC.param <- param.i
     FuncType <- param.i[dim(param.i)[1],]$FuncType
     Slope <- param.i[dim(param.i)[1],]$Slope
     param.i <- as.numeric(param.i[dim(param.i)[1],c(5:7)])
     ITC <- pred.abcd.1(param.i, Age)[[1]]
     ITD <- pred.abcd.1(param.i, Age)[[2]]
  
     ##-----------------------------
     ## initial values and times
     ##-----------------------------
     a <- TTC.param$a
     b <- TTC.param$b
     c <- TTC.param$c
     
#=================================== ================== ================== ==================
#  Modelling the response of the animal to a single perturbation
#=================================== ================== ================== ==================
     
     Data.xy <- Data
     times <- Data.xy$Age.plot
     
     p1.init <- 0.7
     
     yinit <- c(CumFI = ITC[1], p1 = p1.init) #state
     times.ode <- seq(from = times[1], to = times[length(times)], by = .1)
    
     # Find initial values for p1 and maxcompFI1
     st1 <- expand.grid(p1 = seq(0, 1, len = 4),
                        p2 = seq(1, 10, len = 4),
                        p3 = seq(times[1]+5, times[length(times)]-10, len = 4),
                        p4 = seq(0, 0.1, len = 4))

     # # Use NLS2 to select the best initial parameters
     st2 <- nls2(CFI.plot ~ ODE.CFI.obj.nls.1(p1, p2, p3, p4),
                 Data.xy,
                 start = st1,
                 algorithm = "brute-force")

     ##-----------------------------
     ## parameters
     ##-----------------------------
     # Extract values for initial parameters
     p1.init <- as.numeric(coef(st2)[1]); p1.init
     max.compFI1 <- as.numeric(coef(st2)[2]); max.compFI1
     tbeg1 <- as.numeric(coef(st2)[3]); tbeg1
     decay1 <- as.numeric(coef(st2)[4]); decay1
     
     param <- c(p1.init, max.compFI1, tbeg1, decay1, a, b, c)
     yinit <- c(CumFI = ITC[1], p1 = p1.init) #state
     
     ##-----------------------------
     ## solve the model
     ##-----------------------------
     
     yout <- ode(y= yinit,
                 times = times.ode,
                 func = ODE.CFI.1,
                 parms = param,
                 atol = 1e-10,
                 rtol = 1e-10)
     summary(yout)
     
     # plot(yout[,1], yout[,2], type = "l", lwd = 2)
     # points(Age, CFI.obs, type = "p")
     plot(yout[,1], yout[,3], type = "l", lwd = 2)
     
     # ##--------------------------------------------------
     # ## plot when fitting initial parameters to the model
     # ##--------------------------------------------------
     # 
     Time <- times.ode
     onoff <- ifelse(Time>tbeg1, 1, 0)
     ITC.sim <- pred.abcd.1(param.i, Time)[[1]]

     #Difference between actual and target CFI (simulation)
     plot(yout[ , 1], yout[ ,2] -(ITC.sim), type = "l",
          xlab = "Age (days)", ylab = "CFI - TTC(kg)",
          ylim = c(min(CFI.obs - ITC, yout[ ,2] -(ITC.sim)), max(CFI.obs - ITC, yout[ ,2] -(ITC.sim))),
          lwd = 2,
          cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
     points ( Age,  CFI.obs - ITC, type = "p", col = "red")
     # if(unique(Data$lot) == "DC"){
     #   abline(v=113, lty=2)
     #   abline(v=119, lty=2)
     # } else if(unique(Data$lot) == "CD"){
     #   abline(v=134, lty=2)
     #   abline(v=140, lty=2)
     # } else if(unique(Data$lot) == "DD"){
     #   abline(v=113, lty=2)
     #   abline(v=119, lty=2)
     #   abline(v=134, lty=2)
     #   abline(v=140, lty=2)
     # } else{}
     abline(0,0, col = "red")

     # plot(Age, ITC, type = "l", lwd = 2.5, col = "blue")
     # points(yout[,1], yout[,2], type = "l", lwd = 2.5, col = "red")
     # points(Age, CFI.obs, cex = 1.2)
     # abline(v=tbeg1, lty=2)
     
     # ==============  Estimation of parameters ================== ==================
     
     ##-----------------------------
     ## initial values and times
     ##-----------------------------
     
     yinit <- c(CumFI = ITC[1], p1 = p1.init) #state
     times.ode <- seq(from = times[1], to = times[length(times)], by = 1)
     
     ##-----------------------------
     ## parameters
     ##-----------------------------
     
     par.init <- c(p1.init, max.compFI1, tbeg1, decay1)
     
     ##-----------------------------
     ## solve the model
     ##-----------------------------
     
     yout <- ode(y= yinit,
                 times = times.ode,
                 func = ODE.CFI.optim.1,
                 parms = par.init,
                 atol = 1e-10,
                 rtol = 1e-10)
     
     summary(yout)
     
     ##-----------------------------
     #  run the optimization with NLS2
     ##-----------------------------
     # Data.xy = Data
     times <- Data.xy$Age.plot
     ODE.CFI.obj.1(par.init, Data.xy)
     
     #Estimate parameters by Optim and the best initial parameters
     optim.res <- optim(par.init, ODE.CFI.obj.1,
                        lower = c(0, 1, 0, 0),
                        upper = c(1, 1000, 1000, 1000), method="L-BFGS-B",
                        control = list(trace=6),
                        hessian = TRUE)
     
     P.optim <- c(optim.res$par[1], optim.res$par[2], optim.res$par[3], optim.res$par[4])
     P.optim
     
     p1.init <- P.optim[1]
     max.compFI1 <- P.optim[2]
     tbeg1 <- P.optim[3]
     decay1 <- P.optim[4]
     
     yinit <- c(CumFI = ITC[1], p1 = p1.init) #state
     
     RSS <- ODE.CFI.obj.1(P.optim, Data.xy); RSS
     
     # Simulate the ode function
     yout   <- ode(yinit, times, ODE.CFI.optim.1, P.optim)
     
     #Test graph
     Time <- times
     onoff <- ifelse(Time>tbeg1, 1, 0)
     ITC.sim <- pred.abcd.1(param.i, Time)[[1]]
     
     #Difference between actual and target CFI (simulation)
     # png(filename = paste("Graphs/",idc, ".",i ,".",unique(Data$lot), ".1-per",".png",sep=""),
     #     height = 720, width = 1000, units = 'px', type="cairo-png")
     par(mar=c(4.5,4.5,3,1.5))
     plot(yout[ , 1], yout[ ,2] -(ITC.sim), type = "l",
          main = paste("Group ",unique(Data$lot), " ,", "Pig", i, "\nCFI - target CFI", ", RSS = ", round(RSS,2)),
          xlab = "Age (days)", ylab = "Kg", 
          ylim = c(min(CFI.obs - ITC, yout[ ,2] -(ITC.sim)), max(CFI.obs - ITC, yout[ ,2] -(ITC.sim))),
          lwd = 2,
          cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
     points ( Age,  CFI.obs - ITC, type = "p", col = "red")
     # if(unique(Data$lot) == "DC"){
     #   abline(v=113, lty=2)
     #   abline(v=119, lty=2)
     # } else if(unique(Data$lot) == "CD"){
     #   abline(v=134, lty=2)
     #   abline(v=140, lty=2)
     # } else if(unique(Data$lot) == "DD"){
     #   abline(v=113, lty=2)
     #   abline(v=119, lty=2)
     #   abline(v=134, lty=2)
     #   abline(v=140, lty=2)
     # } else{}
     abline(0,0, col = "red")
     par(mfrow = c(1,1))
     # dev.off()
  #-------------------
  # Add perturbation parameters to a data frame
  #-------------------
     
     # Report the estimated parameters
     k1 <- -(p1.init)*100
     k2 <- max.compFI1
     
     ANIMAL_ID <- Data %>% distinct(ANIMAL_ID)
     lot <- Data %>% distinct(lot)
     
     Per.para1 <- data.frame(cbind(lot, ANIMAL_ID, tbeg1, k1, decay1, k2))
     Per.para1$RSS <- RSS
     names(Per.para1) <- c("lot", "ANIMAL_ID","t_start", "k1", "decay1", "k2", "RSS")
     Per.para1
     
     # ------------------------
     # Gather data for all pigs
     # ------------------------
     Per.para <- rbind(Per.para, Per.para1)
     Per.para
     
     #-------------------
     #Prepare for PLOT
     #-------------------
     times.ode <- seq(from = Age[1], to = Age[length(Age)], by = 1)
     Time <- times
     onoff <- ifelse(Time>tbeg1, 1, 0)

     if(FuncType == "LM"){
       ITC <- pred.abcd.0(param.i, times)[[1]]
       ITD <- rep(pred.abcd.0(param.i, times)[[2]], length(times))
     } else if(FuncType == "QDR"){
       ITC <- pred.abcd.1(param.i, times)[[1]]
       ITD <- pred.abcd.1(param.i, times)[[2]]
     } else{
       ITC <- pred.abcd.2(param.i, times)[[1]]
       ITD <- pred.abcd.2(param.i, times)[[2]]
     }

     # #Compensatory feed intake
     CompFI <- (1-yout[, 2]/ITC)*max.compFI1
     #Simulation of CFI
     CFI.sim <- yout[, 2]
     #Impact of perturbation
     p1 <- yout[, 3]
     #Simulation of DFI
     DFI.sim <- (-(onoff*p1) + CompFI + 1)*ITD
     #Ratio of DFI
     Ratio.DFI <- DFI.sim/ITD

     Age.remain <- Data.remain[Data.remain$ANIMAL_ID == i,]$Age
     DFI.remain <- Data.remain[Data.remain$ANIMAL_ID == i,]$Observed_DFI
     CFI.remain <- Data.remain[Data.remain$ANIMAL_ID == i,]$Observed_CFI

     #-------------------
     # PLOT THE GRAPH
     #-------------------
     #DFI
     # png(filename = paste("Graphs/Step4_graphs/",unique(Data$lot),idc,".",i, ".DFI",".png",sep=""),
     #     height = 720, width = 1000, units = 'px', type="cairo-png")
     par(mar=c(4.5,4.5,3,1.5))
     plot(Age, DFI.obs,
          main = paste("Pig", i, "\nResponse to perturbation in DFI"),
          ylab = "Daily feed Intake (kg/d)",
          xlab = "Age (days)",
          ylim = c(min(ITD, DFI.obs, DFI.sim), max(ITD, DFI.obs, DFI.sim)),
          type="p", col ="black",
          cex.main = 1.5, cex = 1.8,
          cex.axis = 1.5, cex.lab = 1.5)
     points(times, ITD, type = "l", col = "blue", lwd = 4)
     points(Age.remain, DFI.remain, type = "p",col = "blue", pch = 19, cex = 1.8)
     points(times, DFI.sim, type = "l", col = "red", lwd = 4)

     #CFI
     # png(filename = paste("Graphs/Step4_graphs/",unique(Data$lot),idc,".",i, ".CFI",".png",sep=""),
     #     height = 720, width = 1000, units = 'px', type="cairo-png")
     par(mar=c(4.5,4.5,3,1.5))
     plot(Age, CFI.obs,
          main = paste("Pig", i, "\nResponse to perturbation in CFI"),
          ylab = "Cumulative feed Intake (kg)",
          xlab = "Age (days)",
          ylim = c(min(ITC, CFI.obs, CFI.sim), max(ITC, CFI.obs, CFI.sim)),
          type="p", col ="black",
          cex.main = 1.5, cex = 1.2,
          cex.axis = 1.5, cex.lab = 1.5)
     points(times, ITC, type = "l", col = "blue", lwd = 3)
     points(Age.remain, CFI.remain, type = "p",col = "blue", pch = 19, cex = 1.2)
     points(times, CFI.sim, type = "l", col = "red", lwd = 3)
     # if(unique(Data$lot) == "DC"){
     #   abline(v=113, lty=2)
     #   abline(v=119, lty=2)
     # } else if(unique(Data$lot) == "CD"){
     #   abline(v=134, lty=2)
     #   abline(v=140, lty=2)
     # } else if(unique(Data$lot) == "DD"){
     #   abline(v=113, lty=2)
     #   abline(v=119, lty=2)
     #   abline(v=134, lty=2)
     #   abline(v=140, lty=2)
     # } else{}
     # dev.off()

     #Compensatory feed intake and perturbation effect
     # png(filename = paste("Graphs/Step4_graphs/",unique(Data$lot),idc,".",i, ".ratio",".png",sep=""),
     #     height = 720, width = 1000, units = 'px', type="cairo-png")
     par(mar=c(4.5,4.5,3,1.5))
     plot(Time, CompFI*100,
          main = paste("Pig", i, ", k1 =", round(k1, 1), ", k2 =", round(k2, 2)),
          ylab = "Change in DFI over time (%)",
          xlab = "Age (days)",
          ylim = c(-100, 100),
          # ylim = c(min(CompFI*100, (0- onoff*(1-p1))*100), max(CompFI*100, (0- onoff*(1-p1))*100)),
          type = "l", lwd = 4, col = "green", cex.main = 1.5,
          cex.axis = 1.5, cex.lab = 1.5)
     points(Time, (0- onoff*(p1))*100, type = "l", lwd = 4, col = "purple")
     # legend("bottomleft", c("k1", "k2"),
     #        col=c("green", "purple"), lty = c(1,1), bty = "n", cex= 1.8, lwd = 2)
     # dev.off()
# }
#   One.Per.para <- Per.para
#   save(One.Per.para, file = "One.Per.para.Rdata")
  