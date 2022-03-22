  #==========================================================================
  #
  # Model ode.pert.14/11 - Jaap and Hieu - 14/11 /2017
  #
  # DESCRIPTION: Model with p1 constant and orignial maxcompFI (without reparametrization)
  #
  #==========================================================================
 

  rm(list=ls(all=TRUE)) #To remove the hisory 
  # dev.off() #To close all graphs

  #==========================================================================

  setwd("C:/Users/Kevin Le/PycharmProjects/Pig Data Black Box")
  # setwd("/Users/nguyenbahieu/Dropbox/INRA/Modelling_Perturbation/Hieu_PHD_shared/Data/Finished-R-scripts/All\ steps/individual/")
  
  #-------------------------------------------------------------------------------
  # Import data file and functions
  #-------------------------------------------------------------------------------
    
   #Packages
    library(deSolve)
    library(nlstools)
    library(nlsMicrobio)
    library(proto)
    library(nls2)
    library(ggplot2)
    library(dplyr)

        
  #-------------------------------------------------------------------------------
  # Import data file and functions
  #-------------------------------------------------------------------------------
  
     load("Data/JRP.Per.detec.RData") #load dataset created in MissingData.step
     load("Data/JRPData.Rdata")
     load("Data/JRPData_TTC.RData")
     load("Data/JRP.DFI.pos1.RData")
     load("Data/JRP.DFI.pos2.RData")
     load("Data/JRP.DFI.neg.RData")
     source("Package/Functions.R") #functions
     source("Package/Step4_functions_X0_problem.R")
     options(digits=3) 

  #===============================================================
  # DATA PREPARATION
  #===============================================================

  #Order number of Animal_ID
  # ID <- unique(as.factor(No.NA.Data$ANIMAL_ID))
  ID <- 6001
  # ID <- c(5775,6001)
      
  #-------------------------------------------------------------------------------
  # Extract data for animal i
  #-------------------------------------------------------------------------------
  AF <- list()
  AF <- NULL
  IDC <- seq_along(ID)
  for (idc in IDC){
    i <- ID[idc]
    Data <- No.NA.Data.1[No.NA.Data.1$ANIMAL_ID == i,]
    CFI.obs <- Data$CFI.plot
    DFI.obs <- Data$DFI.plot
    Age <- Data$Age.plot
  
    #Difference between actual and target CFI (kg)
    res <- fn.res$res[fn.res$ANIMAL_ID == i]

    #Information of TTC function
    # TTC.param <- merge(ITC.param.pos2,ITC.param.pos1,ITC.param.neg,by = "ANIMAL_ID")
    if (i %in% unique(ITC.param.neg$ANIMAL_ID)) {
        TTC.param <- ITC.param.neg[ITC.param.neg$ANIMAL_ID == i,]
        FuncType <- TTC.param$FuncType; FuncType
        Slope <- TTC.param$Slope

      } else if (i %in% unique(ITC.param.pos1$ANIMAL_ID)){
        TTC.param <- ITC.param.pos1[ITC.param.pos1$ANIMAL_ID == i,]
        FuncType <- TTC.param$FuncType; FuncType
        Slope <- TTC.param$Slope

      } else {
        TTC.param <- ITC.param.pos2[ITC.param.pos2$ANIMAL_ID == i,]
        FuncType <- TTC.param$FuncType; FuncType
        Slope <- TTC.param$Slope

      }
    ## ======================= TTC ======================================================
      
    #-------------------------------------------------------------------------------
    # Calulate TTC using abcd function
    #-------------------------------------------------------------------------------

    if(FuncType == "LM"){
      param.i <- as.numeric(TTC.param[dim(TTC.param)[1], 4:5])
      ITC <- pred.abcd.0(param.i, Age)[[1]]
      ITD <- rep(param.i[2], length(Age))

    } else if(FuncType == "QDR"){
      param.i <- as.numeric(TTC.param[dim(TTC.param)[1], 5:7])
      ITC <- pred.abcd.1(param.i, Age)[[1]]
      ITD <- pred.abcd.1(param.i, Age)[[2]]

    } else{
      param.i <- as.numeric(TTC.param[dim(TTC.param)[1], 6:8])
      Xs <- TTC.param[dim(TTC.param)[1],]$Xs
      ITC <- pred.abcd.2(param.i, Age)[[1]]
      ITD <- pred.abcd.2(param.i, Age)[[2]]

    }
    # #Magnitude of the perturbation
    magnitude <- fn.res[fn.res$ANIMAL_ID == i,]
    #   Age_start <-  fn.pertub.table$Start[fn.pertub.table$ANIMAL_ID == i]
    #   Age_end <- fn.pertub.table$End[fn.pertub.table$ANIMAL_ID == i]
    # magnitude <-  filter(magnitude$Age >= Age_start & magnitude$Age <= Age_end)
    # # magnitude   = res.data %>% filter(Age >= pertub.table$Start[1] & Age <= pertub.table$End[1])
    # magnitude <- magnitude %>% filter(res == min(magnitude$res))
    #-------------------------------------------------------------------------------
    #Plot
    #-------------------------------------------------------------------------------

    B <- dev.df[!is.na(dev.df$ppert),]
    #
    plot(Age, res,
        main = paste( "Pig ID:", ID, "\nDifference between CFI and TTC"),
        xlab = "Age (days)",
        ylab = "Amount of difference: CFI - TTC (kg)",
        type = "o", pch = 10, cex = 0.5,
        ylim = c(min(B$dif.CFI, res),
                 max(B$dif.CFI, res)),
       cex.main = 1.5, cex.lab = 1.2)
    abline(0,0, col = "red")
    points(fn.pertub.table$Start, fn.pertub.table$value.start,
          col = "orange", pch=15, cex = 2)
    points(fn.pertub.table$End, fn.pertub.table$value.end,
          col = "green", pch=17, cex = 2)
    points(magnitude$Age, magnitude$res,
           col = "purple", pch=19, cex = 2)

    #Read pertubation data
    fn.table.input <- fn.pertub.table[fn.pertub.table$ANIMAL_ID == i,]
    for (ii in fn.table.input$ppert){
      staging <- fn.table.input[fn.table.input$ppert == ii, ]
      start <- staging$Start
      end <- staging$End
      #Magnitude of the perturbation
      seq_age <- round(seq(start,end))
      magnitude <- magnitude %>% filter(magnitude$Age >= start & magnitude$Age <= end)
      # magnitude   = res.data %>% filter(Age >= pertub.table$Start[1] & Age <= pertub.table$End[1])
      magnitude <- magnitude %>% filter(res == min(magnitude$res))


      #=================================== ================== ================== ==================
      #  Modelling the response of the animal to a single perturbation
      #=================================== ================== ================== ==================

      #------------------------------
      # Linear function for TTC
      #------------------------------
      if(FuncType == "LM"){

        ##-----------------------------
        ## initial values and times
        ##-----------------------------
        a <- TTC.param[dim(TTC.param)[1],]$a
        b <- TTC.param[dim(TTC.param)[1],]$b
        tbeg1 <- start
        tstop1 <- end

        # Create a grid to choose the best initial values for k1 and k2
        Data.xy <- Data
        times <- Data.xy$Age.plot
        yinit <- c(CumFI = ITC[1]) #state
        times.ode <- seq(from = Age[1], to = Age[length(Age)], by = .1)

        st1 <- expand.grid(p1 = seq(0, 1, len = 10),
                           p2 = seq(1, 10, len = 10))

        # # Use NLS2 to select the best initial parameters
        st2 <- nls2(CFI.plot ~ ODE.CFI.obj.nls.0(p1,p2),
                    Data.xy,
                    start = st1,
                    algorithm = "brute-force")

        p1 <- coef(st2)[1]
        max.compFI1 <- coef(st2)[2]



        ##-----------------------------
        ## parameters
        ##-----------------------------

        par.init <- c(p1, max.compFI1)

        ##-----------------------------
        #  run the optimization with NLS2
        ##-----------------------------
        # ODE.CFI.obj.0(par.init, Data.xy)

        #Estimate parameters by Optim and the best initial parameters
        optim.res <- optim(par.init, ODE.CFI.obj.0,
                           lower = c(0, 1), upper = c(1,1000), method="L-BFGS-B",
                           hessian = TRUE)

        P.optim <- c(optim.res$par[1], optim.res$par[2])

        RSS <- ODE.CFI.obj.0(P.optim, Data.xy)

        # Simulate the ode function
        yout   <- ode(yinit, times, ODE.CFI.optim.0, P.optim)

      } else{}

      #----------------------------------
      # Quadratic function for TTC
      #----------------------------------

      if(FuncType == "QDR"){

        ##-----------------------------
        ## initial values and times
        ##-----------------------------
        tbeg1 <- start
        # tbeg1       = pertub.table$Start # tbeg1 corresponds to the start point
        tstop1 <- end
        a <- TTC.param[dim(TTC.param)[1],]$a
        b <- TTC.param[dim(TTC.param)[1],]$b
        c <- TTC.param[dim(TTC.param)[1],]$c

        # Create a grid to choose the best initial values for k1 and k2
        Data.xy <- Data
        times <- Data.xy$Age.plot
        yinit <- c(CumFI = ITC[1]) #state
        times.ode <- seq(from = Age[1], to = Age[length(Age)], by = .1)

        st1 <- expand.grid(p1 = seq(0, 1, len = 10),
                           p2 = seq(1, 10, len = 10))
        #
        # # Use NLS2 to select the best initial parameters
        st2 <- nls2(CFI.plot ~ ODE.CFI.obj.nls.1(p1, p2),
                      Data.xy,
                      start = st1,
                      algorithm = "brute-force")

        p1 <- coef(st2)[1]
        max.compFI1 <- coef(st2)[2]
        ##-----------------------------
        ## parameters
        ##-----------------------------

        param <- c(p1, max.compFI1)

        ##-----------------------------
        ## solve the model
        ##-----------------------------

        # yout <- ode(y = yinit,
        #             times = times.ode,
        #             func = ODE.CFI.1,
        #             parms = param,
        #             atol = 1e-10,
        #             rtol = 1e-10)
        # summary(yout)

        ##--------------------------------------------------
        ## plot when fitting initial parameters to the model
        ##--------------------------------------------------

        Time <- times.ode
        onoff <- ifelse(Time>tbeg1 & Time<tstop1, 1, 0)
        ITC.sim <- pred.abcd.1(param.i, Time)[[1]]

        #Difference between actual and target CFI (simulation)
        # plot(yout[ , 1], yout[ ,2] -(ITC.sim),
        #        xlab = "Age (days)", type = "l", ylab = "CFI - TTC(kg)",
        #        ylim = c(min(CFI.obs - ITC, yout[ ,2] -(ITC.sim)), max(CFI.obs - ITC, yout[ ,2] -(ITC.sim))),
        #        lwd = 2,
        #        cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
        # points ( Age,  CFI.obs - ITC, type = "p", col = "red")
        # abline(v=tbeg1, lty=2)
        # abline(v=tstop1, lty = 2)
        # abline(0,0, col = "red")
        # par(mfrow = c(1,1))
        # dev.off()

        # ==============  Estimation of parameters ================== ==================

        ##-----------------------------
        ## parameters
        ##-----------------------------

        par.init <- c(p1, max.compFI1)

        ##-----------------------------
        #  run the optimization with NLS2
        ##-----------------------------
        ODE.CFI.obj.1(par.init, Data.xy)

        #Estimate parameters by Optim and the best initial parameters
        optim.res <- optim(par.init, ODE.CFI.obj.1,
                       lower = c(0, 1), upper = c(1,1000), method="L-BFGS-B",
                       hessian = TRUE)

        P.optim <- c(optim.res$par[1], optim.res$par[2])

        RSS <- ODE.CFI.obj.1(P.optim, Data.xy)

        # Simulate the ode function
        yout   <- ode(yinit, times, ODE.CFI.optim.1, P.optim)


        } else{}

      #----------------------------------
      # Quadratic-linear function for TTC
      #----------------------------------
      if(FuncType == "QLM"){

        ##-----------------------------
        ## initial values and times
        ##-----------------------------

        tbeg1 <- start # tbeg1 corresponds to the start point
        tstop1 <- end
        a <- TTC.param[dim(TTC.param)[1],]$a
        b <- TTC.param[dim(TTC.param)[1],]$b
        c <- TTC.param[dim(TTC.param)[1],]$c
        Xs <- TTC.param[dim(TTC.param)[1],]$Xs

        # Create a grid to choose the best initial values for k1 and k2
        Data.xy <- Data
        times <- Data.xy$Age.plot
        yinit <- c(CumFI = ITC[1]) #state
        times.ode <- seq(from = Age[1], to = Age[length(Age)], by = .1)

        st1 <- expand.grid(p1 = seq(0, 1, len = 10),
                           p2 = seq(1, 10, len = 10))
        #
        # # Use NLS2 to select the best initial parameters
        st2 <- nls2(CFI.plot ~ ODE.CFI.obj.nls.2(p1, p2),
                Data.xy,
                start = st1,
                algorithm = "brute-force")

        p1 <- coef(st2)[1]
        max.compFI1 <- coef(st2)[2]

        ##-----------------------------
        ## parameters
        ##-----------------------------

        param <- c(p1, max.compFI1)

        ##-----------------------------
        ## solve the model
        ##-----------------------------

        # yout <- ode(y = yinit,
        #             times = times.ode,
        #             func = ODE.CFI.2,
        #             parms = param,
        #             atol = 1e-10,
        #             rtol = 1e-10)
        # summary(yout)


        ##--------------------------------------------------
        ## plot when fitting initial parameters to the model
        ##--------------------------------------------------

        # Time <- seq(from = Age[1], to = Age[length(Age)], by = .1)
        # onoff <- ifelse(Time>tbeg1 & Time<tstop1, 1, 0)
        # ITC.sim <- pred.abcd.2(param.i, Time)[[1]]

        # par(mar=c(4.5,4.5,4.5,1.5))

        #Difference between actual and target CFI (simulation)
           # plot(yout[ , 1], yout[ ,2] -(ITC.sim),
           #      xlab = "Age (days)", type = "l", ylab = "CFI - TTC(kg)",
           #      ylim = c(min(CFI.obs - ITC, yout[ ,2] -(ITC.sim)), max(CFI.obs - ITC, yout[ ,2] -(ITC.sim))),
           #      lwd = 2,
           #      cex.main = 1.5, cex.axis = 1.5, cex.lab = 1.5)
           # points ( Age,  CFI.obs -( ITC), type = "p", col = "red")
           # abline(v=tbeg1, lty=2)
           # abline(v=tstop1, lty = 2)
           # abline(0,0, col = "red")
           # par(mfrow = c(1,1))
           # dev.off()

        # ==============  Estimation of parameters ================== ==================

        ##-----------------------------
        ## parameters
        ##-----------------------------

        par.init <- c(p1, max.compFI1)

        ##-----------------------------
        #  run the optimization with NLS2
        ##-----------------------------
        ODE.CFI.obj.2(par.init, Data.xy)
        # Estimate parameters by Optim and the best initial parameters
        optim.res <- optim(par.init, ODE.CFI.obj.2,
                           lower = c(0, 1), upper = c(1,1000), method="L-BFGS-B",
                           hessian = TRUE)

        P.optim <- c(optim.res$par[1], optim.res$par[2])

        RSS <- ODE.CFI.obj.2(P.optim, Data.xy)

        # Simulate the ode function
        yout   <- ode(yinit, times, ODE.CFI.optim.2, P.optim)

      } else{}

      ###############################################################################
      ##-----------------------------
      #  Prepare for plotting
      ##-----------------------------

      p1 <- P.optim[1]
      max.compFI1 <- P.optim[2]
      # tbeg1 <- P.optim[3]
      # tstop1 <- P.optim[4]

      Time <- times
      onoff <- ifelse(Time>start & Time<end, 1, 0)
      #Target trajectory curve
      ITC; ITD
      #Compensatory feed intake
      CompFI <- (1-yout[, 2]/ITC)*max.compFI1
      #Simulation of DFI
      DFI.sim <- (onoff*(p1-1) + CompFI + 1)*ITD
      #Ratio of DFI
      Ratio.DFI <- DFI.sim/ITD
      dif <- yout[,2]/ITC
      #-----------------------------
      # plot the results
      #-----------------------------
      #Ratio between actual and target CFI
      plot(times, dif,
           main = "Ratio between CFI and ITC_CFI", type = "l",
           xlab = "Age (days)", ylab= "Ratio", col="black",
           # ylim = c(0.94,1),
           lwd = 2)
      abline(v=P.optim[3], lty=2)
      abline(v=P.optim[4], lty=2)

      #Compensatory feed intake
      par(mfrow = c(1,1))
      # plot(Time, onoff, type = "l", ylab = "onoff")
      plot(times, CompFI,
           main = "Compensatory of feed intake", type = "l", ylab= "compensatory (%)", xlab = "Age (days)",
           col="black",
           # ylim = c(1, 1.4),
           lwd = 2)
      abline(v=P.optim[3], lty=2)
      abline(v=P.optim[4], lty=2)

      #Impact of p1
      plot(Time, 0- onoff*(1-p1),
           main = "Impact of p1", type = "l", ylab= "Reduction in DFI (%)", col="black")

      #Ratio between actual and target CFI

      plot(times, yout[,2]/ITC,
           main = "Ratio between CFI and ITC_CFI", type = "l",
           xlab = "Age (days)", ylab= "Ratio", col="black",
           # ylim = c(0.94,1),
           lwd = 2)
      abline(v=P.optim[3], lty=2)
      abline(v=P.optim[4], lty=2)

      #MAIN PLOTS
      par(mar=c(4.5,4.5,4.5,1.5))
        e <- CFI.obs - ITC
      dev.new()
      #dif
      plot(Age, CFI.obs - ITC,
           main = paste("Difference between CFI and ITC_CFI", "\nPig ID:", ID),
           xlab = "Age (days)", ylab = "CFI - ITC (kg)",
           ylim = c(min(CFI.obs - ( ITC),  yout[ ,2] -( ITC)),
                    max(CFI.obs - ( ITC),  yout[ ,2] -( ITC))),
           type = "p", col = "red",
           cex.main = 1.7, cex.axis = 1.5, cex.lab = 1.5)
      abline(0,0, col = "blue", lwd = 2)
      points ( yout[ , 1],  yout[ ,2] -ITC,type = "l", lwd = 2.2)
      abline(v=P.optim[3], lty=3)
      abline(v=P.optim[4], lty=3)
      legend("bottomleft",legend = c("ITC", "Dif-obs", "Dif-simul"), col = c("blue","red", "black"),
             lty = c(1,NA, 1), pch = c(NA, 1, NA), bty = "n", lwd = 3, cex = 1.5, pt.cex = 2.5)
      dev.off()
      par(mfrow = c(1,1))

      #ggplot2
      #Compensatory feed intake and perturbation effect
      AF1 <- data.frame(cbind(times, CompFI*100))
      names(AF1) <- c("Age.plot", "Percent")
      AF1 <- AF1 %>% mutate(Curve = rep("Resilience", length(times)))
      # if (AF1$Percent < 0 ){
      #       AF1$Percent <- 0
      #     }
      AF2 <- data.frame(cbind(times, (0- onoff*(1-p1))*100))
      names(AF2) <- c("Age.plot", "Percent")
      AF2 <- AF2 %>% mutate(Curve = rep("Resistance", length(times)))

      AF3 <- data.frame(cbind(times, rep(0, length(times))))
      names(AF3) <- c("Age.plot", "Percent")
      AF3 <- AF3 %>% mutate(Curve = rep("AA", length(times)))
      AF <- rbind(AF1, AF2, AF3)




      # ggsave(file = paste0("Graphs/Step4_graphs/", Data$ANIMAL_ID, ".", "Ratio",ii, ".png"), width = 6000, height = 3500, units = "px", res=600)
      cols.CFI <- c("Resilience" = "green", "Resistance" = "purple", "AA" = "blue")
      type.CFI <- c("Resilience" = "solid", "Resistance" = "solid", "AA" = "dashed")
      size.CFI <- c("Resilience" = 1.3, "Resistance" = 2.5, "AA" = 1)
      plot_res_model <- ggplot(data = AF, aes(x = Age.plot, y = Percent)) +
        geom_hline(yintercept=0, linetype="dashed", color = "blue", size = 1) +
        geom_line(aes(color = Curve, linetype=Curve, size = Curve)) +
        scale_color_manual(values = cols.CFI, labels = c("Resilience","Resistance","Target")) +
        scale_linetype_manual(values = type.CFI, labels = c("Resilience","Resistance","Target")) +
        scale_size_manual(values = size.CFI, labels = c("Resilience","Resistance","Target")) +
        xlab("Age, d") +
        ylab("Change in daily feed intake (%)") +
        # expand_limits(y=-100) +
        # expand_limits(y=100) +
        scale_y_continuous(breaks=seq(-100, 100, 10)) +
        scale_x_continuous(breaks=seq(60, 240, 20)) +
        ggtitle(paste("Change in DFI of", "the pig", i, "\nduring and after a single perturbation")) +
        theme(plot.title = element_text(hjust = .5)) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black")
        ) +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14,face="bold"))
      ggsave(file = paste0("Graphs/Step4_graphs/Pertubation Model/", Data$ANIMAL_ID, ".", "Ratio",ii, ".png"),plot_res_model,device = "png",  width = 6000, height = 3500, units = "px", dpi=600)
      dev.off()


      #CFI data preparation
      cf1 <- data.frame(cbind(Age, ITC))
      names(cf1) <- c("Age.plot", "CFI.plot")
      cf1 <- cf1 %>% mutate(Curve = rep("Target_CFI", length(Age)))

      cf2 <- data.frame(cbind(yout[ , 1], yout[ ,2]))
      names(cf2) <- c("Age.plot", "CFI.plot")
      cf2 <- cf2 %>% mutate(Curve = rep("Simu_CFI", length(Age)))

      cf <- rbind(cf1, cf2)
      cf3 <- cbind(Data$Age.plot, Data$CFI.plot)
      cf3 <- as.data.frame(cf3)
      names(cf3) <- c("Age.plot", "CFI.plot")

      # ggsave(file = paste0("Graphs/Step4_graphs/", Data$ANIMAL_ID, ".", "Simu_CFI",ii, ".png"), width = 6000, height = 3500, units = "px", res=600)
      cols.CFI <- c("Target_CFI" = "blue", "Simu_CFI" = "red")
      typs.CFI <- c("Target_CFI" = "solid", "Simu_CFI" = "solid")
      size.CFI <- c("Target_CFI" = 1, "Simu_CFI" = 1.8)
      plot_simu_cfi_model <- ggplot(data = cf, aes(x = Age.plot, y = CFI.plot)) +
        geom_point(data = cf3, aes(x = Age, y = CFI.obs), color = "black", shape = 21, size = 3, stroke = 0.5) +
        geom_line(aes(color = Curve, linetype = Curve, size = Curve)) +
        scale_color_manual(values = cols.CFI) +
        scale_linetype_manual(values = typs.CFI) +
        scale_size_manual(values = size.CFI) +
        geom_vline(xintercept=tbeg1, linetype="dashed",
                   color = "black", size = 0.3)+
        geom_vline(xintercept=tstop1, linetype="dashed",
                   color = "black", size = 0.3)+
        xlab("Age, d") +
        ylab("Cumulative Feed Intake, kg") +
        scale_y_continuous(breaks=seq(0, 360, 40)) +
        scale_x_continuous(breaks=seq(60, 240, 20)) +
        ggtitle(paste("Modelling CFI response of", "\nthe pig", i, "to a single perturbation")) +
        theme(plot.title = element_text(hjust = .5)) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black")
        ) +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14,face="bold"))
      ggsave(file = paste0("Graphs/Step4_graphs/Simulation CFI/", Data$ANIMAL_ID, ".", "Simu_CFI",ii, ".png"),plot_simu_cfi_model, device ="png", width = 6000, height = 3500, units = "px", dpi=600)
      # dev.off()


      #DFI data preparation
      df1 <- data.frame(cbind(Age, ITD))
      names(df1) <- c("Age.plot", "DFI.plot")
      df1 <- df1 %>% mutate(Curve = rep("Target_DFI", length(Age)))

      df2 <- data.frame(cbind(times.ode, DFI.sim))
      names(df2) <- c("Age.plot", "DFI.plot")
      df2 <- df2 %>% mutate(Curve = rep("Simu_DFI", length(times.ode)))

      df <- rbind(df1, df2)
      df3 <- cbind(Data$Age.plot, Data$CFI.plot)
      df3 <- as.data.frame(cf3)
      names(df3) <- c("Age.plot", "CFI.plot")

      # ggsave(file = paste0("Graphs/Step4_graphs/", Data$ANIMAL_ID, ".", "Simu_DFI",ii, ".png"), width = 6000, height = 3500, units = "px", res=600)
      cols.DFI <- c("Target_DFI" = "blue", "Simu_DFI" = "red")
      typs.DFI <- c("Target_DFI" = "solid", "Simu_DFI" = "solid")
      size.DFI <- c("Target_DFI" = 1.2, "Simu_DFI" = 2)
      plot_simu_dfi_model <- ggplot(data = df, aes(x = Age.plot, y = DFI.plot)) +
        geom_point(data = df3, aes(x = Age, y = DFI.obs), color = "black", shape = 21, size = 4, stroke = 1.1) +
        geom_line(aes(color = Curve, linetype = Curve, size = Curve)) +
        scale_color_manual(values = cols.DFI) +
        scale_linetype_manual(values = typs.DFI) +
        scale_size_manual(values = size.DFI) +
        geom_vline(xintercept=tbeg1, linetype="dashed",
                   color = "black")+
        geom_vline(xintercept=tstop1, linetype="dashed",
                   color = "black")+
        xlab("Age, d") +
        ylab("Daily Feed Intake, kg/ d") +
        scale_y_continuous(breaks=seq(0, 8, 1)) +
        scale_x_continuous(breaks=seq(60, 240, 20)) +
        ggtitle(paste("Modelling DFI response of", "\nthe pig", i, "to a single perturbation")) +
        theme(plot.title = element_text(hjust = .5)) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black")
        ) +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14,face="bold"))
      ggsave(file = paste0("Graphs/Step4_graphs/Simulation DFI/", Data$ANIMAL_ID, ".", "Simu_DFI",ii, ".png"),plot_simu_dfi_model, device ="png", width = 6000, height = 3500, units = "px", dpi=600)
      # dev.off()

      # #Ratio
      ratio1 <- data.frame(cbind(Time, Ratio.DFI))
      names(ratio1) <- c("Age.plot", "Ratio.plot")
      ratio1 <- ratio1 %>% mutate(Curve = rep("Compen", length(Time)))

      ratio2 <- data.frame(cbind(Time, 1+ onoff*(p1-1)))
      names(ratio2) <- c("Age.plot", "Ratio.plot")
      ratio2 <- ratio2 %>% mutate(Curve = rep("Impact", length(Time)))
      ratio <- rbind(ratio1, ratio2)

      cols.ratio  <- c("Compen" = "blue", "Impact" = "red")
      types <- c("Compen" = "solid", "Impact" = "dotdash")
      ratio_model <- ggplot(data = ratio, aes(x = Age.plot, y = Ratio.plot)) +
        geom_hline(yintercept=1, linetype="solid", color = "black", size = 1) +
        geom_line(aes(color = Curve, linetype = Curve), size = 1.2) +
        scale_color_manual(values = cols.ratio) +
        scale_linetype_manual(values = types) +
        xlab("Age (day)") +
        ylab("Ratio") +
        ggtitle(paste("Ratio between actual and target DFI", "\n of the Pig", i)) +
        theme(plot.title = element_text(hjust = .5)) +
        theme_bw()
      ggsave(file = paste0("Graphs/Step4_graphs/Ratio DFI/", Data$ANIMAL_ID, ".", "Ratio_DFI",ii, ".png"),ratio_model, device="png", width = 6000, height = 3500, units = "px", dpi=600)

      Data1 <- cbind(Age, ITC, ITD, df2$DFI.plot, cf2$CFI.plot, ratio1$Ratio.plot, ratio2$Ratio.plot, CompFI)
      names(Data1) <- c("Age","ITC", "ITD", "Sim.DFI","Sim.CFI", "Ratio.Compen", "Ratio.Impact", "CompFI")
      write.csv2(Data1,file="EAAP_Data3.csv",row.names=FALSE)

      #Legend
      LD <- rbind(AF2, AF3)

      # ggsave(file = paste0("Graphs/Step4_graphs/", Data$ANIMAL_ID, ".", "Legend",ii, ".png"), width = 6000, height = 3500, units = "px", res=600)
      cols.LD <- c("Resistance" = "black", "AA" = "blue")
      type.LD <- c("Resistance" = "solid", "AA" = "dashed")
      size.LD <- c("Resistance" = 1, "AA" = 1)
      change_model <- ggplot(data = LD, aes(x = Age.plot, y = Percent)) +
        geom_hline(yintercept=0, linetype="dashed", color = "blue", size = 1)+
        geom_line(aes(color = Curve, linetype=Curve, size = Curve)) +
        scale_color_manual(values = cols.LD, labels = c("Target_CFI", "Deviation")) +
        scale_linetype_manual(values = type.LD, labels = c("Target_CFI", "Deviation"))+
        scale_size_manual(values = size.LD, labels = c("Target_CFI", "Deviation"))+
        xlab("Age, d") +
        ylab("Change in daily feed intake (%)") +
        expand_limits(y=-100)+
        scale_y_continuous(breaks=seq(-100, 100, 10)) +
        scale_x_continuous(breaks=seq(60, 240, 20)) +
        ggtitle(paste("Change in DFI of", "the pig", i, "\nduring and after a single perturbation")) +
        theme(plot.title = element_text(hjust = .5)) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black")
        ) +
        theme(axis.text=element_text(size=12),
              axis.title=element_text(size=14,face="bold"))
      ggsave(file = paste0("Graphs/Step4_graphs/Legend/", Data$ANIMAL_ID, ".", "Legend",ii, ".png"),change_model, device="png", width = 6000, height = 3500, units = "px", dpi=600)
    }
    # dev.off()
  }
  