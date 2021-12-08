##################################################################              
#               JRP data - Missing data
#               Masoomeh and Hieu 04/12/2017
#               Automatic detection of missing rows
#               
##################################################################

rm(list=ls(all=TRUE)) #To remove the hisory
 # dev.off() #To close all graphs

#-------------------------------------------------------------------------------
# Packages needed for estimaton of Ideal trajectory - nonlinear regression
#-------------------------------------------------------------------------------

  library(GrapheR) 
  library(lattice)
  library(ggplot2)
  library(plotrix)

################################################################  

#Set working directory  

setwd("C:/Users/Kevin Le/PycharmProjects/Pig Data Black Box")

#load data after data treatment 

  load("Data/JRPData.Rdata")

  #===============================================================
  # DATA PREPARATION
  #===============================================================
  
  #Create a new dataframe which stores Data after treatment of missing data
    No.NA.Data <- data.frame(AnimalID=factor(),
                             AGE=double(),
                             DFI=double(),
                             CFI=double(),
                             stringsAsFactors=FALSE)
  
  #Order number of Animal_ID
   ID <- unique(as.factor(No.NA.Data.0$ANIMAL_ID))

  #===============================================================
  # For loop for automatically estimating and distributing missing data of all pigs
  #===============================================================
   
  for(idc in seq_along(ID)){
  
    #----------------------------------------------------
    # extract data for animal i 
    #----------------------------------------------------
    
    i <- ID[idc]
    
    #extract dataset for animal i
    Data <- No.NA.Data.0[No.NA.Data.0$ANIMAL_ID == i,]
    
    #Age vector associated with 
    Age.plot <- Data$Age.plot
    
    #DFI vector associated with animal i
    DFI.plot <- Data$DFI.plot
    
    #CFI vector associated with animal i
    CFI.plot <- Data$CFI.plot
    
    
    #----------------------------------------------------
    # New dataset with removed NA only containing Age, DFI and CFI
    #----------------------------------------------------
    
    JRP_new <- as.data.frame(cbind(Age.plot,DFI.plot, CFI.plot)) 
    
    # to keep a backup of this dataset we rename it
    JRP_new_initial <-JRP_new 
    
  #====================================================
  #           DETECTING MISSING DATA
  #====================================================
    
    #----------------------------------------------------
    # 1. Missing lines detection
    #----------------------------------------------------
    
    #Expected number of rows
    Exp.row <- max(Age.plot) - min(Age.plot)+1 
    
    #observed number of rows
    obs.row <- as.numeric(length(Age.plot))
    
    #number of missing days
    NumMissRow <-  Exp.row - obs.row 
    
    #To locate the missing rows 
    #if age of row_(n+1) - row_n > 1 then miss row
    A1 <- Age.plot[seq_along(Age.plot) -1] #
    A2 <- Age.plot[2:length(Age.plot)]
    A3 <- A2 - A1 
    
    #Day before missing series of data
    MissRow <- as.data.frame(JRP_new[A3!=1,])
    Missdays <- MissRow$Age.plot 
    
    # Extract A3 values different to 1
    A4 <- A3[A3!=1]
    
    #----------------------------------------------------
    # 2.List of ages associated with missing rows
    #----------------------------------------------------
    
    #An empty object (here a  list) that will gather 
    # all missing rows for a given pig will be stored in a list
    MissAgeT <- list() 
    
    #JRP_new with no missing rows
    #This dataframe will be used to compare with JRP_new_Final if needed
    JRP_new_Final <-JRP_new  
    
    # Missing ages in the dataframe
    if (length(A4)>0){
      for (ii in seq_along(A4)){
        MissAgeT[[ii]] <- seq(Missdays[ii]+1,Missdays[ii]+A4[ii]-1,1) 
      }
    }
    MissAgeT
    
  #====================================================
  #          ESTIMATING MISSING DATA 
  #====================================================
    
   #before going to "while" or "for " loop, we create empty lists
   #in which we will store all data for a given animal
    
    b1l <- list() #days before missing rows
    b2l <- list() #days after missing rows
    
    Xl <- list()   #list of ages before and after missing rows
    Yl <- list()   #list of CFI before and after missing rows
    tal <- NULL      #the vector of first missing rows
    dl <- list()   #list of coefficients for delta
    
    #information of linear model will be put in list associeted with that information 
    Par_initl <- list() #Initial parameters for each set of missing rows
    resl <- list() #function optim to estimate parameters
    pred_valuesl <- list()
    deltal <- NULL
    P.dfil <- list() #parameters of the linear function of DFI
    FI.missl <- list()
    Corr.datal <- list()
    Res.RSS <- list()
    
    #----------------------------------------------------
    # 3. While loop
    #----------------------------------------------------
    
    ki <- 1 #first value to start while loop
    
    while(ki< length(MissAgeT)+1 ){
      
      k <- ki
      MissAge <- MissAgeT[[k]] #Missing rows at position k
      
    #----------------------------------------------------
    # 3.1. Extract age and CFI before and after missing data
    #      for the interpolation via linear regression
    #----------------------------------------------------
      
      L <- as.numeric(length(MissAge))
      
      #number of data before = length(missing rows)+1
      b1l[[k]] <- seq(MissAge[1]-(length(MissAge)+1), MissAge[1]-1, 1)
      
      #number of data after missing age
      b2l[[k]] <- seq( MissAge[length(MissAge)]+1,
                       MissAge[length(MissAge)]+
                         (length(MissAge)+1), 1)
      
      #number of missing rows to be replaced at the end
      MissAge2 <- MissAge
      
      b1 <- b1l[[k]] # Ages before missing row(s)
      b2 <- b2l[[k]] # Ages after missing row(s)
      
      # rows associated with these ages b1, b2
      B1 <- JRP_new[JRP_new$Age.plot %in% b1,]
      B2 <- JRP_new[JRP_new$Age.plot %in% b2,]
      
      #CFI associated with these ages b1, b2  
      CFI_b <- B1$CFI.plot
      CFI_a <- B2$CFI.plot
      
      #DFI associated with these ages b1, b2 
      DFI_b <-B1$DFI.plot
      DFI_a <-B2$DFI.plot
      
    #----------------------------------------------------
    #          Regression 
    #----------------------------------------------------
      
      # Put all values of age together
      Xl[[k]] <-c(b1,b2) 
      Yl[[k]] <-c(CFI_b,CFI_a)
      
      tal[k] <- MissAge[length(MissAge)] #from this point delta is valid
      
      X <- Xl[[k]]
      Y <- Yl[[k]]
      ta <- tal[k]
      
      dl[[k]]  <- as.numeric(X>ta)
      d <- dl[[k]] #on/off condition which activates function after missing rows
      
      #################################
      
      #Quadratic function of estimating missing CFI values
      pred_fn <- function(P){
        a0 <- P[1]
        a1 <- P[2]
        a2 <- P[3]
        b0 <- P[4]
        return ( a2*X^2  +  a1*X + a0*(1-d)+b0*d)
      }

      # 2 days before and after for calculation of initial values
      
      X2 <- b1[length(b1)-1] #X[2]
      X3 <- b1[length(b1)] #X[3]
      X4 <- b2[1] #X[4]
      X5 <- b2[2]#X[5]
      
      Y2 <- CFI_b[length(CFI_b)-1]#Y[2]
      Y3 <- CFI_b[length(CFI_b)]#Y[3]
      Y4 <- CFI_a[1] #Y[4]
      Y5 <- CFI_a[2] #Y[5]
      
      #Reparametrize parameters of quadratic function to parameters from 2 days
      #before and after missing rows
      
      a2_init <- -(-Y3*X5+Y3*X4+Y2*X5-Y2*X4+X2*Y4-X2*Y5-X3*Y4+X3*Y5)/
        (X2*X5^2-X2*X4^2-X2^2*X5+X2^2*X4-X3*X5^2+X3*X4^2+X3^2*X5-X3^2*X4)


      a1_init <- (-X5^2*Y3+X5^2*Y2-Y5*X2^2+Y5*X3^2+Y4*X2^2-X4^2*Y2-Y4*X3^2+X4^2*Y3)/
        ((X5-X4)*(X5*X2-X5*X3-X2^2+X3^2+X4*X2-X4*X3))

      a0_init <- -(-X5^2*Y3*X2+X5^2*Y2*X3+X5*Y3*X2^2-
                    X5*Y2*X3^2+Y4*X2^2*X3-Y3*X2^2*X4-
                    X2^2*X3*Y5+X2*Y5*X3^2-Y4*X3^2*X2+Y3*X2*X4^2-Y2*X3*X4^2+Y2*X3^2*X4)/
        (X2*X5^2-X2*X4^2-X2^2*X5+X2^2*X4-X3*X5^2+X3*X4^2+X3^2*X5-X3^2*X4)


      b0_init <- -(-Y4*X2*X5^2+Y4*X3*X5^2+X5^2*Y2*X4-
                    X5^2*Y3*X4+Y4*X2^2*X5-Y4*X3^2*X5-
                    X4^2*Y2*X5+X4^2*Y3*X5-Y5*X2^2*X4+Y5*X2*X4^2-Y5*X3*X4^2+Y5*X3^2*X4)/
        (X2*X5^2-X2*X4^2-X2^2*X5+X2^2*X4-X3*X5^2+X3*X4^2+X3^2*X5-X3^2*X4)

      # initial values calculated reparametrization
      Par_initl[[k]] <- c(a0_init, a1_init, a2_init, b0_init) #initial values of params
      Par_init <- Par_initl[[k]]
      
      
      # Using optim function to minimize Residual Sum of Squares
      # between predicted and observed CFI to estimate function parameters
      
      obj <- function(P){
        output   <- (Y -pred_fn(P))
        obj.func <- t(output)%*%output # Residual sum of squares
        return(obj.func)
      }
      
      resl[[k]] <- optim (Par_init, obj) #function optim to estimate parameters
      res <- resl[[k]]
      
      pred_valuesl[[k]] <- pred_fn(res$par)
      pred_values <- pred_valuesl[[k]]
      
      
      #parameters values
      a0 <- res$par[1]
      a1 <- res$par[2]
      a2 <- res$par[3]
      b0 <- res$par[4]
      deltal[k] <- a0 - b0
      delta <- deltal[k]
      
      # Plot the function
      P1 <- c(a0, a1, a2, b0)
      plot(X, Y)
      lines(X, pred_fn(P1))
      lines(X, a2*X^2+a1*X+a0, col = "red")
      lines(X, a2*X^2+a1*X+b0, col = "blue")
      
    #################################################
      
    #----------------------------------------------------
    #         Attribution of data to missing age 
    #----------------------------------------------------
      
      #Function of DFI (1st derivative of quadratic fucntion)
      DFI_pred <- function(age,P.dfi){
        return (2*P.dfi[2]*age  +  P.dfi[1])
      }

      #parameters of the linear function of DFI (2*a2*age + a1)
      P.dfil[[k]] <- c(a1,a2)
      P.dfi <- P.dfil[[k]]
      
      #sum of estimated DFI for missing days
      CFI_miss <- sum(DFI_pred(MissAge2, P.dfi))
      
      FI.miss <- NULL
      for (kj in seq_along(MissAge2)){
        FI.miss[kj] <- DFI_pred(MissAge2[kj],P.dfi)*delta /CFI_miss}
      FI.missl[[k]] <- FI.miss
      
      #Put estimation of DFI in right places (attributed rows)
      t1 <- as.numeric(JRP_new$Age.plot %in% MissAge2-1)
      probs <- rep(TRUE, length(DFI.plot[!JRP_new$Age.plot %in% MissAge2]))
      ind <- t1[1] + MissAge2/(sum(MissAge2)+.05)
      val.FI <- c(DFI.plot[!JRP_new$Age.plot %in% MissAge2], FI.miss)
      val.age <- c(Age.plot[!JRP_new$Age.plot %in% MissAge2], MissAge2)
      
      # val.FI <- c( DFI.plot, FI.miss)
      # val.age <- c( Age.plot[!c(JRP_new$Age.plot %in% MissAge)],MissAge)
      id  <- c( seq_along(probs), ind)
      
      Corr.DFI <- val.FI[order(id)]
      Corr.Age <- val.age[order(id)]
      Corr.Age <- Corr.Age[order(Corr.Age)]
    #----------------------------------------------------
    #     Recalculate CFI based on new DFI 
    #----------------------------------------------------
      
      Corr.CFI <- NULL
      Corr.CFI[1] <- Corr.DFI[1]
      
      for(j in 2:length(Corr.Age)){
        Corr.CFI[j] <- Corr.CFI[j-1]+Corr.DFI[j]
      }
      
      #New dataframe containing estimated missing data
      Corr.datal[[k]] <- as.data.frame(cbind(Corr.Age, Corr.DFI, Corr.CFI))
      Corr.data <- Corr.datal[[k]]
      colnames(Corr.data) <- c("Age.plot", "DFI.plot", "CFI.plot")
      
    ######################################################################
      
      #RSS and delta
      Res.RSS[[k]] <- as.data.frame(c(res$value,delta) , row.names = c("RSS", "Delta"))
      colnames(Res.RSS) <- "Results"
      # 
      
    ######################################################################
      
    #----------------------------------------------------
    #We replace new data set in the initial data set and we restart 
    #the processus while the condition for ki is respected  
    #----------------------------------------------------
      
      JRP_new <- Corr.data
      DFI.plot <- JRP_new$DFI.plot
      Age.plot <- JRP_new$Age.plot
      CFI.plot <- JRP_new$CFI.plot
      
    #----------------------------------------------------
    #          Plot the estimation 
    #----------------------------------------------------     
      
      #Vector contains estimated missing data
      X1 <- c(b1, MissAge2, b2)
      Y1 <- JRP_new$CFI.plot[JRP_new$Age.plot %in% X1] #CFI data
      
      #Initial Age around missing row
      Age.ini <- JRP_new_initial$Age.plot[JRP_new_initial$Age.plot >= b1[1] & JRP_new_initial$Age.plot <= b2[length(b2)]]
      Z <- JRP_new_initial$DFI.plot[JRP_new_initial$Age.plot %in% Age.ini] #Initial DFI values
      Z1 <- JRP_new$DFI.plot[JRP_new$Age.plot %in% X1] #Estimated DFI values
      
      #To indicate the location of missing rows in the title of graphs
      A <- NULL
      for(i in seq_along(MissAge2)){
        A <- paste(A, MissAge2[i])
      }
      
 #====================================================
 #     Plot CFI graph after missing data estimation
 #====================================================    
      
      #To save plots  you have to create a folder called "Step1_graphs"
      #in the same working directory with this R-script     
      
      # For CFI plot
      # png(filename = paste("Graphs/Step1_graphs/",idc,".",ID[idc],".","ki=", ki,".","CFI_MA_esti",".png",sep=""),
      #     height = 720, width = 1200, units = 'px', type="cairo-png")
      # plot(X1, Y1,
      #      main = paste0("Pig_ID = ", ID[idc], ",", " idc = ", idc, "\nMissing rows:", A),
      #      ylab = "Cummulative feed intake (kg)",
      #      xlab = "Age (days)",
      #      type = "o", lwd = 2, cex.main = 1.7,
      #       pch = 1, cex.lab = 1.5, cex.axis = 1.3 , cex = 2, col = 'blue')
      # points(X, Y,
      #        lwd = 2, cex = 2,  type = "o")
      # legend("bottomright",c("real data", "estimated data"),
      #        pch = c(1,1), col = c("black", "blue"), bty = "n", cex = 1)
      # dev.off()
      
      # For DFI
      # png(filename = paste("Graphs/Step1_graphs/",idc,".",ID[idc],".","ki=", ki,".","DFI_MA_esti",".png",sep=""),
      #     height = 720, width = 1200, units = 'px', type="cairo-png")
      # plot(X1, Z1,
      #      main = paste0("Pig_ID = ", ID[idc], ",", " idc = ", idc, "\nMissing rows:", A),
      #      ylab = "Daily feed intake (kg)",
      #      xlab = "Age (days)",
      #      type = "o", lwd = 2, cex.main = 1.7,
      #      pch = 1, cex.lab = 1.5, cex.axis = 1.3 , cex = 2, col = 'blue')
      # points(Age.ini, Z,
      #        lwd = 2, cex = 2,  type = "o")
      # legend("bottomright",c("real data", "estimated data"),
      #        pch = c(1,1), col = c("black", "blue"), bty = "n", cex = 1)
      # dev.off()
      
      ki <- ki+1 #next series of missing rows
      #JRP_new with no missing rows
      
    } #end of WHILE loop
    
    #==========================================================================
    #  Remove the last day because sometimes animal is fasting before slaughter 
    #==========================================================================
    
    JRP_new <- JRP_new[seq_along(JRP_new$Age.plot) -1,]
  
    #----------------------------------------------------------------
    #  Save data with missing data estimated to dataframe No.NA.Data 
    #----------------------------------------------------------------
    
    ANIMAL_ID <- rep(ID[idc], dim(JRP_new)[1])
    JRP_new_Final <- cbind (ANIMAL_ID,JRP_new)  
    No.NA.Data <- rbind(No.NA.Data , JRP_new_Final )
    
    #====================================================
    #     Plot CFI graph after missing data estimation
    #====================================================
    
    png(filename = paste0("Graphs/Step1_graphs/", idc, ".", ID[idc], ".", "MA_esti_overal", ".png"),
         height = 720, width = 1200, units = 'px', type="cairo-png")
    plot(JRP_new_Final$Age.plot, JRP_new_Final$CFI.plot, 
         main = paste0("Pig_ID = ", ID[idc], ",", " idc = ", idc, "\nEstimating missing data"),
         ylab = "Cummulative feed intake (kg)", 
         xlab = "Age (days)",
         cex.main = 1.7 , pch = 1, cex.lab = 1.5, cex.axis = 1.3 , lwd = 2,  
         col = 'blue', type = "l")
    points(JRP_new_initial$Age.plot, JRP_new_initial$CFI.plot,
           pch = 1, cex = 1.5, col = 'black')
    legend("topleft",c("real data", "estimated data"), 
           col = c("black","blue"), pch = c(1,NA), lwd = c(NA, 1), bty = "n")
    dev.off()
    
    
  } # end of FOR loop

#===============================================================
#Save results to Rdata file
#===============================================================  

save( JRP_NA.0, JRP_NA,              #Step 0: data treatment
      No.NA.Data.0, No.NA.Data,          #Step 1: missing data
      file = "Data/JRPData.Rdata")

   