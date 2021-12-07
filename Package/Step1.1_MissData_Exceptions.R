##################################################################              
#               JRP data - Missing data
#               Masoomeh and Hieu 04/12/2017
#               Automatic detection of missing rows
#               
##################################################################

  rm(list=ls(all=TRUE)) #To remove the hisory 
  dev.off() #To close all graphs

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
  
  load("Data income/JRPData.Rdata")
  


  #===============================================================
  # DATA PREPARATION
  #===============================================================

  #Create a new dataframe which stores Data after treatment of missing data
  No.NA.Data.0 <- data.frame(AnimalID=factor(),
                             AGE=double(),
                             DFI=double(),
                             CFI=double(),
                             stringsAsFactors=FALSE)
   
  #----------------------------------------------------
  # Attribute new easier name to each column
  #----------------------------------------------------
    
  #Order number of Animal_ID
  ID <- unique(JRP_NA$ANIMAL_ID) 
    
  #Pig ID  
  Pig_ID <- JRP_NA$ANIMAL_ID
    
  #Age (days)
  Age <- JRP_NA$AGE
    
  #Daily Feed Intake (kg)
  DFI <- JRP_NA$FEED_INTAKE
    
  #===============================================================
  # For loop for automatically filtering missing data for all pigs
  #===============================================================

  # IDC = c(1:2,4:51, 53:79, 81:119) # animals 3, 52 and 80 are removed because of some problems

   for(idc in seq_along(ID)){
    # idc = 2
      #----------------------------------------------------
      # extract data for animal i
      #----------------------------------------------------
      i <- ID[idc] # 90

      #extract dataset for animal i
      Data <- JRP_NA[JRP_NA$ANIMAL_ID == i,]

      #Age vector associated with
      Age.plot <- Age[Pig_ID == i]

      #DFI vector associated with animal i
      DFI.plot <- DFI[Pig_ID == i]

      # Plot DFI graph

      # png(filename = paste("Graphs/Step0_graphs/",idc,".",ID[idc],".","DFI",".png",sep=""),
      # height = 536 , width = 700, units = 'px', type="cairo-png")
      #par(mar=c(4.5,4.5,4.5,1.5))
      #plot(Age.plot, DFI.plot, # Observation data of DFI
           # main = paste0("Daily Feed Intake measured from automatic feeder", "\nPig ID = ", i, ",", " idc = ", idc),
           # ylab = "Daily Feed Intake (kg)",
           # xlab = "Age (days)",
           # type="p", col ="black", cex.main = 1.7, cex = 1.5,
           # cex.axis = 1.5, cex.lab = 1.5)
      # legend("topleft", "DFI (kg/ day)",
      #        col="black",
      #        pch=1, bty = "n", cex = 1.2)

      #dev.off()
      #----------------------------------------------------
      # Calculate CFI from DFI
      #----------------------------------------------------

      #CFI
      CFI.plot <- NULL
      CFI.plot[1] <- DFI.plot[1]
                for(j in 2:length(Age.plot)){
                  CFI.plot[j] <- CFI.plot[j-1]+DFI.plot[j]
                }

      # Plot CFI graph
      # png(filename = paste("Graphs/Step0_graphs/",idc,".",ID[idc],".","CFI",".png",sep=""),
      #     height = 536 , width = 700, units = 'px', type="cairo-png")
      #par(mar=c(4.5,4.5,4.5,1.5))
      #plot(Age.plot, CFI.plot, # Observation data of DFI
           # main = paste0("Cumulative Feed Intake measured from automatic feeder", "\nPig ID = ", i, ",", " idc = ", idc),
           # ylab = "Cumulative Feed Intake (kg)",
           # xlab = "Age (days)",
           # type="p", col ="black", cex.main = 1.7, cex = 1.5,
           # cex.axis = 1.5, cex.lab = 1.5)
      # legend("topleft", "DFI (kg/ day)",
      #        col="black",
      #        pch=1, bty = "n", cex = 1.2)

      #dev.off()
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
 #          FILTRATING MISSING DATA
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
            # Extract age and CFI before and after missing data
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

            #----------------------------------------------------
            # Exceptional conditions
            #----------------------------------------------------

            #1. When number of data before missing row is not enough for estimation
            # Ex. 88, 89,..., 92, 93 we remove all the rows until missing data (Age now starts at day 92)
            if(length(MissAgeT)== 1 && b1l[[k]][1] < Age.plot[1]){
              JRP_new <- JRP_new[!JRP_new$Age.plot %in% seq(Age.plot[1], MissAge2[length(MissAge2)], by=1),]
              MissAgeT <- list()
              MissAge <- MissAgeT
              MissAge2 <- MissAge
            }

            # Ex. 88, 89,..., 92, 93 and 100, 101,..., 103, 104 we remove all the rows until missing data (Age now starts at day 92)
            #And the first missing value will be 102
            if(length(MissAgeT)>1 && b1l[[k]][1] < Age.plot[1]){
              ToPig_new <- ToPig_new[!ToPig_new$Age.plot %in% seq(Age.plot[1], MissAge2[length(MissAge2)], by=1),]
              MissAgeT[[1]] <- NULL
              MissAge <- MissAgeT[[k]]
              MissAge2 <- MissAge
            }

            #2. when there are not enough given rows between two series of missing rows
            #    Ex. 108, 109, 110, 112, 114, 115, 116
            #    We merge these two missing rows with the rows between
            # when number of missing row in a serie > 1

            if(length(MissAgeT)>1 && k < length(MissAgeT)){   #pig i[91]
              for(ii in 1:(length(MissAgeT)-1)){
                if(k + ii <= length(MissAgeT)){
                  if(MissAgeT[[k+ii]][1] - MissAge[length(MissAge)] <= length(MissAge)+1){
                    (MissAge <- c(MissAge, seq(MissAge[length(MissAge)]+1, MissAgeT[[k+ii]][1]-1, 1), MissAgeT[[k+ii]]))&&
                      (b2l[[k]] <- seq( MissAge[length(MissAge)]+1,MissAge[length(MissAge)] + (length(MissAge)+1), 1)) &&
                      (b1l[[k]] <- seq(MissAge[1]-(length(MissAge)+1), MissAge[1]-1, 1))&&
                      (ki <- k+ii) &&
                      (MissAge2 <- MissAge)
                  }
                }
              }
            }

            #3. When number of data after missing row is not enough for estimation
            # Ex. 179, 180,..., 183, 184 we remove all the rows from missing data (Age is now just until day 180)
            if(b2l[[k]][length(b2l[[k]])] > Age.plot[length(Age.plot)]){
              JRP_new <- JRP_new[!JRP_new$Age.plot %in% seq(MissAge2[1],Age.plot[length(Age.plot)], by=1),]
              MissAgeT[[length(MissAgeT)]] <- NULL
              break
            }


      #----------------------------------------------------
      # Re-calculate CFI after filtrating process
      #----------------------------------------------------

      #Age vector associated with
      Age.plot <- JRP_new$Age.plot[! JRP_new$Age.plot %in% MissAge2]

      #DFI vector associated with animal i
      DFI.plot <- JRP_new$DFI.plot[JRP_new$Age.plot %in% Age.plot]

      #CFI
      CFI.plot <- NULL
      CFI.plot[1] <- DFI.plot[1]
      for(j in 2:length(Age.plot)){
        CFI.plot[j] <- CFI.plot[j-1]+DFI.plot[j]  }

      #----------------------------------------------------
      # New dataset with removed NA only containing Age, DFI and CFI
      #----------------------------------------------------

      JRP_new <- as.data.frame(cbind(Age.plot,DFI.plot, CFI.plot))


      ki <- ki+1 #next series of missing rows
      #JRP_new with no missing rows

      } #end of WHILE loop

      #-------------------------------------------------------------------------------
      #  Save data with missing data estimated to dataframe No.NA.Data.0
      #-------------------------------------------------------------------------------

      ANIMAL_ID <- rep(ID[idc], dim(JRP_new)[1])
      JRP_new_Final <- cbind (ANIMAL_ID,JRP_new)
      No.NA.Data.0 <- rbind(No.NA.Data.0 , JRP_new_Final )

  } # end of FOR loop

    
  #===============================================================
  #Save results to Rdata file
  #===============================================================  

      save( JRP_NA.0, JRP_NA,              #Step 0: data treatment
            No.NA.Data.0,                        #Step 1: missing data
            file = "Data income/JRPData.Rdata")

    