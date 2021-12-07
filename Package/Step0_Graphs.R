##################################################################              
#               JRP-INRA data - STEP 0: Data treatment
#               Masoomeh and Hieu 04/12/2017
#                     
#       Plot the graphs of DFI and CFI
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


  # extract data for animal i 

  i <- ID[idc] # 90
  
  #extract dataset for animal i
  Data <- JRP_NA.0[JRP_NA.0$ANIMAL_ID == i,]
  Data
  #Age vector associated with 
  Age.plot <- Age[Pig_ID == i]
  
  #DFI vector associated with animal i
  DFI.plot <- DFI[Pig_ID == i]
  
  # Plot DFI graph
  
  tiff(file = paste0("Graphs/Step0_graphs/DFI/", idc, ".", Data$ANIMAL_ID, ".", "DFI", ".png"), width = 6000, height = 3500, units = "px", res=600)
  print(ggplot(Data, aes(x = Age.plot, y = DFI.plot)) +
    geom_point(shape = 21, colour = "black", fill = "white", size = 4, stroke = 1.5) + 
    xlab("Age, d") +
    ylab("Daily Feed Intake, kg/ d") +
    scale_y_continuous(breaks=seq(0, 8, 1)) +
    scale_x_continuous(breaks=seq(0, 300, 20)) +
    ggtitle(paste("Daily feed intake measured from an auto-feeder","\nPig ID:", ID[idc])) +
    theme(plot.title = element_text(hjust = .5)) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")
    ) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold")))
  png(filename = paste0("Graphs/Step0_graphs/Observasion_DFI/", idc, ".", ID[idc], ".", "DFI", ".png"),
  height = 536 , width = 700, units = 'px', type="cairo-png")
  par(mar=c(4.5,4.5,4.5,1.5))
  plot(Age.plot, DFI.plot, # Observation data of DFI
        main = paste0("Daily Feed Intake measured from automatic feeder", "\nPig ID = ", i, ",", " idc = ", idc),
        ylab = "Daily Feed Intake (kg)",
        xlab = "Age (days)",
        type="p", col ="black", cex.main = 1.7, cex = 1.5,
        cex.axis = 1.5, cex.lab = 1.5)
  legend("topleft", "DFI (kg/ day)",
         col="black",
         pch=1, bty = "n", cex = 1.2)
  
  dev.off()
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
  tiff(file = paste0("Pig Data Black Box/Graphs/Step0_graphs/CFI/", idc, ".", Data$ANIMAL_ID, ".", "CFI", ".png"), width = 6000, height = 3500, units = "px", res=600)
  print(ggplot(Data, aes(x = Age.plot, y = CFI.plot)) +
    geom_point(shape = 21, colour = "black", fill = "white", size = 3, stroke = 1.2) + 
    xlab("Age, d") +
    ylab("Cumulative Feed Intake, kg") +
    scale_y_continuous(breaks=seq(0, 400, 40)) +
    scale_x_continuous(breaks=seq(0, 300, 20)) +
    ggtitle(paste("Cumulative Feed Intake","\nPig ID:", ID[idc])) +
    theme(plot.title = element_text(hjust = .5)) + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")
    ) +
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=14,face="bold")))
  dev.off()

  
} # end of FOR loop
