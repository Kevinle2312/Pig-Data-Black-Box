##################################################################              
#               JRP-INRA data - STEP 0: Data treatment
#               Masoomeh and Hieu 04/12/2017
#                     
#       Remove missing data and save as RData file
##################################################################

#Here we use the original data set of JRP
rm(list=ls(all=TRUE)) #To remove the hisory 
#dev.off() #To close all graphs


# Set working directory  

setwd("C:/Users/Kevin Le/PycharmProjects/Pig Data Black Box")

 
#===============================================================
# Import the data set and attribute a name to it
#===============================================================

  # Import data to R
  JRP_NA.0<- read.csv("Data income/RFI_JRP.csv", header =TRUE, sep=",", dec=".",fileEncoding="UTF-8-BOM")

  # Replace missing values by NA (not available)
  JRP_NA.0[JRP_NA.0 == '.'] <- NA

#===============================================================
# Remove NAs and correct type of data
#===============================================================  

  JRP_NA<- JRP_NA.0[!is.na(JRP_NA.0$FEED_INTAKE),]
  #Correct type of data
  JRP_NA$ANIMAL_ID <- as.factor(JRP_NA$ANIMAL_ID)
  JRP_NA$AGE <- as.numeric(JRP_NA$AGE)
  JRP_NA$FEED_INTAKE <- as.numeric(as.character(JRP_NA$FEED_INTAKE))

#===============================================================
#Save results to Rdata file
#===============================================================  
  
  save( JRP_NA.0, JRP_NA,                                 #Step 0: data treatment
        file = "Data income/JRPData.Rdata")

