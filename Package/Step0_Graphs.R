  ##################################################################
  #               JRP-INRA data - STEP 0: Data treatment
  #               Masoomeh and Hieu 04/12/2017
  #
  #       Plot the graphs of DFI and CFI
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
  # library(tidyverse)
  library(lubridate)

  ################################################################

  #Set working directory

  setwd("C:/Users/Kevin Le/PycharmProjects/Pig Data Black Box")

  #load data after data treatment

  load("Data/JRPData.Rdata")



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
    Data <- as.data.frame(Data)
    #Age vector associated with
    Age.plot <- Age[Pig_ID == i]

    #DFI vector associated with animal i
    DFI.plot <- DFI[Pig_ID == i]

    # Plot DFI graph

    # tiff(file = paste0("Graphs/Step0_graphs/DFI/", idc, ".", Data$ANIMAL_ID, ".", "DFI", ".png"), width = 6000, height = 3500, units = "px", res=600)
    DFI.fig <- plot_ly(data.frame(Data), x = ~Age.plot, y = ~DFI.plot,
               marker = list(size = 10,
                             color = "#052CA3",
                             line = list(color = 'rgba(152, 0, 0, .8)',
                                         width = 2))) %>%
                layout(title =paste("Daily Feed Intake","\nPig ID:", ID[idc]), plot_bgcolor = "fffff", xaxis = list(title = 'Age (d)'),
                        yaxis = list(title = 'Daily Feed Intake, kg') )
                scale_y_continuous(breaks=seq(0, 400, 40)) %>%
                scale_x_continuous(breaks=seq(0, 300, 20)) %>%
                theme(plot.title = element_text(hjust = .5)) %>%
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black")
                ) %>%
                theme(axis.text=element_text(size=12),
                      axis.title=element_text(size=14,face="bold"))
    saveWidget(ggplotly(DFI.fig), file=paste0("C:/Users/Kevin Le/PycharmProjects/Pig Data Black Box/Graphs/Step0_graphs/Observasion_DFI/", idc, ".", ID[idc], ".", "DFI", ".html"))

    #----------------------------------------------------
    # Calculate CFI from DFI
    #----------------------------------------------------

    #CFI
    CFI.plot <- NULL
    CFI.plot[1] <- DFI.plot[1]
    for(j in 2:length(Age.plot)){
      CFI.plot[j] <- CFI.plot[j-1]+DFI.plot[j]
    }

    dt <- cbind(Age.plot,CFI.plot)

    # Plot CFI graph
    # tiff(file = paste0("Pig Data Black Box/Graphs/Step0_graphs/CFI/", idc, ".", Data$ANIMAL_ID, ".", "CFI", ".png"), width = 6000, height = 3500, units = "px", res=600)
    CFI.Fig <- plot_ly(data.frame(dt), x = ~Age.plot, y = ~CFI.plot,
               marker = list(size = 10,
                             color = "#052CA3",
                             line = list(color = 'rgba(152, 0, 0, .8)',
                                         width = 2))) %>%
                layout(title =paste("Cumulative Feed Intake","\nPig ID:", ID[idc]), plot_bgcolor = "fffff", xaxis = list(title = 'Sepal Length (cm)'),
                        yaxis = list(title = 'Cumulative Feed Intake, kg') )
                scale_y_continuous(breaks=seq(0, 400, 40)) %>%
                scale_x_continuous(breaks=seq(0, 300, 20)) %>%
                theme(plot.title = element_text(hjust = .5)) %>%
                theme(panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black")
                ) %>%
                theme(axis.text=element_text(size=12),
                      axis.title=element_text(size=14,face="bold"))
    # CFI.Fig <- CFI.Fig %>%  layout(title =paste("Cumulative Feed Intake","\nPig ID:", ID[idc]), plot_bgcolor = "#e5ecf6", xaxis = list(title = 'Sepal Length (cm)'),
    #      yaxis = list(title = 'Cumulative Feed Intake, kg') )





    # CFI.Fig <- plot_ly(data.frame(dt),mapping = aes(x = Age.plot , y = CFI.plot)) %>%
    #   geom_point(shape = 21, colour = "black", fill = "white", size = 3, stroke = 1.2) %>%
    #   layout(title =paste("Cumulative Feed Intake","\nPig ID:", ID[idc]), plot_bgcolor = "#e5ecf6", xaxis = list(title = 'Sepal Length (cm)'),
    #      yaxis = list(title = 'Cumulative Feed Intake, kg') )

    saveWidget(ggplotly(CFI.Fig), file=paste0("C:/Users/Kevin Le/PycharmProjects/Pig Data Black Box/Graphs/Step0_graphs/Observasion_CFI/", idc, ".", ID[idc], ".", "CFI", ".html"))
    # dev.off()


  } # end of FOR loop

