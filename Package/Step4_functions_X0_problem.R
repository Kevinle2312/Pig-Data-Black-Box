
  #==========================================================================
  # Function of perturbation
  # Masoomeh Jaap Hieu
  # 15/11/2017
  #==========================================================================


                        #============================ ======== ========
                        # New function discussion Jaap and Hieu 14/11/2017
                        #============================ ======== ========

  #=======================================================================================
  # 1. function for perturbation model  -- Model ode.pert.14/11 - Jaap modified 03/08/2018
  #=======================================================================================

  #------------------------------
  # Linear function for TTC
  #------------------------------

  #Fitting initial parameters to the model
  ODE.CFI.0 <- function(Time, State, Pars) {
  #
  p1 <- Pars[1]  #percentage of decrease of ITD during the 1st perturbation
  max.compFI1 <- Pars[2]
  tbeg1 <- Pars[3]
  tstop1 <- Pars[4]
  a <- Pars[5]
  b <- Pars[6]
  
  
  #
  with(as.list(c(State, Pars)), {
    # 
    x <- Time
    #
    ITC.CFI <- a +b*x  #returns the cubic function
    ITC.DFI <- b
    #
    onoff <- ifelse(Time>tbeg1 & Time<tstop1, 1, 0)
    CompFI <- (1-CumFI/ITC.CFI)*max.compFI1
    #
    dCumFI <- (onoff*(p1-1) + CompFI + 1)*ITC.DFI
    #
    return(list(dCumFI))
  })
}

  #Function to estimate parameters
  ODE.CFI.optim.0 <- function(Time, State, Pars) {
  #
  p1 <- Pars[1]  #percentage of decrease of ITD during the 1st perturbation
  max.compFI1 <- Pars[2]
  # tbeg1 <- Pars[3]
  # tstop1 <- Pars[4]
  # a <- Pars[5]
  # b <- Pars[6]
  # CumFI <- Pars[7]
  # c           = Pars[7]
  # d           = Pars[8]
  
  #
  with(as.list(c(State, Pars)), {
    # 
    x <- Time
    #
    ITC.CFI <- a + b*x  #returns the cubic function
    ITC.DFI <- b
    #
    onoff <- ifelse(Time>tbeg1 & Time<tstop1, 1, 0)
    CompFI <- (1-CumFI/ITC.CFI)*max.compFI1
    #
    dCumFI <- (onoff*(p1-1) + CompFI + 1)*ITC.DFI
    #
    return(list(dCumFI))
  })
}

  #Objective function
  ODE.CFI.obj.0 <- function(P, data){
  data <- Data.xy
  f <- ode(yinit, times, ODE.CFI.optim.0, P)
  optim.f <- sqrt(sum((data$CFI.plot- f[, 2])^2))
  return (optim.f)
}


  #nls function for initial parameter
  ODE.CFI.obj.nls.0 <- function(p1,p2){
    P <- c(p1,p2)
    data <- Data.xy
    f <- ode(yinit, times, ODE.CFI.optim.0, P)
    return ( f[,2])
  }

  #------------------------------
  # Quadratic function for TTC
  #------------------------------

  #Fitting initial parameters to the model
  ODE.CFI.1 <- function(Time, State, Pars) {
    #
    p1 <- Pars[1]  #percentage of decrease of ITD during the 1st perturbation
    max.compFI1 <- Pars[2]
    tbeg1 <- Pars[3]
    tstop1 <- Pars[4]
    a <- Pars[5]
    b <- Pars[6]
    c <- Pars[7]

    #
    with(as.list(c(State, Pars)), {
      #
      x <- Time
      #
      ITC.CFI <- a + b*x + c*x^2   #returns the cubic function
      ITC.DFI <- b+2*c*x
      #
      onoff <- ifelse(Time>tbeg1 & Time<tstop1, 1, 0)
      CompFI <- (1-CumFI/ITC.CFI)*max.compFI1
      #
      dCumFI <- (onoff*(p1-1) + CompFI + 1)*ITC.DFI
      #
      return(list(c(dCumFI)))
    })
  }

  #Function to estimate parameters
  ODE.CFI.optim.1 <- function(Time, State, Pars) {
    #
    p1 <- Pars[1]  #percentage of decrease of ITD during the 1st perturbation
    max.compFI1 <- Pars[2]
    # tbeg1 <- Pars[3]
    # tstop1 <- Pars[4]
    # a           = Pars[5]
    # b           = Pars[6]
    # c           = Pars[7]
    # d           = Pars[8]

    #
    with(as.list(c(State, Pars)), {
      #
      x <- Time
      #
      ITC.CFI <- a + b*x + c*x^2  #returns the cubic function
      ITC.DFI <- b+2*c*x
      #
      onoff <- ifelse(Time>tbeg1 & Time<tstop1, 1, 0)
      CompFI <- (1-CumFI/ITC.CFI)*max.compFI1
      #
      dCumFI <- (onoff*(p1-1) + CompFI +1)*ITC.DFI
      #
      return(list(dCumFI))
    })
  }

  #Objective function
  ODE.CFI.obj.1 <- function(P, data){
    data <- Data.xy
    f <- ode(yinit, times, ODE.CFI.optim.1, P)
    optim.f <- sqrt(sum((data$CFI.plot- f[, 2])^2))
    return (optim.f)
  }

  #nls function for initial parameter
  ODE.CFI.obj.nls.1 <- function(p1,p2){
    P <- c(p1,p2)
    data <- Data.xy
    f <- ode(yinit, times, ODE.CFI.optim.1, P)
    return ( f[,2])
  }

  #----------------------------------
  # Quadratic-linear function for TTC
  #----------------------------------

  #Fitting initial parameters to the model
  ODE.CFI.2 <- function(Time, State, Pars) {
    #
    p1 <- Pars[1]  #percentage of decrease of ITD during the 1st perturbation
    max.compFI1 <- Pars[2]
    tbeg1 <- Pars[3]
    tstop1 <- Pars[4]
    a <- Pars[5]
    b <- Pars[6]
    c <- Pars[7]

    #
    with(as.list(c(State, Pars)), {
    #
    x <- Time
    ind1 <- as.numeric(x < Xs)
    #
    ITC.CFI <- ind1*(a+b*x+c*x^2)+(1-ind1)*((a+b*(Xs)+c*(Xs)^2)+(b+2*c*(Xs))*(x-(Xs)))  #returns the cubic function
    ITC.DFI <- ind1*(b+2*c*x) + (1 - ind1)*(b+2*c*(Xs))
    #
    onoff <- ifelse(Time>tbeg1 & Time<tstop1, 1, 0)
    CompFI <- (1-CumFI/ITC.CFI)*max.compFI1
    #
    dCumFI <- (onoff*(p1-1) + CompFI + 1)*ITC.DFI
    #
    return(list(c(dCumFI)))
    })
  }

  #Function to estimate parameters
  ODE.CFI.optim.2 <- function(Time, State, Pars) {
     #
     p1 <- Pars[1]  #percentage of decrease of ITD during the 1st perturbation
     max.compFI1 <- Pars[2]
     # tbeg1 <- Pars[3]
     # tstop1 <- Pars[4]
     # a           = Pars[5]
     # b           = Pars[6]
     # c           = Pars[7]
     # d           = Pars[8]

     #
     with(as.list(c(State, Pars)), {
     #
     x <- Time
     ind1 <- as.numeric(x < Xs)
     #
     ITC.CFI <- ind1*(a+b*x+c*x^2)+(1-ind1)*((a+b*(Xs)+c*(Xs)^2)+(b+2*c*(Xs))*(x-(Xs)))  #returns the cubic function
     ITC.DFI <- ind1*(b+2*c*x) + (1 - ind1)*(b+2*c*(Xs))
     #
     onoff <- ifelse(Time>tbeg1 & Time<tstop1, 1, 0)
     CompFI <- (1-CumFI/ITC.CFI)*max.compFI1
     #
     dCumFI <- (onoff*(p1-1) + CompFI + 1)*ITC.DFI
     #
     return(list(dCumFI))
     })
   }

  #  objective function
  ODE.CFI.obj.2 <- function(P, data){
     data <- Data.xy
     f <- ode(yinit, times, ODE.CFI.optim.2, P)
     optim.f <- sqrt(sum((data$CFI.plot- f[, 2])^2))
     return (optim.f)
   }

  #nls function for initial parameter
  ODE.CFI.obj.nls.2 <- function(p1,p2){
    P <- c(p1,p2)
    data <- Data.xy
    f <- ode(yinit, times, ODE.CFI.optim.2, P)
    return ( f[,2])
  }

  # ##############################################################################################################################
  # # ESTIMATE ALL PARAMETERS (INCLUDING THE TTC)
  # ##############################################################################################################################
  #
  ODE.CFI.obj.nls <- function(p1, p2, p3, p4){
   P = c(p1, p2, p3, p4)
   data = Data.xy
   f <- ode(yinit, times, ODE.CFI.optim, P)
   # optim.f = sqrt(sum((data$CFI.plot- f[,2])^2))
   return ( f[,2])
 }

  ODE.CFI.obj.nls1 <- function(p4){
   p2 = coef(st2)[2]
   p4 = coef(st2)[4]
   p3 = coef(st2)[3]
   P = c(p1, p2, p3, p4)
   data = Data.xy
   f <- ode(yinit, times, ODE.CFI.optim, P)
   # optim.f = sqrt(sum((data$CFI.plot- f[,2])^2))
   return ( f[,2])
 }

  ODE.CFI.optim.all <- function(Time, State, Pars) {
   #
   p1          = Pars[1]  #percentage of decrease of ITD during the 1st perturbation
   max.compFI1 = Pars[2]
   tbeg1       = Pars[3]
   tstop1      = Pars[4]
   a           = Pars[5]
   b           = Pars[6]
   c           = Pars[7]
   d           = Pars[8]

   #
   with(as.list(c(State, Pars)), {
     #
     ITC.CFI = a+b*Time+c*Time^2+d*Time^3
     ITC.DFI = b+2*c*Time+3*d*Time^2
     onoff = ifelse(Time>tbeg1 & Time<tstop1,1,0)
     compFI = (1-CumFI/ITC.CFI)*max.compFI1

     #
     dCumFI    =   (onoff*(1-p1) + (1 - onoff)) *
       ITC.DFI *(1+max.compFI1*(1 - CumFI/ITC.CFI))
     #
     return(list(c(dCumFI)))
   })
 }


   ##-----------------------------
   ##  objective function
   ##-----------------------------

  ODE.CFI.obj.all <- function(P, data){
   data = Data.xy
   f <- ode(yinit, times, ODE.CFI.optim.all, P)
   optim.f = sqrt(sum((data$CFI.plot- f[,2])^2))
   return (optim.f)
 }