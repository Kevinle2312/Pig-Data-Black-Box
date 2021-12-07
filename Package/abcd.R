##################################################################              
#               ToPig data - Function file
#               Masoomeh and Hieu 04/12/2017
#               
#               File contains all functions needed
#
##################################################################
P <- NULL
pp <- NULL

#===============================================================
# A LINEAR FUNCTION for CFI
#===============================================================

#--------------------------------------------------------------
# Reparametrize a, b by X0 and ylast
#--------------------------------------------------------------
#function used for reparametrization in MAPLE 
# solve({0 = X0*b+a, Ylast = Xlast*b+a}, {a, b});
# a = (Ylast*X0)/(X0 - Xlast), b = -(Ylast)/(X0-Xlast)

abcd.0 <- function(P){
  X0 <- P[1]
  ylast <- P[2]
  a <- (ylast*X0)/(X0 - ylast)
  b <- -(ylast)/(X0-xlast)
  
  
  pp <- as.vector(c(a, b))
  return(pp)
}

#--------------------------------------------------------------
# NLS function
#--------------------------------------------------------------

nls.func.0 <- function(X0, ylast){
  pp <- c(X0, ylast)
  #calculation of a,b, c and d using these new parameters
  
  b <- abcd.0(pp)[2]
  a <- abcd.0(pp)[1]
  
  return (a+b*x)
}

#--------------------------------------------------------------
# Fit new parameters to a linear function of CFI
#--------------------------------------------------------------

pred.func.0 <- function(pr,age){
  #
  X0 <- pr[1]
  ylast <- P[2]
  
  #
  x <- age
  #calculation of a,b, c and d using these new parameters
  b <- abcd.0(pr)[2]
  a <- abcd.0(pr)[1]
  
  #
  results <- list()
  cfi <- a+b*x#returns the cubic function
  dfi <- b
  results[[1]] <- cfi
  results[[2]] <- dfi
  return (results)
}

#--------------------------------------------------------------------------------------------------
# Linear function of CFI curve and its 1st derivative (DFI) with original parameters (only a and b)
#--------------------------------------------------------------------------------------------------

pred.abcd.0 <- function(pr,age){
  #
  a <- pr[1]
  b <- pr[2]
  
  x <- age
  #calculation of a,b, c and d using these new parameters
  #
  results <- list()
  cfi <- a+b*x#returns the cubic function
  dfi <- b
  results[[1]] <- cfi
  results[[2]] <- dfi
  return (results)
}

#===============================================================
# A QUADRATIC FUNCTION for CFI
#===============================================================

#--------------------------------------------------------------
# Reparametrize a, b and c by X0, y2 and ylast
#--------------------------------------------------------------
#function used for reparametrization in MAPLE 
# solve({0=a+b*X0+c*X0**2,
# ylast = a+b*xlast+c*xlast**2},{a,b,c});

abcd.1 <- function(P){
  X0 <- P[1]
  y2 <- P[2]
  ylast <- P[3]
  
  a <- (xlast*ylast-4*y2*xlast+X0*ylast)*X0/(xlast^2-2*X0*xlast+X0^2)
  b <- -(xlast*ylast-4*y2*xlast+3*X0*ylast-4*y2*X0)/((xlast-X0)^2)
  c <- 2*(-2*y2+ylast)/(xlast^2-2*X0*xlast+X0^2)
  
  pp <- as.vector(c(a, b, c))
  return(pp)
}

#--------------------------------------------------------------
# NLS function
#--------------------------------------------------------------

nls.func.1 <- function(X0, y2, ylast){
  pp <- c(X0, y2, ylast)
  #calculation of a,b, c and d using these new parameters
  
  c <- abcd.1(pp)[3]
  b <- abcd.1(pp)[2]
  a <- abcd.1(pp)[1]
  
  return (a+b*x+c*x^2)
}

#--------------------------------------------------------------
# Fit new parameters to a quadratic function of CFI
#--------------------------------------------------------------

pred.func.1 <- function(pr,age){
  #
  X0 <- pr[1]
  y2 <- P[2]
  ylast <- P[3]
  #
  x <- age
  #calculation of a,b, c and d using these new parameters
  c <- abcd.1(pr)[3]
  b <- abcd.1(pr)[2]
  a <- abcd.1(pr)[1]
  
  #
  results <- list()
  cfi <- a+b*x+c*x^2  #returns the cubic function
  dfi <- b+2*c*x
  results[[1]] <- cfi
  results[[2]] <- dfi
  return (results)
}

#--------------------------------------------------------------------------------------------------
# Quadratic function of CFI curve and its 1st derivative (DFI) with original parameters (only a, b and c)
#--------------------------------------------------------------------------------------------------

pred.abcd.1 <- function(pr,age){
  #
  a <- pr[1]
  b <- pr[2]
  c <- pr[3]
  
  x <- age
  #calculation of a,b, c and d using these new parameters
  #
  results <- list()
  cfi <- a+b*x+c*x^2 #returns the cubic function
  dfi <- b+2*c*x
  results[[1]] <- cfi
  results[[2]] <- dfi
  return (results)
}

#===============================================================
# A QUADRATIC-LINEAR FUNCTION for CFI
#===============================================================

#--------------------------------------------------------------
# Reparametrize a, b and c by X0, Xs (breaking point),
# DFIs (DFI at Xs) and CFIs (CFI at Xs)
#--------------------------------------------------------------
#--------------------------------------------------------------
# Reparametrize a, b and c by X0, Xs (breaking point),
# DFIs (DFI at Xs) and CFIs (CFI at Xs)
#--------------------------------------------------------------
#function used for reparametrization in MAPLE 
# solve({
# 0=a+b*X_0+c*X_0**2,
# DFIs=b+2*c*Xs,CFIs=a+b*Xs+c*Xs**2},
# {a,b,c});

abcd.2 <- function(P){
  X0 <- P[1]
  Xs <- P[2]
  DFIs <- P[3]
  CFIs <- P[4]
  
  a <- -X0*(2*CFIs*Xs-CFIs*X0-Xs^2*DFIs+Xs*DFIs*X0)/(Xs^2-2*X0*Xs+X0^2)
  b <- (-Xs^2*DFIs+DFIs*X0^2+2*CFIs*Xs)/(Xs^2-2*X0*Xs+X0^2)
  c <- -(CFIs-Xs*DFIs+X0*DFIs)/(Xs^2-2*X0*Xs+X0^2)
  
  pp <- as.vector(c(a, b, c))
  return(pp)
}

#--------------------------------------------------------------
# NLS function
#--------------------------------------------------------------

nls.func.2 <- function(X0, Xs, DFIs, CFIs){
  pp <- c(X0, Xs, DFIs, CFIs)
  #calculation of a,b, c and d using these new parameters
  
  #-------------------------------------------------------
  # Function to reparameter a, b, c and of cubic function:
  #-------------------------------------------------------  
  c <- abcd.2(pp)[3]
  b <- abcd.2(pp)[2]
  a <- abcd.2(pp)[1]
  
  ind1 <- as.numeric(x < Xs)
  return (ind1*(a+b*x+c*x^2)+(1-ind1)*((a+b*(Xs)+c*(Xs)^2)+(b+2*c*(Xs))*(x-(Xs))))
}

#--------------------------------------------------------------
# Fit new parameters to a quadratic-linear function of CFI
#--------------------------------------------------------------

pred.func.2 <- function(pr,age){
  #
  X0 <- pr[1]
  Xs <- pr[2]
  DFIs <- pr[3]
  CFIs <- pr[4]
  #
  x <- age
  #calculation of a,b, c and d using these new parameters
  c <- abcd.2(pr)[3]
  b <- abcd.2(pr)[2]
  a <- abcd.2(pr)[1]
  
  #
  ind1 <- as.numeric(x < Xs)
  #
  results <- list()
  cfi <- ind1*(a+b*x+c*x^2)+(1-ind1)*((a+b*(Xs)+c*(Xs)^2)+(b+2*c*(Xs))*(x-(Xs)))  #returns the cubic function
  dfi <- ind1*(b+2*c*x) + (1 - ind1)*(b+2*c*(Xs))
  results[[1]] <- cfi
  results[[2]] <- dfi
  return (results)
}

#--------------------------------------------------------------------------------------------------
# Quadratic-linear function of CFI curve and its 1st derivative (DFI) with original parameters (only a, b and c)
#--------------------------------------------------------------------------------------------------

pred.abcd.2 <- function(pr,age){
  #
  a <- pr[1]
  b <- pr[2]
  c <- pr[3]
  
  x <- age
  #calculation of a,b, c and d using these new parameters
  #
  ind1 <- as.numeric(x < Xs)
  #
  results <- list()
  cfi <- ind1*(a+b*x+c*x^2)+(1-ind1)*((a+b*(Xs)+c*(Xs)^2)+(b+2*c*(Xs))*(x-(Xs)))  #returns the cubic function
  dfi <- ind1*(b+2*c*x) + (1 - ind1)*(b+2*c*(Xs))
  results[[1]] <- cfi
  results[[2]] <- dfi
  return (results)
}