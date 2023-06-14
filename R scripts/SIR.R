
## Load deSolve package
library(deSolve)



## Create an SIR function
sir <- function(time, state, parameters) {
  
  with(as.list(c(state, parameters)), {
    
    dS <- -beta * S * I
    dI <-  beta * S * I - gamma * I
    dR <-                 gamma * I
    
    return(list(c(dS, dI, dR)))
  })
}
M1
sir(0,init,parameters)

### Set parameters
## Proportion in each compartment: Susceptible 0.999999, Infected 0.000001, Recovered 0
#init       <- c(S = 1-1e-6, I = 1e-6, R = 0.0)
init       <- c(S = 1-1e-2, I = 1e-2, R = 0.0)
## beta: infection parameter; gamma: recovery parameter
parameters <- c(beta = 1.4247, gamma = 0.14286)
#parameters <- c(beta=0.8,gamma=0.5)
## Time frame
times      <- seq(0, 70, by = 1)

## Solve using ode (General Solver for Ordinary Differential Equations)
out <- ode(y = init, times = times, func = sir, parms = parameters)

## change to data frame
out <- as.data.frame(out)
## Delete time variable
out$time <- NULL
## Show data
head(out, 10)



## Plot
matplot(x = times, y = out, type = "l",
        xlab = "Time", ylab = "Susceptible and Recovered", main = "SIR Model",
        lwd = 1, lty = 1, bty = "l", col = 2:4)

## Add legend
legend(40, 0.7, c("Susceptible", "Infected", "Recovered"), pch = 1, col = 2:4, bty = "n")



sir(0,init,parameters)
parameters[1]*init[2]*init[1]


M1





1*0*parameters[1]



###function by hand 

#personal SIR

personal_sir <- function(init, parameters) {
  #definition of elements
  beta = parameters[1]
  gamma = parameters[2]
  S = init[1]
  I = init[2]
  R = init[3]
  
  
  #SIR function
  dS <- -beta * S * I
  dI <-  beta * S * I - gamma * I
  dR <-                 gamma * I
  
  #return
  res <- c(dS, dI, dR)
  names(res) <- c('dS','dI','dR')
  return(res)
}

#testing one output:
d_sir <- personal_sir(init, parameters)
d_sir

sir(1,init,parameters)



personal_SIR_model <- function(y = init, times = times, func = personal_sir, parms = parameters) {
  
  res <- matrix(NA,nrow = length(times), ncol = 3)
  
  j=1
  #loop for time
  for (el in times){
    y <- y + d_sir
    #append result
    res[j,] <- y
    
    
    d_sir <- func(1,y, parms)
    #y <- y + d_sir
    
    j=j+1
  }
  return(res)
}

out2 <- personal_SIR_model(y = init, times = times, func =sir, parms = parameters)

out2
out2 <- as.data.frame(out2)


## Plot
matplot(x = times, y = out2, type = "l",
        xlab = "Time", ylab = "Susceptible and Recovered", main = "SIR Model",
        lwd = 1, lty = 1, bty = "l", col = 2:4)

## Add legend
legend(40, 0.7, c("Susceptible", "Infected", "Recovered"), pch = 1, col = 2:4, bty = "n")


(out$S - out2$V1)[7]

plot(out$S - out2$V1)

out$S[7]
out2$V1[7]












##lets see if function is iterable 

init       <- c(S = 1-1e-6, I = 1e-6, R = 0.0)
## beta: infection parameter; gamma: recovery parameter
parameters <- c(beta = 1.4247, gamma = 0.14286)

## Time frame
times      <- seq(0, 1, by = 1)
list_SIR <- init
for (i in 1:10){
  out <- ode(y = list_SIR, times = times, func = sir, parms = parameters)
  list_SIR <- c(out[2,c(2:4)])
  print(list_SIR)
}


