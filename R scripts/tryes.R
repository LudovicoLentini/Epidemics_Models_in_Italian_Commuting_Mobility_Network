
Simulation_SEIR <- function(zone0_infected,zone0_id,
                            beta,gamma,sigma,max_sim,OD){
  
  N_k <- rowSums(M1) # population for each cell
  N_k_sum <- sum(N_k) # total population simulated
  
  locs_len <- length(N_k) # number of cells simulated
  
  
  # set up parameters for each cell
  beta_vec <- rep(beta, locs_len) # same transmission rate at each cell
  sigma_vec <- rep(sigma, locs_len) # same incubation-infectious transition rate at each cell
  gamma_vec <- rep(gamma, locs_len) # same recovery rate at each cell
  
  
  # set up the SEIR matrix
  SEIR <- matrix(nrow = locs_len, ncol = 4) # initiate an empty SEIR matrix
  colnames(SEIR) <- c("S", "E", "I", "R") # rename the vectors
  SEIR[, "S"] <- N_k # assign the number of susceptible people in each cell
  SEIR[, "E"] <- 0 # assign the number of exposed people in each cell
  SEIR[, "I"] <- 0 # assign the number of infected people in each cell
  SEIR[, "R"] <- 0 # assign the number of recovered people in each cell
  
  
  N_k_id <- seq(1,locs_len,1)
  
  # first infection
  first_infections <- (N_k_id == zone0_id) * zone0_infected
  first_infections
  # update the SEIR matrix
  SEIR[, "S"] <- SEIR[, "S"] - first_infections
  SEIR[, "I"] <- SEIR[, "I"] + first_infections
  
  # row normalize the SEIR matrix for keeping track of group proportions
  SEIR_n <- SEIR / rowSums(SEIR)
  SEIR_n[is.na(SEIR_n)] <- 0
  
  
  # make copy of the SEIR matrix
  SEIR_sim <- SEIR
  SEIR_nsim <- SEIR_n
  
  #store results
  SEIR_list <- list(SEIR_nsim)
  S_list <- c()
  E_list <- c()
  I_list <- c()
  R_list <- c()
  
  for(i in 1:max_sim){
    # New Exposed
    infected_mat <- replicate(locs_len, SEIR_nsim[, "I"])
    OD_infected <- round(OD * infected_mat) # people who are infected that travel to other locations
    
    inflow_infected <- colSums(OD_infected)
    total_inflow_infected <- sum(inflow_infected)
    #print(paste0("Total infected inflow: ", total_inflow_infected))
    
    
    new_exposed <-
      beta_vec * SEIR_sim[, "S"] * inflow_infected / (N_k + colSums(OD)) + # exposed by contacting with imported infected cases
      beta_vec * SEIR_sim[, "S"] * SEIR_sim[, "I"] / N_k # exposed by contacting with local infected cases
    new_exposed[is.na(new_exposed)] <- 0
    total_new_exposed <- round(sum(new_exposed))
    #print(paste0("New exposed: ", total_new_exposed))
    new_exposed <- ifelse(new_exposed > SEIR_sim[, "S"], SEIR_sim[, "S"], new_exposed) # make sure the N exposed is not bigger than N susceptible
    
    
    
    # New I
    new_infected <- sigma_vec * SEIR_sim[, "E"]
    total_new_infected <- round(sum(new_infected, na.rm = T))
    #print(paste0("New infected: ", total_new_infected))
    
    # New R
    new_recovered <- gamma_vec * SEIR_sim[, "I"]
    total_new_recovered <- round(sum(new_recovered, na.rm = T))
    #print(paste0("New recovered: ", total_new_recovered))
    SEIR_sim[, "S"] <- SEIR_sim[, "S"] - new_exposed
    SEIR_sim[, "E"] <- SEIR_sim[, "E"] + new_exposed - new_infected
    SEIR_sim[, "I"] <- SEIR_sim[, "I"] + new_infected - new_recovered
    SEIR_sim[, "R"] <- SEIR_sim[, "R"] + new_recovered
    SEIR_sim <- ifelse(SEIR_sim < 0, 0, SEIR_sim)
    
    # recompute the normalized SEIR matrix
    SEIR_nsim <- SEIR_sim / rowSums(SEIR_sim)
    SEIR_nsim[is.na(SEIR_nsim)] <- 0
    S <- sum(SEIR_sim[, "S"]) / N_k_sum
    E <- sum(SEIR_sim[, "E"]) / N_k_sum
    I <- sum(SEIR_sim[, "I"]) / N_k_sum
    R <- sum(SEIR_sim[, "R"]) / N_k_sum
    
    
    #store
    S_list <- c(S_list,S)
    E_list <- c(E_list,E)
    I_list <- c(I_list,I)
    R_list <- c(R_list,R)
    
    SEIR_list[[i+1]] <- SEIR_nsim
  }
  SEIR_tot <- data.frame(S_list,E_list,I_list,R_list)
  #SEIR_tot
  
  
  return(SEIR_tot)
}






#random
Dim <- 20
vec <- sample(100,Dim^2, replace=TRUE)

M1 <- matrix(t(vec), nrow=Dim,ncol=Dim)

M1

M1 <- M1*1000
#OD matrix
M1


zone0_infected <- 10 # number of initial infected cases
zone0_id <- sample(Dim,1) # the id of the location where the initial infected cases appeared (this one is in Shinjuku)
beta <- 1.5 # the parameter controlling how often a susceptible-infected contact results in a new exposure
gamma <- 0.1 # the rate an infected recovers and moves into the resistant phase (1/Recovery_Time)
sigma <- 0.2 # the rate at which an exposed person becomes infective (1/Incubation_Time)
max_sim <- 200 # maximum number of day simulated
OD <- M1 # load the OD matrix


df <- Simulation_SEIR(zone0_infected,zone0_id,
                beta,gamma,sigma,max_sim,OD)



times <- seq(1,length(df$S_list),1)
###PLOTTING
matplot(x = times, y = df, type = "l",
        xlab = "Time", ylab = "SEIR", main = "SEIR Model",
        lwd = 2, lty = 1, bty = "l", col = 1:4)

## Add legend
legend(length(times)-30, 0.7, c("Susceptible",'Exposed', "Infected", "Recovered"), pch = 1, col = 1:4, bty = "n")


