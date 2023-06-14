#read table of Comune names
library("readxl")
Comune_Names <- read_excel('Denominazione_Comuni.xlsx')
Comune_Names[Comune_Names$Codice_Istat==15010,]$Denominazione

library(data.table)
library(dplyr)
df <- fread('matrix_pendo2011_10112014.txt')

colnames(df) <- c('TipoRec','TipoRes','ProvR','ComR','Sex','Motivo',
                  'Luogo','ProvL','ComL','Stato','Mezzo','Orario',
                  'Tempo_imp','Stima_num','Num')

df[,"Resid"] <- df$ProvR*1000 + df$ComR
df[,"Dest"] <- df$ProvL*1000 + df$ComL


#cast chr to float
df$Stima_num <- as.numeric(df$Stima_num)


par(mfrow=c(3,2))

barplot(table(df$Mezzo),main='Mezzo')

barplot(table(df$Sex),main='Sesso')
barplot(table(df$Motivo),main='Motivo')

barplot(table(df$Tempo_imp))
barplot(table(df$Orario))
barplot(table(df$Orario))


#Dataset only work people
dfL <- df[df$TipoRec=='L',]

hist(dfL$Stima_num)
#cast chr to float
#dfL$Stima_num <- as.numeric(dfL$Stima_num)

'''
dfL11 <- dfL[dfL$Resid==1001]

a <- dfL11 %>% group_by(dfL11$Dest) %>% 
      summarize(Flussi = sum(Stima_num))
a
'''



### Final Dataframe with Origin, 
#Destination and flux intensityù
df_fin <- dfL %>% group_by(Resid,Dest) %>%
  summarize(Flux = sum(Stima_num),
  )

#LIst of all Comunalities

Comun_list <- unique(df_fin$Resid)
Comun_list2 <- unique(df_fin$Dest)


df_fin



###### GRAPH

library(igraph)
library(dplyr)

n_sample = 500
df_sampled = df_fin

cutting_value = 1000
df_sampled[df_sampled$Flux<cutting_value,]$Flux <- 0
df_sampled <- df_sampled[df_sampled$Flux>0,]

df_sampled <- data.frame(df_sampled)

table(df_sampled$Flux)
df_sampled <- df_sampled [df_sampled$Resid != df_sampled$Dest,]

net <- graph_from_data_frame(df_sampled,directed=TRUE)
#E(net)$weight <- log(df_sampled$Flux)/max(log(df_sampled$Flux)) + 2
E(net)$weight <- (log(df_sampled$Flux)-min(log(df_sampled$Flux)))*0.5

l <- layout_nicely(net)
plot(net,vertex.size=0.1, edge.arrow.size=0.01,vertex.label=NA, edge.width = 0.5, layout =l ) 



as.matrix(table(df_sampled$Resid,df_sampled$Dest))

mat



################generate matrix of one municipality (Milan,15)
df_parz1 <- df_fin[df_fin$Resid>15000,]
df_parz2 <- df_parz1[df_parz1$Resid<16000,]
df_parz3 <- df_parz2[df_parz2$Dest>15000,]
df_Milan <- df_parz3[df_parz3$Dest<16000,]

df_Milan

#install.packages('od')

library(od)


df_Milan$Resid <- as.character(df_Milan$Resid)
df_Milan$Dest <- as.character(df_Milan$Dest)

mat <- od_to_odmatrix(df_Milan)
mat[is.na(mat)] <- 0

mat <- as.matrix(mat)

mat
nrow(mat)
ncol(mat)

heatmap(mat, Colv = NA, Rowv = NA)


rowSums(mat)

min(rowSums(mat))

rowSums(mat)[rowSums(mat) == min(rowSums(mat))]
rowSums(mat)[rowSums(mat) == max(rowSums(mat))]
#colSums(mat)



## handling graph


df_Milan_graph <- df_Milan

cutting_value = 500
df_Milan_graph[df_Milan_graph$Resid==df_Milan_graph$Dest, ]$Flux <- 0 
df_Milan_graph[df_Milan_graph$Flux<cutting_value,]$Flux <- 0
df_Milan_graph <- df_Milan_graph[df_Milan_graph$Flux>0,]


net <- graph_from_data_frame(df_Milan_graph,directed=TRUE)

E(net)$weight <- (log(df_Milan_graph$Flux)-min(log(df_Milan_graph$Flux)))*2

#l <- layout_nicely(net)
l <- layout.graphopt(net)

plot(net,vertex.size=3, edge.arrow.size=0.05,vertex.label=NA, edge.width = E(net)$weight, layout =l) 




################ TRUE SIMULATION 


zone0_infected <- 10 # number of initial infected cases

## the id of the location where the initial infected cases appeared
zone0_id <- sample(rownames(mat),1)  #random sample
#zone0_id <- '15155' #Least populated zone
#zone0_id <- '15146' #Milan 

beta <- 0.6 # the parameter controlling how often a susceptible-infected contact results in a new exposure
gamma <- 0.1 # the rate an infected recovers and moves into the resistant phase (1/Recovery_Time)
sigma <- 0.2 # the rate at which an exposed person becomes infective (1/Incubation_Time)
max_sim <- 200 # maximum number of day simulated
OD <- mat # load the OD matrix


Simulation_SEIR_pend <- function(zone0_infected,zone0_id,
                            beta,gamma,sigma,max_sim,OD){
  
  N_k <- rowSums(OD) # population for each cell
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
  
  
  N_k_id <- rownames(OD)
  
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



df_sim <- Simulation_SEIR_pend(zone0_infected,zone0_id,
                      beta,gamma,sigma,max_sim,OD)

df_sim




###PLOTTING
times <- seq(1,max_sim,1)

matplot(x = times, y = df_sim, type = "l",
        xlab = "Time", ylab = "SEIR", main = paste("SEIR Model, initial: ",zone0_id),
        lwd = 2, lty = 1, bty = "l", col = 1:4)

## Add legend
legend(length(times)-40, 0.7, c("Susceptible",'Exposed', "Infected", "Recovered"), pch = 1, col = 1:4, bty = "n")

which(df_sim$I_list==max(df_sim$I_list))





