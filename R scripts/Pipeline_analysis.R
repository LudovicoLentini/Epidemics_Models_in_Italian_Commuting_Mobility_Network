


################# PIPELINE for provinces
### needed: 
#dataframe of all fluxes
#df_fin
#Comune_Names
#read table of Comune names
library("readxl")
Comune_Names <- read_excel('Denominazione_Comuni.xlsx')
Comune_Names

Province <- read_excel('Dati Provinciali.xlsx')

#Province[Province$`Codice Provincia`==1,]$`Regione/Provincia`


############################## CHOOSE PROVINCE



prov_choose <- 28
prov_choose



province_name <- Province[Province$`Codice Provincia`==prov_choose,]$`Regione/Provincia`
province_name


df_part1 <- subset(df_fin,Resid_Prov %in% prov_choose & Dest_Prov %in% prov_choose)


#df_part <- df_fin[df_fin$Resid_Prov==prov_choose &
#                   df_fin$Dest_Prov== prov_choose,]

df_part <- df_part1


Comune_Names_prov<-Comune_Names[as.numeric(Comune_Names$`Codice Provincia`) == prov_choose,]
Comune_Names_prov


#handle partial dataframe
df_part$Resid <- as.character(df_part$Resid)
df_part$Dest <- as.character(df_part$Dest)



library(od)
mat <- od_to_odmatrix(df_part)
mat[is.na(mat)] <- 0

mat <- as.matrix(mat)

nrow(mat)
ncol(mat)


### if rows and columns of mat are different, 
if (nrow(mat) != ncol(mat)) {
  
  #find same comuni
  diff_comuni <- setdiff(rownames(mat), colnames(mat))
  same_comuni <- setdiff(rownames(mat), diff_comuni)
  
  
  
  print(paste('array of diff comuni',diff_comuni))
  
  #printing out Comunes left out for analysis
  for (el in diff_comuni){
    id_code <- as.numeric(el)
    print(Comune_Names[Comune_Names$Codice_Istat ==id_code ,]$Denominazione)
    print(paste('Population:',sum(mat[el,])))
  }
  
  #use only mat with same comuni
  mat <- mat[same_comuni,same_comuni]
}

nrow(mat)
ncol(mat)



## OD Matrix representation                                             
heatmap(mat, Colv = NA, Rowv = NA)



### GRAPH ANALYSIS  and representation

library(igraph)

#Graph for analysis

df_graph_values <- df_part
net_an <- graph_from_data_frame(df_graph_values,directed=TRUE)

cluster_tr <- round(transitivity(net_an),4)
transitivity(net_an)
#transitivity(net_an, type = "average")

#largest.cliques(net_an)

assortativity.degree(net_an, directed = TRUE)
reciprocity(net_an)

max_deg <- which.max(degree(net_an,mode=c('out')))
max_deg_idx <- names(max_deg)
max_deg_idx
Comune_Names[Comune_Names$Codice_Istat==max_deg_idx,]$Denominazione

which.min(degree(net_an,mode=c('out')))


edge.betweenness.community(net_an)

print(paste('Clustering coefficient is:',round(cluster_tr,3)))

### graph for representation
df_graph <- df_part
cutting_value = 500
df_graph[df_graph$Resid==df_graph$Dest, ]$Flux <- 0 
df_graph[df_graph$Flux<cutting_value,]$Flux <- 0
df_graph <- df_graph[df_graph$Flux>0,]


net <- graph_from_data_frame(df_graph,directed=TRUE)

E(net)$weight <- (log(df_graph$Flux)-min(log(df_graph$Flux)))*2

#edge_attr(net,as.character(Comune_Names$Codice_Istat) , index = E(net)) <- Comune_Names$Denominazione


#l <- layout_nicely(net)
l <- layout.graphopt(net)

plot(net,vertex.size=3, edge.arrow.size=0.05,vertex.label=NA, edge.width = E(net)$weight, layout =l, 
     main=paste('Province of',prov_choose,'. Clustering coeff:', cluster_tr)) 



### Simulation 


################  DEFINE WAY TO CHOOSE COMUNE ACCORDING TO SOMETHING

#store max and min population comune
id_min <- Comune_Names_prov[Comune_Names_prov$Popolazione == min(Comune_Names_prov$Popolazione),]$Codice_Istat
id_max <- Comune_Names_prov[Comune_Names_prov$Popolazione == max(Comune_Names_prov$Popolazione),]$Codice_Istat
id_min
id_max
########### CHOOSE ZONE OF FIRST INFECTED CASES
## the id of the location where the initial infected cases appeared

#zone0_id <- sample(rownames(mat),1)  #random sample


zone0_id <- max_deg_idx   # CHOOSE MAX OR MIN
name_id <- Comune_Names[Comune_Names$Codice_Istat==zone0_id,]$Denominazione

zone0_id
name_id


zone0_infected <- 10 # number of initial infected cases
beta <- 1.4 # the parameter controlling how often a susceptible-infected contact results in a new exposure
gamma <- 0.2 # the rate an infected recovers and moves into the resistant phase (1/Recovery_Time)
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
  #SEIR_list <- list(SEIR_nsim)
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
    
    #SEIR_list[[i+1]] <- SEIR_nsim
    
  }
  SEIR_tot <- data.frame(S_list,E_list,I_list,R_list)
  #SEIR_tot
  
  
  return(SEIR_tot)
}



df_sim <- Simulation_SEIR_pend(zone0_infected,zone0_id,
                               beta,gamma,sigma,max_sim,OD)

#df_sim


###PLOTTING
times <- seq(1,max_sim,1)

#name_id
#name_id

matplot(x = times, y = df_sim, type = "l",
        xlab = "Time", ylab = "SEIR", main = paste("SEIR Model,provincia:",province_name, 'initial:',name_id),
        lwd = 2, lty = 1, bty = "l", col = 1:4)

## Add legend
legend(length(times)-50, 0.7, c("Susceptible",'Exposed', "Infected", "Recovered"), pch = 1, col = 1:4, bty = "n")

which(df_sim$I_list==max(df_sim$I_list))





