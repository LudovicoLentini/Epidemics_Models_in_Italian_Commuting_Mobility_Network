
#dataframe of all fluxes
#df_fin
#Comune_Names
#read table of Comune names
library("readxl")
Comune_Names <- read_excel('Denominazione_Comuni.xlsx')

Province <- read_excel('Dati Provinciali.xlsx')

unique(Comune_Names$`Codice Regione`)
Region_list <- c('Piemonte','Valdaosta','Lombardia','Trentino','Veneto',
                 'Friuli','Liguria','EmiliaRomagna','Toscana',
                 'Umbria','Marche','Lazio','Abruzzo','Molise',
                 'Campania','Puglia','Basilicata','Calabria',
                 'Sicilia','Sardegna')

############################## CHOOSE REgion -> NORTH ITALY REGIONS


#Region <- 'Lazio'
#Region_code <- which(Region_list==Region)

Region_code <- c(1,2,3,4,5,6)

prov_choose<- as.numeric(names(table(Comune_Names[as.numeric(Comune_Names$`Codice Regione`)==Region_code,]$`Codice Provincia`)))
prov_choose


df_part1 <- subset(df_fin,Resid_Prov %in% prov_choose & Dest_Prov %in% prov_choose)
df_part <- df_part1
remove(df_part1)


Comune_Names_prov<-subset(Comune_Names,as.numeric(`Codice Provincia`) %in% prov_choose)
Comune_Names_prov


#handle partial dataframe
df_part$Resid <- as.character(df_part$Resid)
df_part$Dest <- as.character(df_part$Dest)



##########CREATE MATRIX 
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
#heatmap(mat)



### GRAPH ANALYSIS  and representation

library(igraph)

#Graph for analysis

df_graph_values <- df_part
net_an <- graph_from_data_frame(df_graph_values,directed=TRUE)


glob_ef <- global_efficiency(net_an)
print(paste('Global efficiency is:',round(glob_ef,3)))
## SOme indicators (CLEAN AND DECIDE) 
cluster_tr <- round(transitivity(net_an),4)
print(paste('Clustering coefficient is:',round(cluster_tr,3)))

transitivity(net_an, type = "average")


assortativity.degree(net_an, directed = TRUE)
reciprocity(net_an)


#### Comunes with max and min degree
hist(degree(net_an),breaks = 100)

hist(degree(net_an, mode = c('total')),breaks = 100, main = 'Degrees of Nodes (Municipality): North Italy.',
     xlab = 'Degrees')


max(degree(net_an,mode=c('out')))
max_deg <- which.max(degree(net_an,mode=c('out')))
max_deg_idx <- names(max_deg)
max_deg_idx
Comune_Names[Comune_Names$Codice_Istat==max_deg_idx,]$Denominazione


min_deg <- which.min(degree(net_an,mode=c('out')))
min_deg_idx <- names(min_deg)
min_deg_idx
Comune_Names[Comune_Names$Codice_Istat==min_deg_idx,]$Denominazione



### graph for representation
df_graph <- df_part
cutting_value = 500
df_graph[df_graph$Resid==df_graph$Dest, ]$Flux <- 0 
df_graph[df_graph$Flux<cutting_value,]$Flux <- 0
df_graph <- df_graph[df_graph$Flux>0,]

df_graph$Resid
df_graph$Dest

length(union(df_graph$Resid,df_graph$Dest))


Vertex.df <- subset(Comune_Names,Codice_Istat %in% union(df_graph$Resid,df_graph$Dest))
length(Vertex.df$Codice_Istat)


net <- graph_from_data_frame(df_graph,directed=TRUE, vertices = Vertex.df)

E(net)$weight <- (log(df_graph$Flux)-min(log(df_graph$Flux)))*2
V(net)$size <- sqrt(Vertex.df$Popolazione)/80
V(net)$size


cutting_quant <- 0.98
cutting_Population <- quantile(Vertex.df$Popolazione,cutting_quant)
cutting_Population


V(net)$label <- ifelse(Vertex.df$Popolazione>cutting_Population,Vertex.df$Denominazione,'')
V(net)$label

V(net)$label.cex <- 1.2
#V(net)$label.cex

#l <- layout.graphopt(net)
l <- layout.davidson.harel(net)

par(mfrow = c(1, 1))
plot(net,vertex.size=V(net)$size, edge.arrow.size=0.05,vertex.label=V(net)$label,
     vertex.label.degree = -pi/2,
     edge.width = E(net)$weight, layout =l, 
     main=paste('North Italy. Global efficiency', round(glob_ef,3))) 





############# NEW -> INSERT VECTOR OF POPULATION 
#PRONE TO ERRORS: NEED CHECKING

length(rownames(mat))
Populat_vec <- c()
for ( i in 1:length(rownames(mat))){
  code <- as.numeric(rownames(mat)[i])
  pop <- Comune_Names[Comune_Names$Codice_Istat == code,]$Popolazione
  if (identical(pop, numeric(0))){
    print('Error?')
    print(code)
    pop <- 4209
  }
  Populat_vec <- c(Populat_vec,pop)
}

length(rownames(mat))==length(Populat_vec)



####### SIMULTION PIPELINE WITH DOUBLE PLOT
################  DEFINE WAY TO CHOOSE COMUNE ACCORDING TO SOMETHING

#store max and min population comune
id_min <- Comune_Names_prov[Comune_Names_prov$Popolazione == min(Comune_Names_prov$Popolazione),]$Codice_Istat
id_max <- Comune_Names_prov[Comune_Names_prov$Popolazione == max(Comune_Names_prov$Popolazione),]$Codice_Istat
########### CHOOSE ZONE OF FIRST INFECTED CASES
## the id of the location where the initial infected cases appeared
id_min
id_max
max_deg_idx
min_deg_idx

zone0_id <-min_deg_idx

zone0_id <- sample(rownames(mat),1)  #random sample
name_id <- Comune_Names[Comune_Names$Codice_Istat==zone0_id,]$Denominazione

zone0_id
name_id

zone0_infected <- 10 # number of initial infected cases
beta <- 0.5 # the parameter controlling how often a susceptible-infected contact results in a new exposure
gamma <- 0.1 # the rate an infected recovers and moves into the resistant phase (1/Recovery_Time)
sigma <- 0.2 # the rate at which an exposed person becomes infective (1/Incubation_Time)
max_sim <- 300 # maximum number of day simulated
OD <- mat # load the OD matrix
Pop <- Populat_vec #population vector


Simulation_SEIR_pend_tot <- function(zone0_infected,zone0_id,
                                     beta,gamma,sigma,max_sim,OD, Pop){
  
  
  
  #N_k <- rowSums(OD) # population for each cell
  N_k <- Pop
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
  
  
  Ret <- list(SEIR_tot,SEIR_list)
  return(Ret)
}



Ret <- Simulation_SEIR_pend_tot(zone0_infected,zone0_id,
                                beta,gamma,sigma,max_sim,OD,Pop)



## TOTAL RESULTS -> df_sim 
## PARTIAL (EVERY CITY) RESULTS IN df_list

df_sim <- Ret[[1]]
df_list <- Ret[[2]]



########CREATE DF FOR PLOTTING INF CITIES

Extract_SEIR <- function(SEIR_list, locat){
  RES <- matrix(NA,length(SEIR_list),4)
  colnames(RES) <- c('S','E','I','R')
  for (el in colnames(RES)){
    for(i in 1:length(SEIR_list)){
      RES[i,el] <- SEIR_list[[i]][locat,el][[1]]
    }
  } 
  return(data.frame(RES))
}


cutting_quant <- 0.99
cutting_Population <- quantile(Vertex.df$Popolazione,cutting_quant)
cutting_Population

Top_comun <- Comune_Names_prov[Comune_Names_prov$Popolazione>cutting_Population,]
num_top_comun <- length(Top_comun$Codice_Istat)
num_top_comun


##### Handle dataframe. extract only top COmuni

df_top <- list()
for (i in 1:num_top_comun){
  idx_extraction <- which(rownames(mat) == Top_comun$Codice_Istat[i])
  df_top[[i]] <-Extract_SEIR(df_list,idx_extraction) 
}

#Create Df for plotting
df_top_I <- matrix(NA,nrow=max_sim+1,ncol = num_top_comun)
for (i in 1:num_top_comun){
  df_top_I[,i] <- df_top[[i]]$I
}

df_top_I <- data.frame(df_top_I)


####DOUBLE PLOT WITH REAL NUMBER OF INFECTED
df_top_true_pop <- matrix(NA,nrow=max_sim+1, ncol = num_top_comun)
for (i in 1:num_top_comun){
  df_top_true_pop[,i] <- df_top_I[,i]*Top_comun$Popolazione[i]
}
df_top_true_pop <- data.frame(df_top_true_pop)/1000



## DOUBLE PLOTTING
par(mfrow = c(2, 1))

### PLOT 1
times <- seq(1,max_sim,1)

matplot(x = times, y = df_sim, type = "l",
        xlab = "Time", ylab = "SEIR",main = paste("SEIR Model: North Italy"),
        lwd = 2, lty = 1, bty = "l", col = 1:4)

text(max_sim/2,1.1,paste('Start at:',name_id))

## Add legend
legend('topright', c("Susceptible",'Exposed', "Infected", "Recovered"), pch = 1,cex=1, col = 1:4, bty = "n",y.intersp = 0.7)

#PLOT 2
times <- seq(1,max_sim+1,1)
matplot(x = times, y = df_top_true_pop, type = "l",ylim=c(0,max(df_top_true_pop)),
        xlab = "Time", ylab = "N. Infected (1000)", main = paste("Infection at different big cities"),
        lwd = 2, lty = 1, bty = "l", col = 1:num_top_comun)

## Add legend
legend(x=240,y=500,Top_comun$Denominazione,cex=1, col = 1:num_top_comun,pch=1, bty = "n",y.intersp = 0.7)

which(df_sim$I_list==max(df_sim$I_list))
