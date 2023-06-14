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



df_fin$Resid_Prov <- floor(df_fin$Resid / 1000)
df_fin$Resid_Com <- df_fin$Resid - df_fin$Resid_Prov*1000


df_fin$Dest_Prov <- floor(df_fin$Dest / 1000)
df_fin$Dest_Com <- df_fin$Dest - df_fin$Dest_Prov*1000


df_fin








#printing function
Print_fun <- function(prev_c,next_c){
  if (prev_c!=next_c){
    prev_c <- next_c
    print(prev_c)
  }
  return(prev_c)
}


#Maker of Commuting Indicators
###############

#initialize empty lists
Autoff_list <- c()
Autdom_list <- c()
Autflux_list <- c()
Resident_w_list <- c()
Workers_list <- c()

n=length(Comun_list)
x <- 1
print(x)
for (id in Comun_list[1:n]) {
  #Print statements
  y <- floor(id/1000)
  x <- Print_fun(x,y)

  #Define Indicators
  Autflux <- max(0,df_fin[df_fin$Resid==id & df_fin$Dest==id,]$Flux)
  Resident_w <- sum(df_fin[df_fin$Resid==id,]$Flux)
  Workers <- sum(df_fin[df_fin$Dest==id,]$Flux)
  Autoff <- Autflux/Resident_w
  Autdom <- Autflux/Workers
  
  
  #Save indicators
  Autflux_list <- c(Autflux_list,Autflux)
  Resident_w_list <- c(Resident_w_list,Resident_w)
  Workers_list <- c(Workers_list, Workers)
  Autoff_list <- c(Autoff_list,Autoff)
  Autdom_list <- c(Autdom_list,Autdom)
}

#create DF
df_Indic <- data.frame(Comun_list[1:n],Resident_w_list,
                       Workers_list,Autflux_list,
                       Autoff_list,Autdom_list)

df_Indic


write.csv2(df_Indic,file="Pend_Indic.csv")




####
################# PLOTTING GRAPH

library(igraph)
library(dplyr)

n_sample = 5000
df_sampled = df_fin[1:n_sample,]

cutting_value = 70
df_sampled[df_sampled$Flux<cutting_value,]$Flux <- 0
df_sampled <- df_sampled[df_sampled$Flux>0,]


table(df_sampled$Flux)
df_sampled <- df_sampled [df_sampled$Resid != df_sampled$Dest,]

net <- graph_from_data_frame(df_sampled,directed=TRUE)
#E(net)$weight <- log(df_sampled$Flux)/max(log(df_sampled$Flux)) + 2
E(net)$weight <- (log(df_sampled$Flux)-min(log(df_sampled$Flux)))*1.5
hist(E(net)$weight)



l <- layout_nicely(net)
plot(net,vertex.size=0.1, edge.arrow.size=0.01,vertex.label=NA, edge.width = E(net)$weight, layout =l ) 


#plot(zach,vertex.size=10,vertex.label=NA)
#net <- graph_from_data_frame(d=links, vertices=nodes, directed=T)
