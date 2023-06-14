
#test a small matrix 
vec <- c(2,1,2,8,7,8,0,2,1)
M1 <- matrix(t(vec), nrow=3,ncol=3)
#OD matrix
M1

#useful values
Dest_pop <-colSums(M1)
Orig_pop <- rowSums(M1)
Perc_t_mat <- round(t(M1)/colSums(M1),2)
Perc_t_mat
Perc_mat <- round(M1/rowSums(M1),2)
Perc_mat



#infectious rate
Beta_w <- 0.4
Beta_bus <- 0.7
#Beta_car <- 0.1
#Beta_car diversa??

#initialize contagious by Origin
Orig_inf <- c(0.1,0,0)


###notation: 4 phases
# Start of day, After commuting, after day of work,
# after second commuting

# _s means start,
# _d1 means day morn
# _d2 means day end
# _r means return 


##### STEP 1: TRAVEL 
#matrix multip for perc of infected traveling 
Mat_inf_s <- M1 * Orig_inf
Mat_inf_s

#update Matrix of infected in traveling 
Mat_inf_d1 <- Mat_inf_s + Mat_inf_s * Beta_bus
Mat_inf_d1

##### STEP 2 : WORK 
Dest_inf_d1 <- colSums(Mat_inf_d1)
Dest_inf_d1

#update list of infected at work
Incr_w <- Dest_inf_d1*Beta_w
Incr_w
Dest_inf_d2 <- Dest_inf_d1 + Incr_w
Dest_inf_d2

##### STEP 3 : TRAVEL BACK

#distribute the increment proportionally
Incr_r_t <-  Incr_w * Perc_t_mat
Incr_r_t

Mat_inf_d2_t <- t(Mat_inf_d1) + Incr_r_t 
Mat_inf_d2_t


## STEP 4: ARRIVE
Mat_inf_r_t <- Mat_inf_d2_t + Mat_inf_d2_t * Beta_bus
Mat_inf_r_t

#Number of Infected at end of day
Orig_inf_r <- colSums(Mat_inf_r_t)
Orig_inf_r

Inf_i <- Orig_inf_r/rowSums(M1)
Inf_i

Inf_i

Orig_inf_r *Perc_mat
M1
M1 * Inf_i


Update_day <- function(M1,Beta_w,Beta_bus,Orig_inf){

  ##### STEP 1: TRAVEL 
  #matrix multip for perc of infected traveling 
  Mat_inf_s <- M1 * Orig_inf
  #Mat_inf_s
  
  #update Matrix of infected in traveling 
  Mat_inf_d1 <- Mat_inf_s + Mat_inf_s * Beta_bus
  #Mat_inf_d1
  
  ##### STEP 2 : WORK 
  Dest_inf_d1 <- colSums(Mat_inf_d1)
  #Dest_inf_d1
  
  
  #update list of infected at work
  Incr_w <- Dest_inf_d1*Beta_w
  #Incr_w
  Dest_inf_d2 <- Dest_inf_d1 + Incr_w
  #Dest_inf_d2

  ##### STEP 3 : TRAVEL BACK
  Perc_t_mat <- round(t(M1)/colSums(M1),2)
  
  #distribute the increment proportionally
  Incr_r_t <-  Incr_w * Perc_t_mat
  #Incr_r_t
  
  Mat_inf_d2_t <- t(Mat_inf_d1) + Incr_r_t 
  #Mat_inf_d2_t
  
  ## STEP 4: ARRIVE
  Mat_inf_r_t <- Mat_inf_d2_t + Mat_inf_d2_t * Beta_bus
  #Mat_inf_r_t
  
  #Number of Infected at end of day
  Orig_inf_r <- colSums(Mat_inf_r_t)/rowSums(M1)
  
  return(Orig_inf_r)
}


#test function
vec <- c(2,1,2,8,7,8,0,2,1)
M1 <- matrix(t(vec), nrow=3,ncol=3)
M1 <- M1*1000
M1 
Dest_pop <-colSums(M1)
Orig_pop <- rowSums(M1)

#infectious rate
Beta_w <- 0.6
Beta_bus <- 0.8
#Beta_car <- 0.1

#initialize contagious by Origin
Orig_inf <- c(0.001,0,0)

Update_day(M1,Beta_w,Beta_bus,Orig_inf)


##test cycle

Orig_inf <- c(0.000,0.000,0.0001)
print(Orig_inf)

for (i in 1:10){
  inf_i <- Update_day(M1,0.3,0.6,Orig_inf)
  print(inf_i)
  Orig_inf <- inf_i
}





##test while
Orig_inf <- c(0.0001,0.000,0.000)
Update_day(M1,0.01,0.01,Orig_inf)
steps = 0
stop=FALSE
while (min(Orig_inf)<1 && stop==FALSE) {
  inf_i <- Update_day(M1,0.01,0.01,Orig_inf)
  print(inf_i)
  Orig_inf <- inf_i
  steps = steps + 1
  
  if (steps >100){
    stop=TRUE
#    Orig_inf <- 2
  }
}

print(steps)



####################### LOOP FOR DAYS COUNTING VARYING BETA WORK AND BETA BUS -----

sequence <- seq(0.01,1,0.01)
sequence
Res <- matrix(0,nrow = length(sequence),ncol = length(sequence))

i_idx <- 1
for(i in sequence){
  j_idx <- 1
  for (j in sequence){
    
    Orig_inf <- c(0.0001,0.000,0.000)
    steps = 0
    stop = FALSE
    while (min(Orig_inf)<1 && stop==FALSE ) {
      inf_i <- Update_day(M1,i,j,Orig_inf)
      #print(inf_i)
      Orig_inf <- inf_i
      steps = steps + 1
      
      if (steps >1000){stop = TRUE}
    }
    #print(steps)
    Res[i_idx,j_idx] <- steps
    
    j_idx <- j_idx + 1
  }
  i_idx <- i_idx +1
}


Res


