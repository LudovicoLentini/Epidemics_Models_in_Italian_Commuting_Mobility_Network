#test a small matrix 
vec <- c(2,1,1,8,7,8,0,2,1)
M1 <- matrix(t(vec), nrow=3,ncol=3)
M1 <- M1*10
#OD matrix
M1


#useful values
Orig_pop <- rowSums(M1)  # Pop of origin is row sums 
Dest_pop <- colSums(M1)  # Pop of destiantions is col sums

# Percentage by row and columns
Perc_mat <- round(M1/rowSums(M1),2)
Perc_mat      # Perc mat  is how the origin pop splits in the travel

Perc_t_mat <- round(t(M1)/colSums(M1),2)
Perc_t_mat  # perc Mat transposed is how the workers coming in are divided


#infectious rate ? 
Beta_w <- 0.4
Beta_bus <- 0.7
#Beta_car <- 0.1


## FIND A WAY TO START THE EPID.

Define_start_perc <- function(M,idx, N_infected) {
  Orig_pop <- rowSums(M)  # Pop of origin is row sums
  #empty perc
  Perc <- rep(0,nrow(M))
  #fill with value in idx pos ( if N_infected > Orig pop, then set Orig_pop)
  Perc[idx] <- min(N_infected,Orig_pop[idx])
  #Fraction
  Perc <- Perc/Orig_pop
  return(Perc)
}

M1

Orig_inf <- Define_start_perc(M1,idx=1, N_infected = 1)
Orig_inf

##### STEP 1: TRAVEL 
#matrix multip for perc of infected traveling 
Mat_inf_s <- M1 * Orig_inf
Mat_inf_s

#update Matrix of infected in traveling 
Mat_inf_d1 <- Mat_inf_s + Mat_inf_s * Beta_bus
Mat_inf_d1

M1 - Mat_inf_d1


##### STEP 2 : WORK 
Dest_inf_d1 <- colSums(Mat_inf_d1)
Dest_inf_d1



#update list of infected at work
Incr_w <- Dest_inf_d1*Beta_w
Incr_w
Dest_inf_d2 <- Dest_inf_d1 + Incr_w
Dest_inf_d2


#distribute the increment proportionally
Incr_r_t <-  Incr_w * Perc_t_mat
Incr_r_t

Mat_inf_d2_t <- t(Mat_inf_d1) + Incr_r_t 
Mat_inf_d2_t


## STEP 4: ARRIVE (Maybe can be skipped)
#Mat_inf_r_t <- Mat_inf_d2_t + Mat_inf_d2_t * Beta_bus
Mat_inf_r_t <- Mat_inf_d2_t
Mat_inf_r_t

#check if no one of this values is negative
check_M <- M1 - t(Mat_inf_r_t)
sum(check_M<0)

#Number of Infected at end of day
Orig_inf_r <- colSums(Mat_inf_r_t)
Orig_inf_r

Inf_i <- Orig_inf_r/rowSums(M1)
Inf_i



Update_day_v2 <- function(M1,Beta_w,Beta_bus,Orig_inf){
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
  
  
  #distribute the increment proportionally
  Incr_r_t <-  Incr_w * Perc_t_mat
  Incr_r_t
  
  Mat_inf_d2_t <- t(Mat_inf_d1) + Incr_r_t 
  Mat_inf_d2_t
  
  
  ## STEP 4: ARRIVE (Maybe can be skipped)
  #Mat_inf_r_t <- Mat_inf_d2_t + Mat_inf_d2_t * Beta_bus
  Mat_inf_r_t <- Mat_inf_d2_t
  Mat_inf_r_t
  
  #check if no one of this values is negative
  if (sum((M1 - t(Mat_inf_r_t))<0)>0){
    Mat_ret <- Orig_inf
  }
  
  #Number of Infected at end of day
  Orig_inf_r <- colSums(Mat_inf_r_t)
  #Orig_inf_r
  Inf_i <- Orig_inf_r/rowSums(M1)
  Mat_ret <- Inf_i
  return(Mat_ret)
}



Orig_inf <- Define_start_perc(M1,idx=1, N_infected = 1)
Orig_inf
Update_day_v2(M1,Beta_w,Beta_bus,Orig_inf)
