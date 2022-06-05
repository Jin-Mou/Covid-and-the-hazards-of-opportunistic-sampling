### Work Group 43 by Yixin Wang and Jin Mou
### Github repo address https://github.com/JinMou-JM/WG43P2.git and invitation has been sent


## SEIR is a function to simulate the COVID process. 
## S(Susceptible, 0) -> E(Exposed, 1) -> I(Infected, 2) 
## -> R(Recovery or serious disease, 3)
SEIR <- function(n = 5500000, ini_E = 10, days = 150, lambda = 0.4/n, prob_I = 1/3, prob_R = 1/5) {
  ## n = population size; lambda = overall viral infectivity parameter
  ## ini_E = initial state of exposed people
  ## days = the overall simulating days
  ## prob_I = daily probability of entering the infectious state
  ## prob_R = daily probability of leaving the infectious state
  ## beta = a relative contact rate with other people
  beta <- rlnorm(n,0,0.5); beta <- beta/mean(beta)
  beta_threshold <- sort(beta)[n * 0.1] ## the maximal beta among first 10% lowest beta
  
  sample_index <- sample(1:n, 0.001*n, replace = FALSE) ## randomly select 0.1% of the population
  
  x <- rep(0, n) ## initialize to susceptible state
  x[sample(1:n, 10, replace = FALSE)] <- 1 ## randomly choosing 10 people with E state
  S <- E <- I <- R <- rep(0, days) ## initialize the number of people with each state every day
  E[1] <- ini_E; S[1] <- n - ini_E ## initialize the number of people with states on day 1
  
  ## new_I_1 = number of new infections each day
  ## new_I_2 = number of new infections among the 10% of the population with the lowest beta each day
  ## new_I_3 = number of new infections in a random sample of 0.1% of the population
  new_I_1 <- new_I_2 <- new_I_3 <- rep(0, days) ## initialize the number matrix
  new_I_1[1] <- 0; new_I_2[1] <- 0; new_I_3[1] <- 0 ## assume there is not any infection at day 1
  
  ## loop over days
  for (d in 2:days) { 
    u <- runif(n) ## used to justify whether the state of people change
    x[x==2 & u<=prob_R] <- 3 ## I -> R with probability prob_R
    
    ## number of new infections in each sample
    new_I_1[d] <- length(x[x==1 & u<=prob_I]) 
    new_I_2[d] <- length(x[x==1 & u<=prob_I & beta<=beta_threshold]) 
    random_sample <- x[sample_index] ## the state of random_sample
    corresponding_beta <- beta[sample_index] ## the corresponding beta of random sample
    corresponding_u <- u[sample_index] ## the corresponding u of random sample
    new_I_3[d] <- length(random_sample[random_sample==1 & corresponding_u<=prob_I])
    
    x[x==1 & u<=prob_I] <- 2 ## E -> I with probability prob_I
    prob_E <- lambda * beta * sum(beta[x == 2]) ## prob_E is the probability of S -> E
    x[x==0 & u<=prob_E] <- 1 ## S -> E with probability prob_E
    
    ## recording the number of people with each state everyday
    S[d] <- sum(x == 0); E[d] <- sum(x == 1)
    I[d] <- sum(x == 2); R[d] <- sum(x == 3)
  }
  ## combine results to a list
  list(S=S, E=E, I=I, R=R, new_I_1=new_I_1, new_I_2=new_I_2, new_I_3=new_I_3, days=1:days)
}

par(mfcol=c(3,4))##set plot window up for multiple plots

drawing <- function(){  
  ##a function for plot
  list_data <- SEIR() ## run simulation
  plot(list_data$new_I_1/10,ylim=c(0,max(list_data$new_I_1)/10),xlab="day",ylab="N",col='black') 
  ## curve for the standardized whole population
  lines(list_data$days, list_data$new_I_1/10)##connect the dots
  ywholemaxbf <- list_data$new_I_1[which.max(list_data$new_I_1)]##find the peak value's y(the number of infectors)
  ywholemax <- ywholemaxbf/10##standardize the peak value's y
  xwholemaxbf <- list_data$days[which.max(list_data$new_I_1)]##find the peak value's day
  xwholemax <- xwholemaxbf ## get the peak value'day
  text(55,13000,c("total population"))##label for the curve
  text(xwholemax,ywholemax,c("totalmax"))##label for the peak
  
  points(list_data$days, list_data$new_I_2,col='red') 
  ## curve for the crowd who have downloaded ZOE
  lines(list_data$days, list_data$new_I_2, col='red')##connect the dots
  yzoemaxbf <- list_data$new_I_2[which.max(list_data$new_I_2)]##find the peak value's y(the number of infectors)
  yzoemax <-  yzoemaxbf##get the peak value's y
  xzoemaxbf <- list_data$days[which.max(list_data$new_I_2)]##find the peak value's day
  xzoemax <-  xzoemaxbf##get the peak value's day
  text(110,7000,c("ZOE"),col='red')##label for the curve
  text(xzoemax,yzoemax+1500,c("ZOEmax"),col='red')##label for the peak
  
  points(list_data$days,list_data$new_I_3*100,col='blue')
  ## curve for the sample
  lines(list_data$days,list_data$new_I_3*100,col='blue')##connect the dots
  ysamplemaxbf <- list_data$new_I_3[which.max(list_data$new_I_3)]##find the peak value's y(the number of infectors)
  ysamplemax <-  ysamplemaxbf*1000##standardize the peak value's y
  xsamplemaxbf <- list_data$days[which.max(list_data$new_I_3)]##find the peak value's day
  xsamplemax <- xsamplemaxbf##get the peak value's day
  text(60,3000,c("sample"),col='blue')##label for the curve
  text(xsamplemax,ysamplemax,c("samplemax"),col='blue')##label for the peak
}

replicate(10,drawing())

##It's easy to find that the peak value of infectors who have downloaded ZOE always appears later than the 
##whole population, at the same time, the peak value of the former is also significantly lower than the latter.
##This is probably because the former is more cautious about Covid and pay more attention to protect themselves,
##which leads to the result that in the preliminary stage, the situation among this crowd is much better compared
##with the whole population.
##However, as time goes by, due to the reason that we don't consider the situation that people get infected
##more than once, it's clear that in the later stage, the infectors are always the crowd who are cautious about 
##Covid although they might not download ZOE, and that's the reason why the daily infection trajectories of the 
##former and the latter are similar and even coincide during this period.



