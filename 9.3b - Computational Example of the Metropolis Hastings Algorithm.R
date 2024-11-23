# CHAPTER   9 - Metropolis-Hastings Algorithm

# PROGRAM   9.3b - Computational Example of the Metropolis Hastings Algorithm

# SUMMARY:
#
# Define the probability density function 
#       f_x(x) = 1/sqrt(2 \pi) exp{- 1/2 x^2}
#
# Define the real valued, measurable function
#       g(x)   = exp^{ -1/2 (x - 3)^2} + exp^{ -1/2 (x - 6)^2}
#
# In this program we approximate the integral 
#       E[ g(x) ] = \int_R g(x) f_X(x) dx
# by applying Monte Carlo integration to the results of the Metropolis-Hastings
# algorithm. 
#
# Herein, we use the proposal density:
#       q(x,y) = 1/sqrt(2 \pi) exp{- 1/2 (y-x)^2}  i.e. N( X_{t-1}, 1 )
# to generate proposed new states for the chain.
#
# The program is structured as follows:
#
#     Section 1
#     Defining a function / macro that will facilitate the operations of the
#     Metropolis-Hastings algorithm.
#
#     Section 2
#     Running the function / macro produced in the previous section and 
#     analyzing results.
#
#     Section 3
#     Performing Monte Carlo integration using the output of the Metropolis
#     Hastings algorithm
#
#     Section 4
#     Comparing the distribution of the post burn-in elements of the chain to a
#     normal distribution produced in R



# Clear all variables before starting
rm(list = ls())


# SECTION 1 #
#------------------------------------------------------------------------------#
# Defining a function / macro that will facilitate the operations of the
# Metropolis-Hastings algorithm. Parameter definitions:
    # n     = Number of iterations of the chain
    # Start = Initial values of the chain (i.e. the value of X^(0) )

set.seed(12345)

# 1.1 - Defining the function / macro 
MH_Algoritm <- function(n, start){
  
  # 1.1a - Create lists
  chain <- list()         # Define an empty list called chain
  x <- rep(NA,n)          # Create an empty list called 'x' of length n
  accepted <- rep(NA,n)   # Create an empty list called 'accepted' of length n
  
  
  # 1.1b - Initialize the chain
  x[1] <- start
  accepted[1] <- 1
  
  
  # 1.1c - Obtain the elements of the sequence {X^(t)}_{t=2}^n
  for (i in 2:n){
    
    # 1.1c1 - Draw a proposed state y from the proposal density 
    #         q(x,.) = N(x_{t-1}, 1)
    y <- rnorm(1, mean=x[i-1], sd=1)
    
    # 1.1c2 - Define the acceptance ratio
    r <- exp( -1/2*(y^2 - x[i-1]^2) )
    
    # 1.1c3 - Define the acceptance probability
    aprob <- min(1,r)
    
    # 1.1c4 - Determine whether proposed state y is accepted or rejected.
    # If acceptance probability >= uniform random variable, then state y is accepted.
    if ( aprob >= runif(1) ){
      x[i] <- y
      accepted[i]<- 1
    } 
    # If acceptance probability <  uniform random variable, then state y is rejected.
    else {
      x[i] <- x[i-1]
      accepted[i] <- 0
    }
  }
  
  # 1.1d - Return list of states and accept-reject occurrences
  return(list(x=x, accepted=accepted))
}

      # There are two outputs of this function / macro:
      # $x        -  A list of all the elements produced in the chain
      # $accepted -  A binary list of whether each proposed state was accepted or 
      #              rejected in each iteration



# SECTION 2 #
#------------------------------------------------------------------------------#
# Running the function / macro produced in the previous section and analyzing results

# 2.1a - Define burn-in and number of elements of the chain to use after the
# burn-in has occurred.
burnin  <- 1000
mcmc    <- 9000


# 2.1b - Running the function / macro defined in the previous section. 
# Results are saved in a list object called MH_Output
MH_Output <- MH_Algoritm(burnin + mcmc, start = 20)



# 2.1c - Specify the acceptance rate.
acceptance_rate <- mean(MH_Output$accepted)
acceptance_rate


# 2.1d - Create a time series plot of all elements in the chain
all_elements <- (1):(burnin + mcmc)
plot(all_elements, MH_Output$x[all_elements], type = "l",
     main = "Time series plot", ylab = "x")


# 2.1e - Create a time series plot of all elements in the chain AFTER the burn-in
want_elements <- (burnin + 1):(burnin + mcmc)
plot(want_elements, MH_Output$x[want_elements], type = "l",
     main = "Time series plot", ylab = "x")


# 2.1f - Create autocorrelation plot of all elements in the chain after the burn-in
acf(MH_Output$x[want_elements], main = "Autocorrelation plot")






# SECTION 3 #
#------------------------------------------------------------------------------#
# Performing Monte Carlo integration using the output of the Metropolis-Hastings
# algorithm

# 3.1a - Defining the function g(x) and create a graph of the function
g <- function(x){ exp(-(x-3)^2/2)+exp(-(x-6)^2/2) }
curve(g,from=0, to=10, xlab="Function",ylab="",lwd=2)


# 3.1b - Refining the elements of the chain produced by the Metropolis-Hastings
# algorithm to just the elements after the burn-in
post_burn <- tail(MH_Output$x, mcmc)

# 3.1c - Calculate \bar{g}(x) as an approximation to E[g(x)]
barg = sum(g(post_burn)) / mcmc
barg

  # [1] 0.07771766
  # This is very close to the true value of the integral = 0.0746





# SECTION 4 #
#------------------------------------------------------------------------------#
# Comparing the distribution of the post burn-in elements of the chain to a
# normal distribution produced in R

# 4.1 - Producing two histograms showing:
#       a)  The distribution of post burn-in elements of the chain
#       b)  The distribution of values drawn from the rnorm command

  # Create two histograms side by side
  par(mfrow=c(1,2))
  
  # Create histogram of the post burn-in elements of the chain
  hist(post_burn, prob=T, main="MH sample", xlab="x")
  lines(density(post_burn))
  
  # Create histogram of the distribution of values drawn from the rnorm command
  set.seed(123456)
  xrnorm<-rnorm(mcmc,0,1)
  hist(xrnorm, prob=T, main="Sample from rnorm()", xlab="x")
  lines(density(xrnorm))

  
  
# 4.2 - Conducting a two sided Kolmogorov-Smirnov test to test whether the data
# from either samples is drawn from the same distribution.
ks.test(post_burn, xrnorm)
  
    # p-value = 0.747  =>  No evidence against the null hypothesis that the
    # data from either sample are drawn from the same distribution
