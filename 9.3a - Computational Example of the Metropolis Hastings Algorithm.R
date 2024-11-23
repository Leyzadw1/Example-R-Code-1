# CHAPTER   9 - Metropolis-Hastings Algorithm

# PROGRAM   9.3a - Computational Example of the Metropolis Hastings Algorithm

# SUMMARY:
#
# Define the probability density function 
#       f_x(x) = A x e^(-x)   if x \in [0, 100] = \Omega
#                0            otherwise
# Here A is a normalizing constant.
#
# In this program we approximate the integral 
#       E[ X ] = \int_\Omega x f_X(x) dx
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




# Clear all variables before starting
rm(list = ls())



# SECTION 1 #
#------------------------------------------------------------------------------#
# Defining a function / macro that will facilitate the operations of the
# Metropolis-Hastings algorithm. Parameter definitions

  # n = number of elements of the chain
  # start = initial state of the chain
  # proposeal_sd = sigma

# 1.1 - Defining the function / macro 
MetHat <- function(n, start, proposal_sd){
  
    # 1.1a - Create blank vectors
    x <- rep(NA, n)
    accepted <- rep(NA, n)
    
    
    # 1.1b - Initialize the chain
    x[1] <- start
    accepted[1] <- 1
    
    
    # 1.1c - Obtain the elements of the sequence {X^(t)}_{t=2}^n
    for (i in 2:n){
      
      
        # 1.1c1 - Draw a proposed state y from the proposal density 
        #         q(x,.) = N(x_{t-1}, proposal_sd)
        y <- rnorm(1, mean=x[i-1], sd=proposal_sd)
        
        # 1.1c2 - Define f(y)
        if (y >= 0 && y <= 100) fy <- y * exp(-y) else fy <- 0
        
        # 1.1c3 - Define f(x)
        if (x[i-1] >= 0 && x[i-1] <= 100) fx <- x[i-1] * exp(- x[i-1]) else fx <- 0
        
        # 1.1c3 - Define acceptance ratio
        r <- fy / fx
        
        # 1.1c4 - Define the acceptance probability
        aprob <- min(1,r)
        
        # 1.1c5 - Define uniform random variable
        u <- runif(1)
        
        # 1.1c6 - Determine if proposed state is accepted or rejected
        # If acceptance probability >= uniform random variable, accept proposed state
        # If acceptance probability <  uniform random variable, reject proposed state
        if (aprob >= u) x[i] <- y
        if (aprob <  u) x[i] <- x[i-1]
        
        # 1.1c7 - Append result to accepted vector
        # If acceptance probability >= uniform random variable, accept proposed state
        # If acceptance probability <  uniform random variable, reject proposed state
        if (aprob >= u) accepted[i] <- 1
        if (aprob <  u) accepted[i] <- 0
      
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
set.seed(12345)
MH_Output <- MetHat(burnin + mcmc, start = 25, proposal_sd = 5)



# 2.1c - Specify the acceptance rate.
acceptance_rate <- mean(MH_Output$accepted)
acceptance_rate



# 2.1d - Create a time series plot of ALL elements in the chain
all_elements <- (0 + 1):(burnin + mcmc)
plot(all_elements, MH_Output$x[all_elements], type = "l",
     main = "Time series plot", ylab = "x")


# 2.1e - Create a time series plot of all0 elements in the chain AFTER burn-in
want_elements <- (burnin + 1):(burnin + mcmc)
plot(want_elements, MH_Output$x[want_elements], type = "l",
     main = "Time series plot", ylab = "x")



# 2.1f - Create autocorrelation plot of all elements in the chain AFTER burn-in
acf(MH_Output$x[want_elements], main = "Autocorrelation plot")



# 2.1g - Refining the elements of the chain produced by the Metropolis-Hastings
# algorithm to just the elements after the burn-in
post_burn <- tail(MH_Output$x, mcmc)



# 2.1h - Creating a histogram of the values of the chain produced by the 
# Metropolis-Hastings algorithm in relation to the target probability density
# function

  # Graph target pdf
  f <- function(x) {  (1 - 101*exp(-100) ) * x * exp(-x)   }
  xvals <- seq(0.01, 25, by = 0.01)
  fvals <- f(xvals)

  # Construct histogram
  hist(post_burn, nclass = 16, col = "grey",
       freq = FALSE, xlim = c(0, 10),  ylim = c(0, 0.4),
       xlab = "", 
       main = "Histogram of Observations from MH Algorithm")
  lines(xvals, fvals)


  
  
  

# SECTION 3 #
#------------------------------------------------------------------------------#
# Performing simple Monte Carlo integration





# 3.1b - Calculate approximate value of E[X]
g1 <- function(x){ x }
barg1 = sum(g1(post_burn)) / mcmc

  Expx <- barg1
  Expx
  
  # E[X] \approx 1.997814
  # This is very close to the true value of 2.0000


  
