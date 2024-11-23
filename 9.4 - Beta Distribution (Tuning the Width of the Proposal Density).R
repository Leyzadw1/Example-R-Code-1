# CHAPTER   9 - Metropolis-Hastings Algorithm

# PROGRAM   9.4 - Beta Distribution (Tuning the Width of the Proposal Density)

# SUMMARY   
# In this program we use the Metropolis-Hastings algorithm to produce a sample 
# of observations from the Beta(2,4) distribution which has pdf
#
#   f_X(x)  = 20 x (1-x)^3    if  x \in [0,1]
#           = 0               if  x \notin [0,1]
#
# To achieve this we use the proposal density:
#       q(x,y) = 1/[ sqrt(2 \pi) \sigma ] exp{- 1/2 * (y-x)^2 / \sigma^2 }  
#       i.e. y ~ N( X_{t-1}, \sigma^2 )
#
# The effect of changing \sigma on the rate of convergence of the chain can
# be seen by modifying the 'proposal_sd' parameter in Section 2.1b.
#
#     Section 1:
#     Defining a function / macro that will facilitate the operations of the
#     Metropolis-Hastings algorithm.
#
#     Section 2:
#     Running the function / macro produced in the previous section and 
#     analyzing results.
#
#     Section 3:
#     Comparing the distribution of the post burn-in elements of the chain to an
#     exact iid sample drawn using the rbeta command




# Clear all variables before starting
rm(list = ls())




# SECTION 1 #
#------------------------------------------------------------------------------#
# Defining a function / macro that will facilitate the operations of the
# Metropolis-Hastings algorithm. Parameter definitions

    # n = Number of iterations / transitions
    # a = Parameter relating to the Beta(a,b) distribution
    # b = Parameter relating to the Beta(a,b) distribution
    # Start = Initial values of the chain (i.e. the value of X^(0) )
    # proposal_sd = Standard deviation of the proposal density 
    #               q = N( X_{t-1}, \sigma^2 )


# 1.1 - Defining the function / macro 
betaMH <- function(n, a, b, start, proposal_sd){
  
  # 1.1a - Create lists
  xbeta <- list()         # Define an empty list called xbeta
  x <- rep(NA,n)          # Create an empty list called 'x' of length n
  accepted <- rep(NA,n)   # Create an empty list called 'accepted' of length n
  
  
  # 1.1b - Initialize the chain
  x[1] <- start
  accepted[1] <- 1
  
  
  # 1.1c - Obtain the elements of the sequence {X^(t)}_{t=2}^n
  for (i in 2:n){
    
      # 1.1c1 - Draw a proposed state y from the proposal density 
      #         q(x,.) = N(x_{t-1}, proposal_sd)
      y <- rnorm(1, mean=x[i-1], sd=proposal_sd)
      
      # 1.1c2 - Define the acceptance ratio
      r <- dbeta(y,a,b) / dbeta(x[i-1],a,b)
      
        # dbeta(x,a,b) evaluates the beta(a,b) pdf at y, e.g.
        # dbeta(x,a,b)  = K x^{a-1} (1-x)^{b-1}   if  x in [0,1]
        #               = 0                       if  x notin [0,1]
        # where K = \Gamma(a+b) / [ \Gamma(a) \Gamma(b) ]
        # We could code this ourselves manually, but why bother.
      
      
      # 1.1c3 - Define the acceptance probability
      aprob <- min(1,r)
      
      # 1.1c4 - Determine whether proposed state y is accepted or rejected.
      # If uniform random variable < acceptance probability, then state y is accepted.
      if ( aprob >= runif(1) ){
        x[i] <- y
        accepted[i]<- 1
      } 
      # If uniform random variable >= acceptance probability, then state y is rejected.
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
burnin  <- 500
mcmc    <- 5000


# 2.1b - Running the function / macro defined in the previous section. 
# Results are saved in a list object called MH_Output
set.seed(12345)
MH_Output <- betaMH(burnin + mcmc, 2, 4, start = 0.5, proposal_sd = 0.5)


# 2.1c - Specify the acceptance rate.
acceptance_rate <- mean(MH_Output$accepted)
acceptance_rate


# 2.1d - Create a time series plot of all elements in the chain AFTER the burn-in
want_elements <- (burnin + 1):(burnin + mcmc)
plot(want_elements, MH_Output$x[want_elements], type = "l",
     main = "Time series plot", ylab = "x")


# 2.1e - Create autocorrelation plot of all elements in the chain after the burn-in
acf(MH_Output$x[want_elements], main = "Autocorrelation plot")


# 2.1f - Refining the elements of the chain produced by the Metropolis-Hastings
# algorithm to just the elements after the burn-in
post_burn <- tail(MH_Output$x, mcmc)




# SECTION 3 #
#------------------------------------------------------------------------------#
# Comparing the distribution of the post burn-in elements of the chain to an
# exact iid sample drawn using the rbeta command

# 4.1 - Producing two histograms showing:
#       a)  The distribution of post burn-in elements of the chain
#       b)  The distribution of values drawn from the rbeta command

# Create two histograms side by side
par(mfrow=c(1,2))

# Create histogram of the post burn-in elements of the chain
hist(post_burn, prob=T, main="MH sample", xlab="x")
lines(density(post_burn))

# Create histogram of the distribution of values drawn from the rbeta command
set.seed(123456)
xrbeta <- rbeta(mcmc, 2, 4)
hist(xrbeta, prob=T, main="Sample from rbeta()", xlab="x")
lines(density(xrbeta))











