# dmultinomial_explorations.R
# examine features of the R stats version of the multinomial distribution

# start with simple case of flipping a coin
# what is probability of getting three heads and no tails?
pred <- c(0.5, 0.5) # fifty fifty shot of heads and tails
dmultinom(x = c(3,0), prob = pred) # 0.125 as expected
dmultinom(c(3,0), prob = pred) # can drop "x = " 
dmultinom(c(3,0), pred) # but cannot drop "pred = " due to size being expected after "x = "; perhaps this was the problem?

# this version of the multinomial is rounding non-integers
dmultinom(c(3,0), prob = pred)    # 0.125
dmultinom(c(3,0.49), prob = pred) # 0.125
dmultinom(c(3,0.5), prob = pred)  # 0.25
dmultinom(c(3,1), prob = pred)    # 0.25
dmultinom(c(3,1.49), prob = pred) # 0.25
dmultinom(c(3,1.5), prob = pred)  # 0.3125
dmultinom(c(3,2), prob = pred)    # 0.3125
dmultinom(c(3,2.49), prob = pred) # 0.3125

# what about RTMB's version of dmultinom? same
RTMB::dmultinom(c(3,0), prob = pred)    # 0.125
RTMB::dmultinom(c(3,0.49), prob = pred) # 0.125
RTMB::dmultinom(c(3,0.5), prob = pred)  # 0.25
RTMB::dmultinom(c(3,1), prob = pred)    # 0.25
RTMB::dmultinom(c(3,1.49), prob = pred) # 0.25
RTMB::dmultinom(c(3,1.5), prob = pred)  # 0.3125
RTMB::dmultinom(c(3,2), prob = pred)    # 0.3125
RTMB::dmultinom(c(3,2.49), prob = pred) # 0.3125

# so probably fine to just go with default dmultinom for now
# this shows why want to get away from multinomial distribution in stock assessments 
# (we don't end up with integers very often when multiply ESS by probs)


# can computed the log likelihood as ESS*sum(obs * log(pred)) (ignoring constants)
# where ESS is effective sample size, obs and pred are vectors of proportions that each sum to one
calc_multinom <- function(ESS, obs, pred){
  ESS * sum(obs * log(pred))
}

pred3 <- c(0.5, 0.3, 0.2)
my_prop1 <- seq(1,0.1, -0.01)
res <- rep(NA, length(my_prop1))
for (i in seq_along(my_prop1)){
  obs_props <- c(my_prop1[i], 1-my_prop1[i], 0) # changing observed proportions for bins 1 and 2, bin 3 is always zero
  res[i] <- calc_multinom(3, obs_props, pred3)
}
res # get different values for each value of my_prop1, so no rounding happening
plot(my_prop1, res)
