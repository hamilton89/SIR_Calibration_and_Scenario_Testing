# This model is an assignment for the Coursera Course "Interventions and Calibration" 
# available from: https://www.coursera.org/learn/interventions-and-calibration?

# In this assignment, students estimate parameter values from data on the prevalence of infection
# and then test the impact of a vaccination policy.

# packages
require(deSolve)
require(ggplot2)

# plot data
data <- read.csv("Graphics_and_Data/m3_nb3_data.csv")
plot(data$time, data$number_infected)

# initial state values
initial_state_values <- c(S = 499, I = 1, R = 0)

# times
times <- seq(from = 0, to = 200, by = 1) 

# SIR model function
SIR_fn <- function(time, state, parameters) {
    with(as.list(c(state, parameters)), {
        N  <- S+I+R
        dS <- -beta*S*I/N
        dI <- beta*S*I/N-gamma*I
        dR <- gamma*I
     
        return(list(c(dS, dI, dR)))
    })
}

# distance function
SIR_SSQ <- function(parameters, dat) {  
    # calculate model output using SIR function with ode()  
    beta <- parameters[1]
    gamma <- parameters[2]
    
    output <- as.data.frame(ode(y = initial_state_values
                              , times = times            
                              , func = SIR_fn            
                              , c(beta = beta, gamma = gamma))                  
    ) 
 
    # select elements where results$time is in data$time
    SSQ <- sum((output$I[output$time %in% dat$time]-dat$number_infected)^2)
    return(SSQ)
}

# optim function
optim(par = c(0.1, 0.1),        
      fn = SIR_SSQ,         
      dat = data) 

optimised

# Simulate the model with the estimated best-fitting parameter values
parameters <- c(beta = 0.16,
                gamma = 0.05)

output <- as.data.frame(ode(y = initial_state_values, 
                            times = times, 
                            func = SIR_fn,
                            parms = parameters))

# PLOT OF THE MODEL FIT

ggplot() +
  geom_line(data = output, aes(x = time, y = I)) +                              
  geom_point(data = data, aes(x = time, y = number_infected), col = "red") +  
  xlab("Time (days)")+                                              
  ylab("Prevalence of infection") +                                 
  labs(title = paste("Model fit to the epidemic curve with beta =", parameters["beta"], 
                     "and gamma =", parameters["gamma"]))

# confirm 92% hypothesis
initial_state_values <- c(S = 0.08*499,    # only 8% of the population are susceptible
                          I = 1,       
                          R = 0.92*499)    # 92% of the population are vaccinated/immune

output <- as.data.frame(ode(y = initial_state_values, 
                            times = times, 
                            func = SIR_fn,
                            parms = parameters))

# PLOT OF THE MODEL FIT

ggplot() +
  geom_line(data = output, aes(x = time, y = I)) +                              
  xlab("Time (days)")+                                              
  ylab("Prevalence of infection") + 
  labs(title = "Vaccine coverage of 92% with an all-or-nothing vaccine with 75% efficacy") +
  ylim(c(0,500))




