#Import packages
using CSV
using DataFrames
using Plots
using Statistics
using Random
using Distributions

#Set seed of random number generator
Random.seed!(1234)

#Set parameters
vehicles = 25
horizon = 1440
mean_arrival = 60
samples = 5

#Initialize empty scenario matrix
scenarios = zeros(Int,vehicles, samples)
#Initialize exponential distribution with mean
expdist = Exponential(mean_arrival)
#Generate scenarios
for s = 1:samples
    #Start with time 1 for each sample
    time = 0
    for v = 1:vehicles
        #For each vehicle draw a random number based on exponential distribution (Poisson process)
        arrival = rand(expdist)
        #Calculate time of arrival by advancing time
        time = time+round(arrival)
        #Save arrival time
        #@printf("u %0.4f -> t%i\n", arrival, time )
        scenarios[v,s] = min(horizon, time)
    end
end

#Plot result
plot(scenarios, xlabel="Vehicle", ylabel="Time", labels=["Scenario 1" "Scenario 2" "Scenario 3" "Scenario 4" "Scenario 5"], legend=:bottomright, markershape = :auto)
