#Import packages
using CSV
using DataFrames
using Plots
using Statistics
using Random
using Distributions

#Set seed value for random number generator
Random.seed!(12345)

#Initialize normal distribution with mean and standard deviation
d = Normal(10.8208,1.4184)

#Draw 1000 random samples from the normal distribution
scenarios = rand(d, 1000)
histogram(scenarios, title="1000 samples", label="Samples")

#Alternative solution
#Sample uniform random numbers in [0,1]
norm_u = rand(Float64, 1000)
#Evaluate CDF by using the quantile function of the normal distribution
#The quantile function is the inverse CDF
scenarios_norm = [quantile(d, norm_u[x]) for x in 1:1000]
histogram(scenarios_norm, title="1000 samples", label="Samples")
