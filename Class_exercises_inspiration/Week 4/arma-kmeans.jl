#Import packages
using CSV
using DataFrames
using Plots
using Statistics
using Random
using Printf
using Clustering
#Function to evaluate ARMA model
#yt1 = Y[t-1] (value at period t-1)
#et = epsilon[t] (error term at period t)
#et1 = epsilon[t-1] (error term at period t-1)
function ARMA(yt1,et,et1)
    return (0.75*yt1+et+0.5et1)
end

#Initialize seed for random number generator
Random.seed!(12345)

#Number of time steps in scenario
steps = 24
#Number of scenarios
scenarios = 100
#Last values of time series
y0=15
e0=0.0126

#Allocate data structure for samples
samples = zeros(Float64, steps, scenarios)
#Initialize normal distribution (here with standard deviation = 1)
norm = Normal(0.0,1.0)
#Sample all error terms for all scenarios and time steps
errors = rand(norm, (steps, scenarios))

#Add current initial values to data structur
samples = vcat(y0, samples)
errors = vcat(e0, errors)
print(sample)
print(errors)
println()

#Simulate ARIMA model using the defined function ARMA above
for s=1:scenarios
    for t=1:steps
        samples[t+1,s] = ARMA(samples[t,s],errors[t+1,s],errors[t,s])
        @printf("S%i t%i: yt-1=%0.4f et=%0.4f et-1=%0.4f --> yt=%0.4f\n" , s,t,samples[t,s], errors[t+1,s],errors[t,s], samples[t+1,s])
    end
end


#Add time index to plot data
plotdata = hcat(range(0,steps;step=1), samples)
#Plot scenarios
plot(plotdata[:,1],plotdata[:,2:scenarios+1],legend=false, color=:lightgray, xlabel="Time", ylabel="Y_t", label=["Scenario 1" "Scenario 2" "Scenario 3"])

#Number of clusters
num_cluster = 5
#Apply kmeans to sample data
result = kmeans(samples, num_cluster; display=:iter)
#Cluster centroids from kmeans
M = result.centers

#Add centroids to plot
plotdata2 = hcat(range(0,steps;step=1), M)
plot!(plotdata2[:,1],plotdata2[:,2:num_cluster+1], legend=false, color=:auto)
