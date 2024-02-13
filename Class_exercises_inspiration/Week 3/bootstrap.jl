#Import packages
using CSV
using DataFrames
using Plots
using Statistics
using Random
using Printf


Random.seed!(1234)

#Read demand data from demand.csv file and plot it
df = CSV.read("demand.csv", DataFrame)
data = Array(df)
histogram(data, label="Observations", xlabel="Demand MWh", ylabel="Frequency")
plot(data, label="Observations", xlabel="Hour", ylabel="Demand MWh")

#Sample length
n = length(data)

#Number of samples to draw
num_samples = 1000
sample_length = n

samples = zeros(Float32,sample_length, num_samples)
for i = 1:num_samples
    for j=1:sample_length
        #Draw random observation from original sample
        r = rand(1:n)
        samples[j,i] = data[r]
    end
end

#Plot of the first sample
histogram(samples[:,1], label="Sample 1",  xlabel="Demand MWh", ylabel="Frequency")
plot(samples[:,1], xlabel="Hour", ylabel="Demand MWh", label="Sample")
