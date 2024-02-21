#Import packages
using CSV
using DataFrames
using Plots
using Statistics
using Random
using Printf
using Clustering

#Set random seed value
Random.seed!(1234)

#Read demand data from demand.csv file
df = CSV.read("Week 4/demand1.csv", DataFrame)
data = Array(df)
histogram(data, label="Observations", xlabel="Demand MWh", ylabel="Frequency")
plot(data, label="Observations", xlabel="Hour", ylabel="Demand MWh")

#Plot for a limited part of the x-axis
plot(data, label="Observations", xlabel="Hour", ylabel="Demand MWh", xlim=(336,504))

#Sample length
n = length(data)
#Block size
blocksize = 24

#Split sample in blocks of length blocksize
blocks = [data[i:min(i + blocksize - 1, end)] for i in 1:blocksize:length(data)]

#Number of samples
num_samples = 100
sample_length = n

samples = zeros(Float64,sample_length, num_samples)
for i = 1:num_samples
    t = 1
    while t < sample_length
        #Draw random block
        r = rand(1:length(blocks))
        #Copy block content to sample
        for j = 0:blocksize-1
            samples[t+j,i] = blocks[r][j+1]
        end
        t = t+ blocksize
    end

end

#Plot sample 1
plot(samples[:,1], label="Sample",  xlabel="Hours", ylabel="Demand MWh")

#Calculate distance matrix (euclidean distance)
D = zeros(Float64, num_samples,num_samples)
for i = 1:num_samples
    for j = 1:num_samples
        D[i,j] = sqrt(sum((samples[l,i]-samples[l,j])*(samples[l,i]-samples[l,j])  for l = 1:sample_length))
    end
end

#Number of scenarios to selected
num_reduced = 5
#Initialize equiprobable probabilities
probabilities = repeat([1.0/num_samples], 1, num_samples)[1,:]

#Include fast forward selection and apply it
include("fast-forward-selection.jl")
result = FastForwardSelection(D, probabilities, num_reduced)
#Resulting probabilities
print(result[1])
#Selected scenario indices
print(result[2])

#Plot initial samples
plot(samples, legend=false, xlim=(336,504), color=:lightgray)

#Plot reduced scenarios
reduced_data = zeros(Float64, sample_length, num_reduced)
for i = 1:num_reduced
    reduced_data[:,i] = samples[:,result[2][i]]
end
plot!(reduced_data, legend=false, xlim=(336,504), color=:auto)
