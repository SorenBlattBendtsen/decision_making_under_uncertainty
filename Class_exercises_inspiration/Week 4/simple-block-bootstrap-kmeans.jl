#Import packages
using CSV
using DataFrames
using Plots
using Statistics
using Random
using Printf
using Clustering

Random.seed!(1234)

#Read demand data from demand.csv file
df = CSV.read("demand.csv", DataFrame)
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
clustsbb = plot(samples[:,1], label="Sample",  xlabel="Hours", ylabel="Demand MWh")

#Number of clusters
num_cluster = 5
#Apply k-means
result = kmeans(samples, num_cluster; display=:iter)
M = result.centers # get the cluster centers

a = assignments(result) # get the assignments of points to clusters

#Plot results
plot(samples, legend=false, xlim=(336,504), color=:lightgray)
plot!(M, legend=false, xlim=(336,504), color=:auto)
