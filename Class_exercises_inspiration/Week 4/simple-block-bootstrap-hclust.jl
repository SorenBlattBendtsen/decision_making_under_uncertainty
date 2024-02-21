#Import packages
using CSV
using DataFrames
using Plots
using Statistics
using Random
using Printf


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
plot(samples[:,1], label="Sample",  xlabel="Hours", ylabel="Demand MWh")

D = zeros(Float64, num_samples,num_samples)
for i = 1:num_samples
    for j = 1:num_samples
        D[i,j] = sqrt(sum((samples[l,i]-samples[l,j])*(samples[l,i]-samples[l,j])  for l = 1:sample_length))
    end
end

num_cluster = 5
result = hclust(D; linkage=:complete)
plot(reverse(result.heights), xlim=(1,20), legend=false, ylabel= "Distance of merge (based on linkage)", xlabel="# clusters")



a = cutree(result; k=4) # get the assignments of points to clusters
plot(samples, legend=false, xlim=(336,504), color=:lightgray)
