#Import packages
using CSV
using DataFrames
using Plots
using Statistics
using Random
using Printf


# set random seed
Random.seed!(1234)

#Read demand data from demand.csv file
df = CSV.read("../03-Scenario generation/Exercise/elprices2020.csv", DataFrame)
# transform data to array
data = Array(df)

# plot distribution
histogram(data, label="Observations", xlabel="Price DKK/MWh", ylabel="Frequency")
# plot time series
plot(data, label="Observations", xlabel="Hour", ylabel="Price DKK/MWh")

#Sample length
n = length(data)
#Block size
blocksize = 168

#Split sample in blocks of length blocksize
blocks = [data[i:min(i + blocksize - 1, end)] for i in 1:blocksize:n]

#Number of sampled_scenarios
num_sampled_scenarios = 250
sample_length = 168*4

# initialise sample array of zeros
sampled_scenarios = zeros(Float64,sample_length, num_sampled_scenarios)

# for each sample i
for i = 1:num_sampled_scenarios
    t = 1
    # for each time step t
    while t < sample_length
        # Draw random block
        r = rand(1:length(blocks))
        # Copy block content to sample
        # j runs over time steps within block
        for j = 0:blocksize-1
            sampled_scenarios[t+j,i] = blocks[r][j+1]
        end
        t = t+ blocksize
    end

end

#Plot sample 1
plot(sampled_scenarios, label="Sample",  xlabel="Hours", ylabel="Demand MWh", legend=false)

#######################

#Calculate distance matrix (euclidean distance)
D = zeros(Float64, num_sampled_scenarios,num_sampled_scenarios)
for i = 1:num_sampled_scenarios
    for j = 1:num_sampled_scenarios
        D[i,j] = sqrt(sum((sampled_scenarios[l,i]-sampled_scenarios[l,j])*(sampled_scenarios[l,i]-sampled_scenarios[l,j])  for l = 1:sample_length))
    end
end

#Number of scenarios to selected
num_reduced = 10
#Initialize equiprobable probabilities
probabilities = repeat([1.0/num_sampled_scenarios], 1, num_sampled_scenarios)[1,:]

#Include fast forward selection and apply it
include("fast-forward-selection.jl")
result = FastForwardSelection(D, probabilities, num_reduced)
#Resulting probabilities
print(result[1])
#Selected scenario indices
print(result[2])

#Plot initial sampled_scenarios
figure = plot(sampled_scenarios, legend=false, xlim=(336,504), color=:lightgray)

#Plot reduced scenarios
reduced_data = zeros(Float64, sample_length, num_reduced)
for i = 1:num_reduced
    reduced_data[:,i] = sampled_scenarios[:,result[2][i]]
end
plot!(reduced_data, legend=false, xlim=(336,504), color=:auto)
savefig(figure, "scenred-blockbootstrap-ffs-result.pdf")
