#Import packages
using CSV
using DataFrames
using Plots
using Statistics
using Random
using Printf
using Distributions
using Clustering

# Set seed value for random number generator
Random.seed!(12345)

# Initialize triangular distributions
dist_wheat = TriangularDist(2.0,3.0,2.4)
dist_corn = TriangularDist(2.4,3.6,2.8)
dist_beets = TriangularDist(16.0,24.0,19.0)

# number of scenarios, we'd like to sample
scen_num = 100

# Sample uniform random numbers in [0,1]
norm_u = rand(Float64, scen_num)

# Evaluate CDF by using the quantile function of the normal distribution
# The quantile function is the inverse CDF
scenarios_wheat = [quantile(dist_wheat, norm_u[x]) for x in 1:scen_num]
scenarios_corn = [quantile(dist_corn, norm_u[x]) for x in 1:scen_num]
scenarios_beets = [quantile(dist_beets, norm_u[x]) for x in 1:scen_num]

# plot marginal sample distributions
histogram(scenarios_wheat, title="100 sampled_scenarios", label="Wheat")
histogram!(scenarios_corn,  label="Corn")
histogram!(scenarios_beets,  label="Sugar Beets")

# check correlations
scatter(scenarios_wheat, scenarios_corn, xlab = "Wheat", ylab = "Corn", legend = false)
scatter(scenarios_wheat, scenarios_beets, xlab = "Wheat", ylab = "Sugar Beets", legend = false)
scatter(scenarios_corn, scenarios_beets, xlab = "Corn", ylab = "Sugar Beets", legend = false)

############################

sampled_scenarios = zeros(Float64 , 3, scen_num)
sampled_scenarios[1,:] = scenarios_wheat
sampled_scenarios[2,:] = scenarios_corn
sampled_scenarios[3,:] = scenarios_beets


#Number of clusters
num_cluster = 5
#Apply kmeans to sample data
result = kmeans(sampled_scenarios, num_cluster; display=:iter)
#Cluster centroids from kmeans
M = result.centers

histogram(M[1,:], title="100 sampled_scenarios", label="Wheat")
histogram!(M[2,:],  label="Corn")
histogram!(M[3,:],  label="Sugar Beets")

#Assignments of data points to clusters
assigned = assignments(result) # get the assignments of points to clusters
new_probabilities = zeros(Float64, num_cluster)
for a in assigned
    new_probabilities[a] = new_probabilities[a] + (1/length(assigned))
end

println("Scenarios:")
for i=1:length(M[:,1])
    for j=1:length(M[i,:])
        @printf("%0.2f & ", M[i,j])
    end
    println("\\\\")
end

print(new_probabilities)
