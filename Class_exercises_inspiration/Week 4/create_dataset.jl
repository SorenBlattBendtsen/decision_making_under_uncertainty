#Import packages
using CSV
using DataFrames
using Plots
using Statistics
using Random
using Distributions

Random.seed!(256)
# norm_dist = Normal(0,1)
#
# num_clusters = 4
# centers_x = [2.0, 5.0, 9.0, 1.5]
# centers_y = [2.5, 1.5, 9.0, 8.0]

# norm_dist = Normal(0,0.5)
#
# num_clusters = 4
# centers_x = [3.0, 7.0, 9.0, 1.5]
# centers_y = [2.5, 1.0, 8.0, 8.0]

norm_dist = Normal(0,0.5)

num_clusters = 4
centers_x = [3.0, 7.0, 9.0, 1.5]
centers_y = [3.5, 1.0, 8.0, 8.0]


datapoints = 15
dataset_x = zeros(Float16, datapoints*num_clusters)
dataset_y = zeros(Float16, datapoints*num_clusters)
deviation = rand(norm_dist, (num_clusters, datapoints,2))
idx = 1
for i = 1:num_clusters
    for j = 1:datapoints
        dataset_x[idx] = centers_x[i] + deviation[i,j,1]
        dataset_y[idx] = centers_y[i] + deviation[i,j,2]
        global idx = idx + 1
    end
end
scatter(dataset_x, dataset_y, series_annotations = text.(1:length(dataset_x), :bottom), border=:none,axis=nothing,legend= false,linetype = :scatter, markershape = :auto)

#initplot = scatter(dataset_x, dataset_y, series_annotations = text.(1:length(dataset_x), :bottom), border=:none,axis=nothing,legend= false,linetype = :scatter, markershape = :auto)
initplot = scatter(dataset_x, dataset_y,  border=:none,axis=nothing,legend= false,linetype = :scatter, markershape = :auto)

savefig(initplot, "plots\\Init.pdf")
