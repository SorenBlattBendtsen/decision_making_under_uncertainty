

include("create_dataset.jl")

using Clustering


dataset = transpose(hcat(dataset_x, dataset_y))
init_k = [4,13,53,51]

init_centers = zeros(2, num_clusters)
for j=1:num_clusters
    init_centers[:,j] = dataset[:,init_k[j]]
end

D = zeros(Float64, length(dataset_x),length(dataset_x))
for i = 1:length(dataset_x)
    for j = 1:length(dataset_x)
        D[i,j] = sqrt(((dataset_x[i]-dataset_x[j])*(dataset_x[i]-dataset_x[j])) + ((dataset_y[i]-dataset_y[j])*(dataset_y[i]-dataset_y[j])) )
    end
end

next = false
for i=0:10

    print(init_centers)
    result = kmedoids(D, num_clusters; maxiter=i, display=:iter, init=init_k)

    a = assignments(result) # get the assignments of points to clusters
    c = counts(result) # get the cluster sizes
    M = result.medoids # get the cluster centers

    clusters = hcat(a,hcat(dataset_x, dataset_y))
    colors = [:green, :purple, :blue, :orange]
    plotfig = scatter(clusters[:,2], clusters[:,3], group=clusters[:,1],color_palette=colors,border=:none,axis=nothing,markershape= :circle,markercolor=:match,legend= false)

    centroids = zeros(num_clusters,2)
    for j=1:num_clusters
        centroids[j,:] = dataset[:,M[j]]
    end
    centroids = hcat(range(1,length(M);step=1),centroids)
    scatter!(centroids[:,2],centroids[:,3], group=centroids[:,1],color_palette=colors,border=:none,axis=nothing,markershape= :star,markersize = 7,markercolor=:match,legend= false, title="Iteration $i")
    savefig(plotfig,"plots\\kmedoids-it$i.pdf")


    if next
        break
    end
    if result.converged
       global next = true
    end
end
