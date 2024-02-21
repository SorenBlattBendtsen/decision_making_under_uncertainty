

include("create_dataset.jl")

using Clustering


dataset = transpose(hcat(dataset_x, dataset_y))
init_k = [4,12,53,51]

init_centers = zeros(2, num_clusters)
for j=1:num_clusters
    init_centers[:,j] = dataset[:,init_k[j]]
end

for i=0:10
    print(init_centers)
    result = kmeans(dataset, num_clusters; maxiter=i, display=:iter, init=init_k)

    a = assignments(result) # get the assignments of points to clusters
    c = counts(result) # get the cluster sizes
    M = result.centers # get the cluster centers
    print(M)

    clusters = hcat(a,hcat(dataset_x, dataset_y))
    colors = [:green, :purple, :blue, :orange]
    plotfig = scatter(clusters[:,2], clusters[:,3], group=clusters[:,1],color_palette=colors,border=:none,axis=nothing,markershape= :circle,markercolor=:match,legend= false)

    centroids = hcat(range(1,length(M[1,:]);step=1),transpose(M))
    scatter!(centroids[:,2],centroids[:,3], group=centroids[:,1],color_palette=colors,border=:none,axis=nothing,markershape= :star,markersize = 7,markercolor=:match,legend= false, title="Iteration $i")
    savefig(plotfig,"plots\\kmeans-it$i.pdf")
end
