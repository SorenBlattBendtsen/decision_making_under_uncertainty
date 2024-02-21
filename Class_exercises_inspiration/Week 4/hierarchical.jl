

include("create_dataset.jl")

using Clustering

D = zeros(Float64, length(dataset_x),length(dataset_x))
for i = 1:length(dataset_x)
    for j = 1:length(dataset_x)
        D[i,j] = sqrt(((dataset_x[i]-dataset_x[j])*(dataset_x[i]-dataset_x[j])) + ((dataset_y[i]-dataset_y[j])*(dataset_y[i]-dataset_y[j])) )
    end
end


result = hclust(D; linkage=:average)


c = length(result.heights[:,1])
for i in 1: c

    a = cutree(result; k=i)

    affected = []
    pos = []
    if c-i+2 <= length(result.merges[:,1])
        pos = result.merges[c-i+2,:]
    end
    while length(pos) > 0
        first = pos[1]
        deleteat!(pos, 1)
        if first < 0
            append!(affected, first)
        else
            for elem in result.merges[first,:]
                if elem >=0
                    append!(pos, elem)
                else
                    append!(affected, elem)
                end
            end
        end
    end
    println(affected)
    clustering = hcat(a,hcat(dataset_x, dataset_y))
    plotfig = scatter(clustering[:,2], clustering[:,3],group=clustering[:,1],border=:none,axis=nothing,markershape= :circle,markercolor=:auto,legend= false, title="$i -> $(i-1) clusters")

    affected_points = zeros(Float64, length(affected), 2)

    for i in 1:length(affected)
        id=-1*affected[i]
        affected_points[i,1] = dataset_x[id]
        affected_points[i,2] = dataset_y[id]
    end
    scatter!(affected_points[:,1], affected_points[:,2], markershape=:xcross, markerfacecolor="none",color = :black)

    #plotfig = scatter(clustering[:,2], clustering[:,3],group=clustering[:,1],border=:none,axis=nothing,markershape= :auto,markercolor=:auto,legend= false, title="$i clusters")
    savefig(plotfig,"plots\\hclust-$i.pdf")
end




elbowlimfig = plot(reverse(result.heights), xlim=(1,10), legend = false, ylabel= "Distance of merge (based on linkage)", xlabel="# cluster")
savefig(elbowlimfig,"plots\\hclust-elbowlim.pdf")

elbowfig = plot(reverse(result.heights), legend=false, ylabel= "Distance of merge (based on linkage)", xlabel="# clusters")
savefig(elbowfig,"plots\\hclust-elbow.pdf")


assignment = cutree(result; k=4)
clustering_4 = hcat(assignment,hcat(dataset_x, dataset_y))
plot4fig = scatter(clustering_4[:,2], clustering_4[:,3],group=clustering_4[:,1],border=:none,axis=nothing,markershape= :circle,markercolor=:auto,legend= false, title="4 cluster")
savefig(plot4fig,"plots\\hclust-4cluster.pdf")
