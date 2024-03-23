
# Create and populate the “non-anticipativity” sets by iterating through all
# scenario pairs and, for each pair, checking up to which stage they share history.
function create_non_anticipativity_sets(scenarios)
    non_anticipativity_sets = []
    num_scenarios = length(scenarios)
    for i in 1:num_scenarios
        for j in i+1:num_scenarios
            shared_history = 0
            min_length = min(length(scenarios[i]), length(scenarios[j]))
            for t in 1:min_length
                if scenarios[i][t] == scenarios[j][t]
                    shared_history = t
                else
                    break
                end
            end
            push!(non_anticipativity_sets, (i, j, shared_history))
        end
    end
    return non_anticipativity_sets
end


# Example scenarios
scenarios = [[1.2, 2.3, 3.4], [1.2, 2.3, 4.5], [5.6, 6.7, 7.8]]

# Call the function
sets = create_non_anticipativity_sets(scenarios)

# Expected output should indicate:
# - Scenarios 1 and 2 share history up to period 2.
# - Scenarios 1 and 3 have no shared history.
# - Scenarios 2 and 3 have no shared history.
# Print the sets
println(sets)
