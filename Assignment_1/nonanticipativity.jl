
# Create and populate the “non-anticipativity” sets by iterating through all
# scenario pairs and, for each pair, checking up to which stage they share history.
function create_non_anticipativity_sets(number_of_sim_periods)
    non_anticipativity_sets = []
    for s1 in 1:number_of_sim_periods
        for s2 in 1:number_of_sim_periods
            shared_history = 0
            for t in 1:number_of_sim_periods
                if s1[t] == s2[t]
                    shared_history = t
                else
                    break
                end
            end
            push!(non_anticipativity_sets, (s1, s2, shared_history))
        end
    end
    return non_anticipativity_sets
end