using JuMP
using Gurobi
using Distributions
using Random

include("V2_02435_multistage_problem_data.jl")
include("fast-forward-selection.jl")
include("V2_price_process.jl")
include("nonanticipativity.jl")

function generate_and_reduce_scenarios(current_prices, num_samples, reduced_samples, lookahead)
    all_scenarios = Array{Float64,3}(undef, num_samples, length(current_prices), lookahead)
    reduced_price_scenarios = Array{Float64,3}(undef, reduced_samples, length(current_prices), lookahead)
    reduced_probabilities = zeros(reduced_samples, lookahead)
    
    for t in 1:lookahead
        for s in 1:num_samples
            for w in 1:length(current_prices)
                all_scenarios[s,w,t] = sample_next(current_prices[w])
            end
        end

        D = zeros(num_samples, num_samples)
        for i in 1:num_samples
            for j in 1:num_samples
                D[i,j] = sqrt(sum((all_scenarios[i,:,t] - all_scenarios[j,:,t]).^2))
            end
        end

        probabilities, selected_indices = FastForwardSelection(D, fill(1.0 / num_samples, num_samples), reduced_samples)

        for i in 1:reduced_samples
            reduced_price_scenarios[i,:,t] = all_scenarios[selected_indices[i],:,t]
            reduced_probabilities[i,t] = probabilities[i]
        end
    end
    return reduced_price_scenarios, reduced_probabilities
end

function make_multistage_here_and_now_decision(number_of_sim_periods, tau, current_stock, current_prices)
    number_of_warehouses, W, cost_miss, cost_tr, warehouse_capacities, transport_capacities, initial_stock, number_of_simulation_periods, sim_T, demand_trajectory = load_the_data()

    lookahead = min(5-tau+1, 3)
    num_samples = 100
    reduced_samples = 10
    price_scenarios, probabilities = generate_and_reduce_scenarios(current_prices, num_samples, reduced_samples, lookahead)

    model = Model(Gurobi.Optimizer)
    @variable(model, 0 <= x_wt[w in W, t in 1:lookahead, s in 1:reduced_samples])
    @variable(model, 0 <= z_wt[w in W, t in 1:lookahead, s in 1:reduced_samples])
    @variable(model, 0 <= y_send_wqt[w in W, q in W, t in 1:lookahead, s in 1:reduced_samples])
    @variable(model, 0 <= y_receive_wqt[w in W, q in W, t in 1:lookahead, s in 1:reduced_samples])
    @variable(model, 0 <= m_wt[w in W, t in 1:lookahead, s in 1:reduced_samples])

    non_anticipativity_sets = create_non_anticipativity_sets(price_scenarios)

    for t in 1:lookahead
        for s in 1:reduced_samples
            @constraint(model, demand_balance[w in W, t, s],
                        x_wt[w,t,s] - z_wt[w,t,s] + (t == 1 ? initial_stock[w] : z_wt[w,t-1,s])
                        + sum(y_receive_wqt[w,q,t,s] for q in W)
                        - sum(y_send_wqt[w,q,t,s] for q in W)
                        + m_wt[w,t,s] == demand_trajectory[w,t])
    
            @constraint(model, storage_capacity[w in W, t, s],
                        z_wt[w,t,s] <= warehouse_capacities[w])
    
            @constraint(model, transportation_capacity[w in W, q in W, t, s],
                        y_send_wqt[w,q,t,s] <= transport_capacities[w,q])
    
            @constraint(model, send_receive_consistency[w in W, q in W, t, s],
                        y_send_wqt[w,q,t,s] == y_receive_wqt[q,w,t,s])
    
            @constraint(model, storage_before_transfer[w in W, q in W, t, s],
                        y_send_wqt[w,q,t,s] <= (t == 1 ? initial_stock[w] : z_wt[w,t-1,s]))
        end
    end

    for (i, j, shared_history) in non_anticipativity_sets
        for t in 1:shared_history
            @constraint(model, [w in W], x_wt[w, t, i] == x_wt[w, t, j])
            @constraint(model, [w in W], z_wt[w, t, i] == z_wt[w, t, j])
            @constraint(model, [w in W, q in W], y_send_wqt[w, q, t, i] == y_send_wqt[w, q, t, j])
            @constraint(model, [w in W, q in W], y_receive_wqt[w, q, t, i] == y_receive_wqt[w, q, t, j])
        end
    end

    @objective(model, Min, sum(probabilities[t,s]*(
        sum(current_prices[w] * x_wt[w,t,s] for w in W) 
        + sum(cost_miss[w] * m_wt[w,t,s] for w in W) 
        + sum(cost_tr[w,q] * y_send_wqt[w,q,t,s] for w in W, q in W)) 
        for s in 1:reduced_samples, t in 1:lookahead))

    optimize!(model)

    if termination_status(model) == MOI.OPTIMAL
        x_decision = value.(x_wt[:,1,:])
        z_decision = value.(z_wt[:,1,:])
        send_decision = value.(y_send_wqt[:,:,1,:])
        receive_decision = value.(y_receive_wqt[:,:,1,:])
        m_decision = value.(m_wt[:,1,:])

        return x_decision, send_decision, receive_decision, z_decision, m_decision
    else
        println("No solution found")
        return nothing
    end
end
