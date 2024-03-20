# File: 02435_multistage_policy.jl

using JuMP
using Gurobi
using Distributions
using Random
using Plots
# Set seed for reproducibility
Random.seed!(1234)
# Load data and functions from other files
include("V2_Assignment_A_codes/V2_02435_multistage_problem_data.jl")
include("fast-forward-selection.jl")
include("V2_Assignment_A_codes/V2_price_process.jl")

# Define the function to generate price scenarios
function generate_price_scenarios(current_prices, number_of_warehouses, T_lookahead, S_samples)
    # Initialize the output matrix
    price_scenarios = zeros(Float64, number_of_warehouses, T_lookahead, S_samples)
    
    sqrt_S = ceil(Int, sqrt(S_samples))
    if T_lookahead == 5 
        for w in 1:number_of_warehouses
            for s in 1:S_samples
                for t in 1:T_lookahead
                    if t == 1
                        price_scenarios[w, t, s] = current_prices[w] # Using the first scenario from t = 2
                    elseif t == 2
                        if s <= ceil(Int, sqrt(sqrt(sqrt_S)))
                            price_scenarios[w, t, s] = sample_next(current_prices[w])
                        else
                            price_scenarios[w, t, s] = price_scenarios[w, t, s-ceil(Int, sqrt(sqrt(sqrt_S)))]
                        end 
                        #continue # Scenarios already filled in the previous loop
                    elseif t == 3
                        if s <= ceil(Int, sqrt(sqrt_S))
                            price_scenarios[w, t, s] = sample_next(price_scenarios[w, t-1, s])
                        else
                            price_scenarios[w, t, s] = price_scenarios[w, t, s-ceil(Int, sqrt(sqrt_S))]
                        end 
                    elseif t == 4
                        if s <= sqrt_S
                            price_scenarios[w, t, s] = sample_next(price_scenarios[w, t-1, s])
                        else
                            price_scenarios[w, t, s] = price_scenarios[w, t, s-sqrt_S]
                        end 
                    else
                        price_scenarios[w, t, s] = sample_next(price_scenarios[w, t-1, s]) # Using the same scenarios from t = 3
                    end
                end
            end
        end
    elseif T_lookahead == 4
        for w in 1:number_of_warehouses
            for s in 1:S_samples
                for t in 1:T_lookahead
                    if t == 1
                        price_scenarios[w, t, s] = current_prices[w] # Using the first scenario from t = 2
                    elseif t == 2
                        if s <= ceil(Int, sqrt(sqrt_S))
                            price_scenarios[w, t, s] = sample_next(current_prices[w])
                        else
                            price_scenarios[w, t, s] = price_scenarios[w, t, s-ceil(Int,sqrt(sqrt_S))]
                        end 
                        #continue # Scenarios already filled in the previous loop
                    elseif t == 3
                        if s <= sqrt_S
                            price_scenarios[w, t, s] = sample_next(price_scenarios[w, t-1, s])
                        else
                            price_scenarios[w, t, s] = price_scenarios[w, t, s-sqrt_S]
                        end 
                        #continue # Scenarios already filled in the previous loop
                    else
                        price_scenarios[w, t, s] = sample_next(price_scenarios[w, t-1, s]) # Using the same scenarios from t = 3
                    end
                end
            end
        end
    elseif T_lookahead == 3
        for w in 1:number_of_warehouses
            for s in 1:S_samples
                for t in 1:T_lookahead
                    if t == 1
                        price_scenarios[w, t, s] = current_prices[w] # Using the first scenario from t = 2
                    elseif t == 2
                        if s <= sqrt_S
                            price_scenarios[w, t, s] = sample_next(current_prices[w])
                        else
                            price_scenarios[w, t, s] = price_scenarios[w, t, s-sqrt_S]
                        end 
                        #continue # Scenarios already filled in the previous loop
                    else
                        price_scenarios[w, t, s] = sample_next(price_scenarios[w, t-1, s]) # Using the same scenarios from t = 3
                    end
                end
            end
        end
    elseif T_lookahead == 2
        for w in 1:number_of_warehouses
            for s in 1:S_samples
                for t in 1:T_lookahead
                    if t == 1
                        price_scenarios[w, t, s] = current_prices[w] 
                    else
                        price_scenarios[w, t, s] = sample_next(current_prices[w])
                    end
                end
            end
        end
    else # T == 1
        for w in 1:number_of_warehouses
            for s in 1:S_samples
                for t in 1:T_lookahead
                    price_scenarios[w, t, s] = current_prices[w] 
                end
            end
        end
    end
        
    return price_scenarios
end

# Define the function to reduce the number of scenarios using FFS
function reduce_scenarios(price_scenarios, S_reduced, T_lookahead)
    # Calculate distance matrix (euclidean distance) for scenario reduction
    B = size(price_scenarios, 3)
    D = zeros(Float64, B, B)
    for i = 1:B
        for j = 1:B
            D[i,j] = sqrt(sum((price_scenarios[l,end,i]-price_scenarios[l,end,j])*(price_scenarios[l,end,i]-price_scenarios[l,end,j])  for l = 1:number_of_warehouses))
        end
    end

    # Initialize equiprobable probabilities
    probabilities = repeat([1.0/B], 1, B)[1,:]
    # Use fast forward selection and apply it to get the reduced scenarios and updated probabilities
    result = FastForwardSelection(D, probabilities, S_reduced)
    # Resulting probabilities
    reduced_probabilities = result[1]
    # Selected scenario indices
    reduced_scenario_indices = result[2]
    # Filter price_scenarios based on reduced_scenario_indices
    reduced_price_scenarios = zeros(Float64, number_of_warehouses, T_lookahead, S_reduced)
    
    for i = 1:S_reduced
        reduced_price_scenarios[:,:,i] = price_scenarios[:,:,reduced_scenario_indices[i]]
    end

    return reduced_price_scenarios, reduced_probabilities, reduced_scenario_indices
end

# Define the function to find the indices of the same price scenarios per stage, used for non-anticipativity constraints
function same_price_indices_per_stage(reduced_price_scenarios, T_lookahead, S_reduced)
    # Initialize the output matrix
    same_price_indices = zeros(Int, T_lookahead, S_reduced)
    for t in 1:T_lookahead
        for s in 1:S_reduced
            same_price_indices[t, s] = findall(x -> x == reduced_price_scenarios[1, t, s], reduced_price_scenarios[1, t, :])[1]
        end
    end
    return same_price_indices
end


# Define the function to plot the reduced price scenarios
function plot_reduced_scenarios(price_scenarios, reduced_scenario_indices, num_reduced, number_of_warehouses, W)
    p_wt_scenarios_t = transpose(price_scenarios[:,3,:])
    plot(price_scenarios[:,3,:], xlabel="Warehouse", ylabel="Price", title="Reduced price scenarios", color=:lightgray, legend=false, alpha=0.8, xticks=(1:number_of_warehouses, W))
    # plot reduced scenarios
    reduced_data = zeros(Float64, number_of_warehouses, num_reduced)
    for i = 1:num_reduced
        reduced_data[:, i] = price_scenarios[:,3,reduced_scenario_indices[i]]
    end
    plot!(reduced_data, legend=false, color=:auto)
end 

# Define the function to make a multistage here and now decision
function make_multistage_here_and_now_decision(number_of_sim_periods, tau, current_stock, current_prices)
    # Get relevant input data
    S_samples = 1500
    S_reduced = 45
    T_lookahead = min(number_of_sim_periods-tau+1, 5)
    W = collect(1:number_of_warehouses)
    #T = collect(tau:tau+T_lookahead-1)
    T = collect(1:T_lookahead)
    S = collect(1:S_reduced)
    price_scenarios = generate_price_scenarios(current_prices, number_of_warehouses, T_lookahead, S_samples)
    reduced_price_scenarios, reduced_probabilities, reduced_scenario_indices = reduce_scenarios(price_scenarios, S_reduced, T_lookahead)
    same_price_indices = same_price_indices_per_stage(reduced_price_scenarios, T_lookahead, S_reduced)
    # Define the optimization model
    model_2 = Model(Gurobi.Optimizer)
    # variables
    @variable(model_2, 0 <= x_wts[w in W, t in T, s in S]) # Coffee ordered at warehouse w
    @variable(model_2, 0 <= z_wts[w in W, t in T, s in S]) # Coffee stored at warehouse w
    @variable(model_2, 0 <= y_send_wqts[w in W, q in W, t in T, s in S]) # Coffee sent from warehouse w to warehouse q
    @variable(model_2, 0 <= y_receive_wqts[w in W, q in W, t in T, s in S]) # Coffee received at warehouse w from warehouse q
    @variable(model_2, 0 <= m_wts[w in W, t in T, s in S]) # Missing demand at warehouse w

    @objective(model_2, Min, 
                            sum(reduced_probabilities[s] 
                            * (sum(cost_tr[w,q]*y_send_wqts[w,q,t,s] for w in W, q in W, t in T) 
                            + sum(cost_miss[w]*m_wts[w,t,s] for w in W, t in T)
                            + sum(reduced_price_scenarios[w,t,s]*x_wts[w,t,s] for w in W, t in T))
                            for s in S))

    # constraints
    # Demand balance
    @constraint(model_2, demand_1[w in W, s in S],
                x_wts[w,1,s] - z_wts[w,1,s] + current_stock[w]
                + sum(y_receive_wqts[w,q,1,s] for q in W)
                - sum(y_send_wqts[w,q,1,s] for q in W)
                + m_wts[w,1,s] == current_demands[w])

    @constraint(model_2, demand_t[w in W, t in T[2:end], s in S],
                x_wts[w,t,s] - z_wts[w,t,s] + z_wts[w,t-1,s]
                + sum(y_receive_wqts[w,q,t,s] for q in W)
                - sum(y_send_wqts[w,q,t,s] for q in W)
                + m_wts[w,t,s] == current_demands[w])

    # Storage capacity
    @constraint(model_2, storage_capacity[w in W, t in T, s in S],
                z_wts[w,t,s] <= warehouse_capacities[w])

    # Transportation capacity
    @constraint(model_2, transportation_capacity[w in W, q in W, t in T, s in S],
                y_send_wqts[w,q,t,s] <= transport_capacities[w,q])

    # Constraint saying the amount sent from w to q is the same as received at q from w
    @constraint(model_2, send_receive_1[w in W, q in W, t in T, s in S], 
                y_send_wqts[w,q,t,s] == y_receive_wqts[q,w,t,s])

    #Storage one day before transfer, initial stock in t=1
    @constraint(model_2, storage_1[w in W, q in W, s in S],
                y_send_wqts[w,q,1,s] <= current_stock[w])

    @constraint(model_2, storage_t[w in W, q in W, t in T[2:end], s in S],
                y_send_wqts[w,q,t,s] <= z_wts[w,t-1,s])

    # Non anticipativity constraints for stage all stages
    @constraint(model_2, non_anticipativity_x[w in W, t in T, s in S],
                x_wts[w,t,s] == x_wts[w,t,same_price_indices[t,s]])

    @constraint(model_2, non_anticipativity_z[w in W, t in T, s in S],
                z_wts[w,t,s] == z_wts[w,t,same_price_indices[t,s]])

    @constraint(model_2, non_anticipativity_y_send[w in W, q in W, t in T, s in S],
                y_send_wqts[w,:,t,s] .== y_send_wqts[w,:,t,same_price_indices[t,s]])

    @constraint(model_2, non_anticipativity_y_receive[w in W, q in W, t in T, s in S],
                y_receive_wqts[w,:,t,s] .== y_receive_wqts[w,:,t,same_price_indices[t,s]])

    @constraint(model_2, non_anticipativity_m[w in W, t in T, s in S],
                m_wts[w,t,s] == m_wts[w,t,same_price_indices[t,s]])


    optimize!(model_2)

    # Print the objective value and the optimal solution for the first stage variables
    if termination_status(model_2) == MOI.OPTIMAL
        println("Optimal solution found")
        #println("Time period: ", tau)
        println("Total number of decision variables: ", num_variables(model_2))
        println("Objective value: ", objective_value(model_2))
        #println("Here and now decision:")
        # Store here and now decisions
        x = value.(x_wts[:,1,1])
        send = value.(y_send_wqts[:,:,1,1])
        receive = value.(y_receive_wqts[:,:,1,1])
        z = value.(z_wts[:,1,1])
        m = value.(m_wts[:,1,1])
        # for w in W
        #     println("Warehouse ", w)
        #     println("x_wts: ", x[w])
        #     println("z_wts: ", z[w])
        #     println("y_send_wqts: ", send[w,:])
        #     println("y_receive_wqts: ", receive[w,:])
        #     println("m_wts: ", m[w])
        # end
    end
    return x, send, receive, z, m
end 

# Example usage
# number_of_warehouses, W, cost_miss, cost_tr, warehouse_capacities, transport_capacities, initial_stock = load_the_data()
# current_prices = zeros(number_of_warehouses)
# for w in W
#     current_prices[w] = round.(rand(Uniform(0,10)), digits=2)
# end
# current_demands = 4*ones(number_of_warehouses)
# tau = 1
# number_of_sim_periods = 5
# current_stock = initial_stock

# x = Dict()
# send = Dict()
# receive = Dict()
# z = Dict()
# m = Dict()

# @time begin
#     x, send, receive, z, m = make_multistage_here_and_now_decision(number_of_sim_periods, tau, current_stock, current_prices)
# end 