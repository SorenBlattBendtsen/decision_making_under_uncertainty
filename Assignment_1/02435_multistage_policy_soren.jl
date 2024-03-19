# File: 02435_multistage_policy.jl

using JuMP
using Gurobi
using Distributions
using Random
using Plots
# Make sure this script is saved in the same directory as your other .jl files
include("V2_Assignment_A_codes/V2_02435_multistage_problem_data.jl")
include("fast-forward-selection.jl")
include("V2_Assignment_A_codes/V2_price_process.jl")



function generate_price_scenarios(current_prices, W, T, S)
    # Initialize the output matrix
    price_scenarios = zeros(W, T, S)
    
    sqrt_S = ceil(Int, sqrt(S))
    if T == 3
        for w in 1:W
            for s in 1:S
                for t in 1:T
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
    elseif T == 2
        for w in 1:W
            for s in 1:S
                for t in 1:T
                    if t == 1
                        price_scenarios[w, t, s] = current_prices[w] 
                    else
                        price_scenarios[w, t, s] = sample_next(current_prices[w])
                    end
                end
            end
        end
    else # T == 1
        for w in 1:W
            for s in 1:S
                for t in 1:T
                    price_scenarios[w, t, s] = current_prices[w] 
                end
            end
        end
    end
        
    return price_scenarios
end

# Example usage:
#price_scenarios = generate_price_scenarios(current_prices, number_of_warehouses, lookahead, num_samples)

function reduce_scenarios(price_scenarios, num_reduced)
    # Calculate distance matrix (euclidean distance) for scenario reduction
    S = size(price_scenarios, 3)
    D = zeros(Float64, S, S)
    for i = 1:S
        for j = 1:S
            D[i,j] = sqrt(sum((price_scenarios[l,end,i]-price_scenarios[l,end,j])*(price_scenarios[l,end,i]-price_scenarios[l,end,j])  for l = 1:number_of_warehouses))
        end
    end

    # Initialize equiprobable probabilities
    probabilities = repeat([1.0/S], 1, S)[1,:]
    # Use fast forward selection and apply it to get the reduced scenarios and updated probabilities
    result = FastForwardSelection(D, probabilities, num_reduced)
    # Resulting probabilities
    reduced_probabilities = result[1]
    # Selected scenario indices
    reduced_scenario_indices = result[2]
    # Filter price_scenarios based on reduced_scenario_indices
    reduced_price_scenarios = zeros(Float64, number_of_warehouses, lookahead, num_reduced)
    for i = 1:num_reduced
        reduced_price_scenarios[:,:,i] = price_scenarios[:,:,reduced_scenario_indices[i]]
    end

    return reduced_price_scenarios, reduced_probabilities, reduced_scenario_indices
end
# example usage
num_reduced = 256
#reduced_price_scenarios, reduced_probabilities, reduced_scenario_indices = reduce_scenarios(price_scenarios, num_reduced)

# for 1 warehouse at t=3, plot the original scenarios in gray and the reduced scenarios in color on top
# transpose price_scenarios for plotting
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

function multi_stage_program(current_prices, number_of_warehouses, T_lookahead, S_samples, S_reduced)
    W = collect(1:number_of_warehouses)
    T = collect(1:T_lookahead)
    S = collect(1:S_reduced)
    price_scenarios = generate_price_scenarios(current_prices, number_of_warehouses, T_lookahead, S_samples)
    reduced_price_scenarios, reduced_probabilities, reduced_scenario_indices = reduce_scenarios(price_scenarios, S_reduced)
    model_2 = Model(Gurobi.Optimizer)
    # variables
    @variable(model_2, 0 <= x_wts[w in W, t in T, s in S]) # Coffee ordered at warehouse w
    @variable(model_2, 0 <= z_wts[w in W, t in T, s in S]) # Coffee stored at warehouse w
    @variable(model_2, 0 <= y_send_wqts[w in W, q in W, t in T, s in S]) # Coffee sent from warehouse w to warehouse q
    @variable(model_2, 0 <= y_receive_wqts[w in W, q in W, t in T, s in S]) # Coffee received at warehouse w from warehouse q
    @variable(model_2, 0 <= m_wts[w in W, t in T, s in S]) # Missing demand at warehouse w
    # Count and print total number of decision variables
    println("Total number of decision variables: ", num_variables(model_2))

    @objective(model_2, Min, 
                            sum(reduced_probabilities[s] 
                            * (sum(cost_tr[w,q]*y_send_wqts[w,q,t,s] for w in W, q in W, t in T) 
                            + sum(cost_miss[w]*m_wts[w,t,s] for w in W, t in T)
                            + sum(reduced_price_scenarios[w,t,s]*x_wts[w,t,s] for w in W, t in T))
                            for s in S))

    # constraints
    # Demand balance
    @constraint(model_2, demand_1[w in W, s in S],
                x_wts[w,1,s] - z_wts[w,1,s] + initial_stock[w]
                + sum(y_receive_wqts[w,q,1,s] for q in W)
                - sum(y_send_wqts[w,q,1,s] for q in W)
                + m_wts[w,1,s] == demand[w,1])

    @constraint(model_2, demand_t[w in W, t in T[2:end], s in S],
                x_wts[w,t,s] - z_wts[w,t,s] + z_wts[w,t-1,s]
                + sum(y_receive_wqts[w,q,t,s] for q in W)
                - sum(y_send_wqts[w,q,t,s] for q in W)
                + m_wts[w,t,s] == demand[w,t])

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
                y_send_wqts[w,q,1,s] <= initial_stock[w])

    @constraint(model_2, storage_t[w in W, q in W, t in T[2:end], s in S],
                y_send_wqts[w,q,t,s] <= z_wts[w,t-1,s])

    optimize!(model_2)

    # Print the objective value and the optimal solution for the first stage variables
    if termination_status(model_2) == MOI.OPTIMAL
        #println("Solve time: ", MOI.get(model_2, MOI.SolveTime()))
        println("Optimal solution found")
        println("Objective value: ", objective_value(model_2))
        println("Period 1:")
        for w in W
            println("Warehouse ", w)
            println("x_wts: ", mean(value.(x_wts[w,1,:])))
            println("z_wts: ", mean(value.(z_wts[w,1,:])))
            println("y_send_wqts: ", mean(value.(y_send_wqts[w,:,1,:])))
            println("y_receive_wqts: ", mean(value.(y_receive_wqts[w,:,1,:])))
            println("m_wts: ", mean(value.(m_wts[w,1,:])))
        end
    end
end 


number_of_warehouses, W, cost_miss, cost_tr, warehouse_capacities, transport_capacities, initial_stock = load_the_data()
current_prices = zeros(number_of_warehouses)
for w in W
    current_prices[w] = round.(rand(Uniform(0,10)), digits=2)
end
tau=1
T_lookahead = min(5-tau+1, 3)
demand = 4*ones(number_of_warehouses, T_lookahead)
S_samples = 2500
S_reduced = 256

multi_stage_program(current_prices, number_of_warehouses, T_lookahead, S_samples, S_reduced)