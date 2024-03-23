#Import packages
using JuMP
using Gurobi
using Printf
using Distributions
using Random

# Set seed for reproducibility
Random.seed!(1234)

# Load data and functions from other files
include("V2_Assignment_A_codes/V2_02435_two_stage_problem_data.jl")
include("V2_Assignment_A_codes/V2_simulation_experiments.jl")
include("V2_Assignment_A_codes/V2_price_process.jl")

function calculate_expected_price(current_prices, number_of_warehouses, T_lookahead, S_samples)
    # Initialize the output matrix
    price_scenarios = zeros(Float64, number_of_warehouses, T_lookahead, S_samples)
    expected_price = zeros(Float64, number_of_warehouses, T_lookahead)
    # Depending on the number of stages, fill the price scenarios matrix accordingly.
    # Always use the price from the previous stage to sample the next price
    sqrt_S = ceil(Int, sqrt(S_samples))
    if T_lookahead == 5 
        for w in 1:number_of_warehouses
            for s in 1:S_samples
                for t in 1:T_lookahead
                    if t == 1
                        price_scenarios[w, t, s] = current_prices[w] # Stage 1 price is the current price for all scenarios
                    elseif t == 2
                        if s <= ceil(Int, sqrt(sqrt(sqrt_S)))
                            price_scenarios[w, t, s] = sample_next(current_prices[w]) # Sample a new price for the first sqrt(sqrt(sqrt(S))) scenarios
                        else
                            price_scenarios[w, t, s] = price_scenarios[w, t, s-ceil(Int, sqrt(sqrt(sqrt_S)))] # Use the same price as the previous sqrt(sqrt(sqrt(S))) scenarios
                        end 
                    elseif t == 3
                        if s <= ceil(Int, sqrt(sqrt_S))
                            price_scenarios[w, t, s] = sample_next(price_scenarios[w, t-1, s]) # Sample a new price for the first sqrt(sqrt(S)) scenarios
                        else
                            price_scenarios[w, t, s] = price_scenarios[w, t, s-ceil(Int, sqrt(sqrt_S))] # Use the same price as the previous sqrt(sqrt(S)) scenarios
                        end 
                    elseif t == 4
                        if s <= sqrt_S
                            price_scenarios[w, t, s] = sample_next(price_scenarios[w, t-1, s]) # Sample a new price for the first sqrt(S) scenarios
                        else
                            price_scenarios[w, t, s] = price_scenarios[w, t, s-sqrt_S] # Use the same price as the previous sqrt(S) scenarios
                        end 
                    else
                        price_scenarios[w, t, s] = sample_next(price_scenarios[w, t-1, s]) # Sample a new unique price for all scenarios
                    end
                end
            end
        end
    elseif T_lookahead == 4
        for w in 1:number_of_warehouses
            for s in 1:S_samples
                for t in 1:T_lookahead
                    if t == 1
                        price_scenarios[w, t, s] = current_prices[w] 
                    elseif t == 2
                        if s <= ceil(Int, sqrt(sqrt_S))
                            price_scenarios[w, t, s] = sample_next(current_prices[w])
                        else
                            price_scenarios[w, t, s] = price_scenarios[w, t, s-ceil(Int,sqrt(sqrt_S))]
                        end 
                    elseif t == 3
                        if s <= sqrt_S
                            price_scenarios[w, t, s] = sample_next(price_scenarios[w, t-1, s])
                        else
                            price_scenarios[w, t, s] = price_scenarios[w, t, s-sqrt_S]
                        end 
                    else
                        price_scenarios[w, t, s] = sample_next(price_scenarios[w, t-1, s]) 
                    end
                end
            end
        end
    elseif T_lookahead == 3
        for w in 1:number_of_warehouses
            for s in 1:S_samples
                for t in 1:T_lookahead
                    if t == 1
                        price_scenarios[w, t, s] = current_prices[w] 
                    elseif t == 2
                        if s <= sqrt_S
                            price_scenarios[w, t, s] = sample_next(current_prices[w])
                        else
                            price_scenarios[w, t, s] = price_scenarios[w, t, s-sqrt_S]
                        end 
                    else
                        price_scenarios[w, t, s] = sample_next(price_scenarios[w, t-1, s])
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

    for w in 1:number_of_warehouses
        for t in 1:T_lookahead
            expected_price[w,t] = sum(price_scenarios[w,t,s] for s in 1:S_samples) / S_samples
        end
    end
        
    return expected_price
end

function make_EV_here_and_now_decision(number_of_sim_periods, tau, current_stock, current_prices, current_demands)
    S_samples = 1500
    T_lookahead = min(number_of_sim_periods-tau+1, 3)
    W = collect(1:number_of_warehouses)
    T = collect(1:T_lookahead)
    
    expected_price = calculate_expected_price(current_prices, number_of_warehouses, T_lookahead, S_samples)

    # Define the model
    model_ev = Model(Gurobi.Optimizer)

    # Variables
    @variable(model_ev, 0 <= x_wt[w in W, t in T]) # Coffee ordered at warehouse w in period t
    @variable(model_ev, 0 <= z_wt[w in W, t in T]) # Coffee stored at warehouse w in period t
    @variable(model_ev, 0 <= y_send_wqt[w in W, q in W, t in T]) # Coffee sent from warehouse w to warehouse q in period t
    @variable(model_ev, 0 <= y_receive_wqt[w in W, q in W, t in T]) # Coffee received at warehouse w from warehouse q in period t
    @variable(model_ev, 0 <= m_wt[w in W, t in T]) # Missing demand at warehouse w in period t
   
    @objective(model_ev, Min, 
        sum(expected_price[w,t] * x_wt[w,t] for w in W, t in T) +
        sum(cost_miss[w] * m_wt[w,t] for w in W, t in T) +
        sum(cost_tr[w,q] * y_send_wqt[w,q,t] for w in W, q in W, t in T))

    # Define the constraints
    # Demand hour 1
    @constraint(model_ev, demand_t0[w in W],
                x_wt[w,1] - z_wt[w,1] + current_stock[w] + sum(y_receive_wqt[w,q,1] for q in W) - sum(y_send_wqt[w,q,1] for q in W) + m_wt[w,1] == current_demands[w])

    # Demand hour t
    @constraint(model_ev, demand_t[w in W, t in T[2:end]],
                x_wt[w,t] - z_wt[w,t] + z_wt[w,t-1] + sum(y_receive_wqt[w,q,t] for q in W) - sum(y_send_wqt[w,q,t] for q in W) + m_wt[w,t] == current_demands[w])

    # Storage capacity
    @constraint(model_ev, storage_capacity[w in W, t in T],
                z_wt[w,t] <= warehouse_capacities[w])

    # Transportation capacity
    @constraint(model_ev, transportation_capacity[w in W, q in W, t in T],
                y_send_wqt[w,q,t] <= transport_capacities[w,q])

    # Constraint saying the amount sent from w to q is the same as received at q from w
    @constraint(model_ev, send_receive[w in W, q in W, t in T],
                y_send_wqt[w,q,t] == y_receive_wqt[q,w,t])

    # Storage one day before transfer, initial stock in t=1
    @constraint(model_ev, storage_t1[w in W, q in W],
                y_send_wqt[w,q,1] <= current_stock[w])

    @constraint(model_ev, storage_before_transfer[w in W, q in W, t in T[2:end]],
                y_send_wqt[w,q,t] <= z_wt[w,t-1])

    # Solve the model
    optimize!(model_ev)

    # Return the optimal solution
    if termination_status(model_ev) == MOI.OPTIMAL
        println("Optimal solution found")
        system_cost = objective_value(model_ev)
        println("Total system cost: ", system_cost)
        println("Number of decision variables: ", num_variables(model_ev))
        # Store here and now decisions
        x = value.(x_wt[:,1])
        send = value.(y_send_wqt[:,:,1])
        receive = value.(y_receive_wqt[:,:,1])
        z = value.(z_wt[:,1])
        m = value.(m_wt[:,1])

        # Return the relevant results
        return x, send, receive, z, m
    else
        println("No solution found")
        return nothing
    end 
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
#     x, send, receive, z, m = make_EV_here_and_now_decision(number_of_sim_periods, tau, current_stock, current_prices, current_demands)
# end 
