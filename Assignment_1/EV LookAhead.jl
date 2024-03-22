#Import packages
using JuMP
using Gurobi
using Printf
using Distributions
using Random
using BenchmarkTools

# Set seed for reproducibility
Random.seed!(1234)

# Load data from 02435_two_stage_problem_data.jl function load_the_data()
include("V2_Assignment_A_codes/V2_02435_two_stage_problem_data.jl")
number_of_warehouses, W, cost_miss, cost_tr, warehouse_capacities, transport_capacities, initial_stock, number_of_simulation_periods, sim_T, degradation_factor = load_the_data()

# Generate random data for day 1
include("V2_Assignment_A_codes/V2_simulation_experiments.jl")
number_of_experiments, Expers, price_trajectory = simulation_experiments_creation(number_of_warehouses, W, number_of_simulation_periods)
current_prices = price_trajectory[1,:,1]
demand = 4

# Include the file containing the price process function
include("V2_Assignment_A_codes/V2_price_process.jl")


# Function that takes the inital price and calculates the expected price for 1000 samples using sample_next
function calculate_expected_prices(current_prices,number_of_sim_periods)
    num_samples = if number_of_sim_periods == 3
        1500
    elseif number_of_sim_periods in [4, 5]
        1500
    elseif number_of_sim_periods ==2
        1200
    else
        error("Unsupported simulation period")
    end

    expected_price = zeros(number_of_warehouses)  # Initialize expected prices matrix
    for w in 1:number_of_warehouses
        initial_prices = current_prices[w]  # Extract initial prices as a scalar or vector
        expected_price_tmp = 0.0  # Initialize as a scalar
        for i in 1:num_samples
            price_sample = sample_next(initial_prices)  # Use initial prices as input to sample_next
            expected_price_tmp += price_sample
        end
        expected_price_tmp /= num_samples  # Average over samples
        expected_price[w] = expected_price_tmp  # Assign expected prices back to the main matrix
    end
    return expected_price
end

function make_EV_here_and_now_decision(number_of_sim_periods, tau, current_stock, current_prices)
    expected_price = calculate_expected_prices(current_prices, number_of_sim_periods)

    # Define the model
    model_ev = Model(Gurobi.Optimizer)

    # Variables
    @variable(model_ev, 0 <= x_wt[w in W, t in 1:number_of_sim_periods]) # Coffee ordered at warehouse w in period t
    @variable(model_ev, 0 <= z_wt[w in W, t in 1:number_of_sim_periods]) # Coffee stored at warehouse w in period t
    @variable(model_ev, 0 <= y_send_wqt[w in W, q in W, t in 1:number_of_sim_periods]) # Coffee sent from warehouse w to warehouse q in period t
    @variable(model_ev, 0 <= y_receive_wqt[w in W, q in W, t in 1:number_of_sim_periods]) # Coffee received at warehouse w from warehouse q in period t
    @variable(model_ev, 0 <= m_wt[w in W, t in 1:number_of_sim_periods]) # Missing demand at warehouse w in period t

    println("Number of decision variables: ", num_variables(model_ev))
   
    @objective(model_ev, Min, 
        sum(current_prices[w] * x_wt[w,t] for w in W, t in 1:tau) +
        sum(cost_miss[w] * m_wt[w,t] for w in W, t in 1:tau) +
        sum(cost_tr[w,q] * y_send_wqt[w,q,t] for w in W, q in W, t in 1:tau) +
        sum(expected_price[w] * x_wt[w,t] for w in W, t in (tau + 1):number_of_sim_periods[end]) +
        sum(cost_miss[w] * m_wt[w,t] for w in W, t in (tau + 1):number_of_sim_periods[end]) +
        sum(cost_tr[w,q] * y_send_wqt[w,q,t] for w in W, q in W, t in (tau + 1):number_of_sim_periods[end]))

    # Define the constraints
    # Demand hour 1
    @constraint(model_ev, demand_t0[w in W],
                x_wt[w,1] - z_wt[w,1] + current_stock[w] + sum(y_receive_wqt[w,q,1] for q in W) - sum(y_send_wqt[w,q,1] for q in W) + m_wt[w,1] == demand)

    # Demand hour t
    @constraint(model_ev, demand_t[w in W, t in 2:number_of_sim_periods[end]],
                x_wt[w,t] - z_wt[w,t] + z_wt[w,t-1] + sum(y_receive_wqt[w,q,t] for q in W) - sum(y_send_wqt[w,q,t] for q in W) + m_wt[w,t] == demand)

    # Storage capacity
    @constraint(model_ev, storage_capacity[w in W, t in 1:number_of_sim_periods],
                z_wt[w,t] <= warehouse_capacities[w])

    # Transportation capacity
    @constraint(model_ev, transportation_capacity[w in W, q in W, t in 1:number_of_sim_periods],
                y_send_wqt[w,q,t] <= transport_capacities[w,q])

    # Constraint saying the amount sent from w to q is the same as received at q from w
    @constraint(model_ev, send_receive[w in W, q in W, t in 1:number_of_sim_periods],
                y_send_wqt[w,q,t] == y_receive_wqt[q,w,t])

    # Storage one day before transfer, initial stock in t=1
    @constraint(model_ev, storage_t1[w in W, q in W],
                y_send_wqt[w,q,1] <= current_stock[w])

    @constraint(model_ev, storage_before_transfer[w in W, q in W, t in 2:number_of_sim_periods[end]],
                y_send_wqt[w,q,t] <= z_wt[w,t-1])

    # Solve the model
    optimize!(model_ev)

    # Return the optimal solution
    if termination_status(model_ev) == MOI.OPTIMAL
        println("Optimal solution found")
        system_cost = objective_value(model_ev)
        println("Total system cost: ", system_cost)
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

# sim_periods = [2, 3, 4, 5]
# for number_of_sim_periods in sim_periods
#     println("Simulation Period: ", number_of_sim_periods)
#     current_stock = initial_stock
#     make_EV_here_and_now_decision(number_of_sim_periods, current_stock, current_prices)
#     # x = Dict()
#     # send = Dict()
#     # receive = Dict()
#     # z = Dict()
#     # m = Dict()

#     # @time begin
#     #     x, send, receive, z, m = make_multistage_here_and_now_decision(number_of_sim_periods, tau, current_stock, current_prices)
#     # end 
#     println("----------------------------------------")
# end

# Example usage

tau = 1
number_of_sim_periods = 5
current_stock = initial_stock

x = Dict()
send = Dict()
receive = Dict()
z = Dict()
m = Dict()

@time begin
    x, send, receive, z, m = make_EV_here_and_now_decision(number_of_sim_periods, tau, current_stock, current_prices)
end 