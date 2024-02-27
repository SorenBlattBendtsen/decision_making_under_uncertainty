#1e-Experiments of all Methods
# import packages
using JuMP
using Gurobi
using Printf
using Plots
using Random

# random seed for reproducability
Random.seed!(1234)

# Load data from 02435_two_stage_problem_data.jl function load_the_data()
include("V2_Assignment_A_codes/V2_02435_two_stage_problem_data.jl")
number_of_warehouses, W, cost_miss, cost_tr, warehouse_capacities, transport_capacities, initial_stock, number_of_simulation_periods, sim_T, degradation_factor = load_the_data()

# Constant demand for all warehouses and all periods
demand = 4*ones(number_of_warehouses, number_of_simulation_periods)

include("V2_Assignment_A_codes/V2_simulation_experiments.jl")
number_of_experiments, Expers, price_trajectory = simulation_experiments_creation(number_of_warehouses, W, number_of_simulation_periods)

# Generate random data for day 1
price = rand(Uniform(0,10),10,10)

include("V2_Assignment_A_codes/V2_price_process.jl")
# 1000 scenarios, 10 reduced scenarios
S = 1000
num_reduced = 10


#Load the optimality in hindsight program
include("1c_optimality_in_hindsight.jl")
result = Calculate_OiH_solution(price)

system_cost, x_values, z_values, y_send_values, y_receive_values = Calculate_OiH_solution(price)

#Load the Expected Value Program
include("1b Expected Value.jl")
expected_price = calculate_expected_prices(price)
result_EV = make_EV_here_and_now_decision(price)

system_cost, prices_day_1, orders_day_1, storage_day_1, send_receive_day_1,
expected_prices_day_2, orders_day_2, storage_day_2, send_receive_day_2 = result_EV


#Load the Two-Stage Program
include("1d_two_stage.jl")
scen_gen = Scenario_generation(price, S, number_of_warehouses, num_reduced, plots=true)

num_reduced_values = [5, 20, 50]

# Dictionaries to store results for each num_reduced value
system_costs_dict = Dict()
prices_day_1_dict = Dict()
orders_day_1_dict = Dict()
storage_day_1_dict = Dict()
send_receive_day_1_dict = Dict()
prices_day_2_dict = Dict()
orders_day_2_dict = Dict()
storage_day_2_dict = Dict()
send_receive_day_2_dict = Dict()

# Loop through each num_reduced value and call the function
for num_reduced in num_reduced_values
    println("Running for num_reduced = $num_reduced")
    
    result_TS = Make_Stochastic_here_and_now_decision(p_wt, S, num_reduced)
    
    if result !== nothing
        system_cost, prices_day_1, orders_day_1, storage_day_1, send_receive_day_1,
        prices_day_2, orders_day_2, storage_day_2, send_receive_day_2 = result_TS
        
        # Store results 
        system_costs_dict[num_reduced] = system_cost
        prices_day_1_dict[num_reduced] = prices_day_1
        orders_day_1_dict[num_reduced] = orders_day_1
        storage_day_1_dict[num_reduced] = storage_day_1
        send_receive_day_1_dict[num_reduced] = send_receive_day_1
        prices_day_2_dict[num_reduced] = prices_day_2
        orders_day_2_dict[num_reduced] = orders_day_2
        storage_day_2_dict[num_reduced] = storage_day_2
        send_receive_day_2_dict[num_reduced] = send_receive_day_2
        
    else
        println("No solution found for num_reduced = $num_reduced")
    end
end


