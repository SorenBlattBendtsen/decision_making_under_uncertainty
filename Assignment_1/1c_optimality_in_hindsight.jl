# import packages
using JuMP
using Gurobi
using Printf

# Load data from 02435_two_stage_problem_data.jl function load_the_data()
include("02435_Assignment_A_codes/02435_two_stage_problem_data.jl")
number_of_warehouses, W, cost_miss, cost_tr, warehouse_capacities, transport_capacities, initial_stock, number_of_simulation_periods, sim_T, degradation_factor = load_the_data()
