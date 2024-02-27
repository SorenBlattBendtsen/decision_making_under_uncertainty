#1e-Experiments of all Methods
# Important comment:
# When you run this, you might get an error regarding the Fast Forward Selection file.
# If that is the case, you should open the file fast-forward-selection.jl and run each function one by one.
# And then run this file again.


# import packages
using JuMP
using Gurobi
using Printf
using Plots
using Random

# random seed for reproducability
Random.seed!(1234)

# Include files for input data
include("V2_Assignment_A_codes/V2_02435_two_stage_problem_data.jl")
number_of_warehouses, W, cost_miss, cost_tr, warehouse_capacities, transport_capacities, initial_stock, number_of_simulation_periods, sim_T, degradation_factor = load_the_data()

include("V2_Assignment_A_codes/V2_simulation_experiments.jl")
number_of_experiments, Expers, price_trajectory = simulation_experiments_creation(number_of_warehouses, W, number_of_simulation_periods)

include("V2_Assignment_A_codes/V2_price_process.jl")

# Constant demand for all warehouses and all periods
demand = 4*ones(number_of_warehouses, number_of_simulation_periods)

# intial price, 100 experiments
price = zeros(number_of_experiments, number_of_warehouses, number_of_simulation_periods)
for i in 1:number_of_experiments
    for j in 1:number_of_warehouses
        for t in 1:number_of_simulation_periods
            price[i,j,t] = rand(Uniform(0,10))
        end
    end
end

# 1000 scenarios, 10 reduced scenarios
S = 1000
num_reduced_values = [5, 20, 50]

# Include 1b, 1c and 1d:
include("1c_optimality_in_hindsight.jl")
include("1b Expected Value.jl")
include("1d_two_stage.jl")
# Do for loop,
oih_system_costs = zeros(number_of_experiments)
ev_system_costs = zeros(number_of_experiments)
ts_system_costs = zeros(number_of_experiments, length(num_reduced_values))
for i in 1:number_of_experiments
    print("Running experiment $i")
    print("\n")
    print("Running Optimality in Hindsight")
    result_oih = Calculate_OiH_solution(price[i,:,:])
    oih_system_costs[i] = result_oih[1]
    print("Running Expected Value")
    result_ev = make_EV_here_and_now_decision(price[i,:,1])
    ev_system_costs[i] = result_ev[1]
    print("Running Two-Stage")
    j = 1
    for num_reduced in num_reduced_values
        result_ts = Make_Stochastic_here_and_now_decision(price[i,:,1], S, num_reduced)
        ts_system_costs[i, j] = result_ts[1]
        j += 1
    end 
end

# results
function plot_results(oih_system_costs, ev_system_costs, ts_system_costs)
    # Optimality in Hindsight
    histogram(oih_system_costs, bins=20, label="OiH", xlabel="System Cost", ylabel="Frequency", title="Optimality in Hindsight")
    # Average line
    avg_oih = round(mean(oih_system_costs), digits=2)
    vline!([avg_oih], label="Mean", color="red")
    annotate!([(avg_oih, 13, text("Mean = $avg_oih", 10, :left))])
    savefig("Assignment_1/Plots/1e_OiH.png")

    # Expected Value
    histogram(ev_system_costs, bins=20, label="EV", xlabel="System Cost", ylabel="Frequency", title="Expected Value")
    # Average line
    avg_ev = round(mean(ev_system_costs), digits=2)
    vline!([avg_ev], label="Mean", color="red")
    annotate!([(avg_ev, 10, text("Mean = $avg_ev", 10, :left))])
    savefig("Assignment_1/Plots/1e_EV.png")

    # Two-Stage
    histogram(ts_system_costs[:,1], bins=20, label="S=5", xlabel="System Cost", ylabel="Frequency", title="Two-Stage 5 Scenarios")
    # Average line
    avg_ts_1 = round(mean(ts_system_costs[:,1]), digits=2)
    vline!([avg_ts_1], label="Mean", color="red")
    annotate!([(avg_ts_1, 9, text("Mean = $avg_ts_1", 10, :right))])
    savefig("Assignment_1/Plots/1e_TS_5.png")

    histogram(ts_system_costs[:,2], bins=20, label="S=20", xlabel="System Cost", ylabel="Frequency", title="Two-Stage 20 Scenarios")
    # Average line
    avg_ts_2 = round(mean(ts_system_costs[:,2]), digits=2)
    vline!([avg_ts_2], label="Mean", color="red")
    annotate!([(avg_ts_2, 9, text("Mean = $avg_ts_2", 10, :right))])
    savefig("Assignment_1/Plots/1e_TS_20.png")

    histogram(ts_system_costs[:,3], bins=20, label="S=50", xlabel="System Cost", ylabel="Frequency", title="Two-Stage 50 Scenarios")
    # Average line
    avg_ts_3 = round(mean(ts_system_costs[:,3]), digits=2)
    vline!([avg_ts_3], label="Mean", color="red")
    annotate!([(avg_ts_3, 8, text("Mean = $avg_ts_3", 10, :left))])
    savefig("Assignment_1/Plots/1e_TS_50.png")
end 
plot_results(oih_system_costs, ev_system_costs, ts_system_costs)