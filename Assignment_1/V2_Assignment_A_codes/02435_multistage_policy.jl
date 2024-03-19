# Combined Multi-Stage Stochastic Problem Policy in Julia

# Initialization and Data Loading
using JuMP, Gurobi
using Distributions
include("V2_02435_multistage_problem_data.jl")
include("V2_price_process.jl")

# Load data and define constants
number_of_warehouses, W, cost_miss, cost_tr, warehouse_capacities, transport_capacities, initial_stock, number_of_simulation_periods, sim_T, demand_trajectory = load_the_data()
const NUM_SCENARIOS = 1000
const NUM_WAREHOUSES = number_of_warehouses
const NUM_STAGES = number_of_simulation_periods
const DISCRETE_PRICE_LEVELS = [0.0, 2.5, 5.0, 7.5, 10.0]
const TARGET_NUM_SCENARIOS = 50

# Scenario Generation and Discretization
function generate_scenarios(num_scenarios, num_stages, num_warehouses)
    price_scenarios = zeros(num_warehouses, num_stages, num_scenarios)
    for s in 1:num_scenarios
        for t in 1:num_stages
            for w in 1:num_warehouses
                price_scenarios[w, t, s] = t == 1 ? initial_stock[w] : sample_next(price_scenarios[w, t-1, s])
            end
        end
    end
    return price_scenarios
end
initial_price_scenarios = generate_scenarios(NUM_SCENARIOS, NUM_STAGES, NUM_WAREHOUSES)

# Function to discretize all prices in the scenarios
function discretize_all_prices(prices::Array{Float64,3}, discrete_levels::Array{Float64,1})
    num_warehouses, num_stages, num_scenarios = size(prices)
    discretized_prices = copy(prices)
    
    for w in 1:num_warehouses
        for t in 1:num_stages
            for s in 1:num_scenarios
                # Find the closest discrete level for each price
                _, nearest_index = findmin(abs.(discrete_levels .- prices[w, t, s]))
                discretized_prices[w, t, s] = discrete_levels[nearest_index]
            end
        end
    end
    
    return discretized_prices
end


discretized_price_scenarios = discretize_all_prices(initial_price_scenarios, DISCRETE_PRICE_LEVELS)

# Including Fast-Forward Selection for Scenario Reduction
include("fast-forward-selection.jl")  # Make sure this path is correct

# Function to calculate distance matrix for scenario reduction
function calculate_distance_matrix(prices::Array{Float64,3})
    _, _, num_scenarios = size(prices)
    D = zeros(num_scenarios, num_scenarios)
    for i in 1:num_scenarios
        for j in 1:num_scenarios
            if i != j
                # Calculating the Euclidean distance between the price trajectories of two scenarios
                distance = 0.0
                for t in 1:size(prices, 2)  # Iterate over all time stages
                    for w in 1:size(prices, 1)  # Iterate over all warehouses
                        distance += (prices[w, t, i] - prices[w, t, j])^2
                    end
                end
                D[i, j] = sqrt(distance)
            end
        end
    end
    return D
end


distance_matrix = calculate_distance_matrix(discretized_price_scenarios)
initial_probabilities = fill(1.0 / NUM_SCENARIOS, NUM_SCENARIOS)
new_probabilities, selected_indices = FastForwardSelection(distance_matrix, initial_probabilities, TARGET_NUM_SCENARIOS)
reduced_price_scenarios = discretized_price_scenarios[:, :, selected_indices]

# Optimization Problem Setup
model = Model(Gurobi.Optimizer)
@variable(model, x[w=1:NUM_WAREHOUSES, t=1:NUM_STAGES, s=1:TARGET_NUM_SCENARIOS] >= 0)
@variable(model, z[w=1:NUM_WAREHOUSES, t=1:NUM_STAGES, s=1:TARGET_NUM_SCENARIOS] >= 0)
@variable(model, y_send[w=1:NUM_WAREHOUSES, q=1:NUM_WAREHOUSES, t=1:NUM_STAGES, s=1:TARGET_NUM_SCENARIOS] >= 0)
@variable(model, m[w=1:NUM_WAREHOUSES, t=1:NUM_STAGES, s=1:TARGET_NUM_SCENARIOS] >= 0)

# Constraints Definition

# Inventory Balance Constraints
@constraint(model, balance[w=1:NUM_WAREHOUSES, t=2:NUM_STAGES, s=1:TARGET_NUM_SCENARIOS], 
    z[w, t, s] == z[w, t-1, s] + x[w, t, s] - sum(y_send[w, q, t, s] for q=1:NUM_WAREHOUSES) + 
    sum(y_send[q, w, t, s] for q=1:NUM_WAREHOUSES) - m[w, t, s])

# Demand Fulfillment Constraints
@constraint(model, demand_fulfillment[w=1:NUM_WAREHOUSES, t=1:NUM_STAGES, s=1:TARGET_NUM_SCENARIOS],
    x[w, t, s] + z[w, t, s] + sum(y_send[q, w, t, s] for q=1:NUM_WAREHOUSES) - 
    sum(y_send[w, q, t, s] for q=1:NUM_WAREHOUSES) - m[w, t, s] == demand_trajectory[w, t])

# Transportation Capacity Constraints
@constraint(model, transportation[w=1:NUM_WAREHOUSES, q=1:NUM_WAREHOUSES, t=1:NUM_STAGES, s=1:TARGET_NUM_SCENARIOS],
    y_send[w, q, t, s] <= transport_capacities[w, q])

# Storage Capacity Constraints
@constraint(model, storage[w=1:NUM_WAREHOUSES, t=1:NUM_STAGES, s=1:TARGET_NUM_SCENARIOS],
    z[w, t, s] <= warehouse_capacities[w])

# Non-Negativity and Initial Conditions
@constraint(model, initial_inventory[w=1:NUM_WAREHOUSES, s=1:TARGET_NUM_SCENARIOS], 
    z[w, 1, s] == initial_stock[w])  # Setting initial inventory levels

# Ensuring decisions for coffee sent between the same warehouse are zero
@constraint(model, no_self_send[w=1:NUM_WAREHOUSES, t=1:NUM_STAGES, s=1:TARGET_NUM_SCENARIOS],
    y_send[w, w, t, s] == 0)

# Ensuring initial orders and unmet demands at the start are zero (if necessary)
@constraint(model, initial_orders_and_shortages[s=1:TARGET_NUM_SCENARIOS], 
    sum(x[w, 1, s] for w=1:NUM_WAREHOUSES) + sum(m[w, 1, s] for w=1:NUM_WAREHOUSES) == 0)


# Objective Function and Non-Anticipativity Constraints

# Objective Function: Minimize total expected costs
@objective(model, Min, 
    sum(new_probabilities[s] * (cost_miss[w] * m[w, t, s] + 
    sum(cost_tr[w, q] * y_send[w, q, t, s] for q = 1:NUM_WAREHOUSES) +
    sum(cost_tr[q, w] * y_send[q, w, t, s] for q = 1:NUM_WAREHOUSES) +
    x[w, t, s] * reduced_price_scenarios[w, t, s]) for w = 1:NUM_WAREHOUSES, t = 1:NUM_STAGES, s = 1:TARGET_NUM_SCENARIOS))

# Dummy function for generating non-anticipativity sets based on reduced scenarios
function generate_non_anticipativity_sets(num_stages, num_scenarios, reduced_scenarios)
    # Initialize a dictionary to store non-anticipativity sets for each stage
    non_anticipativity_sets = Dict{Int, Array{Set{Int}, 1}}()
    
    # For each stage, determine which scenarios share the same information up to that stage
    for t in 1:num_stages
        sets_at_t = Array{Set{Int}, 1}()
        already_grouped = Set{Int}()
        
        for s in 1:num_scenarios
            if s in already_grouped
                continue
            end
            # Start a new group with the current scenario
            current_set = Set{Int}([s])
            for other_s in s+1:num_scenarios
                # Compare reduced_scenarios[s] and reduced_scenarios[other_s] up to stage t
                # If they match, add other_s to the current_set
                if all(reduced_scenarios[:, 1:t, s] .== reduced_scenarios[:, 1:t, other_s])
                    push!(current_set, other_s)
                    push!(already_grouped, other_s)
                end
            end
            push!(sets_at_t, current_set)
        end
        non_anticipativity_sets[t] = sets_at_t
    end
    return non_anticipativity_sets
end

# Generate non-anticipativity sets based on the reduced scenarios
non_anticipativity_sets = generate_non_anticipativity_sets(NUM_STAGES, TARGET_NUM_SCENARIOS, reduced_price_scenarios)


# Non-Anticipativity Constraints
# Ensuring decisions made based on the same information up to time t are the same across all scenarios that share this information
for t in 1:NUM_STAGES
    for na_set in non_anticipativity_sets[t]  # Assuming non_anticipativity_sets have been prepared as discussed
        for w in 1:NUM_WAREHOUSES
            # Get a reference scenario from the set
            ref_scenario = first(collect(na_set))
            # Apply non-anticipativity by forcing decision variables to be equal for all scenarios in the set
            for other_scenario in setdiff(collect(na_set), [ref_scenario])
                @constraint(model, x[w, t, other_scenario] == x[w, t, ref_scenario])
                @constraint(model, m[w, t, other_scenario] == m[w, t, ref_scenario])
                for q in 1:NUM_WAREHOUSES
                    @constraint(model, y_send[w, q, t, other_scenario] == y_send[w, q, t, ref_scenario])
                end
            end
        end
    end
end


# Optimization Execution and Solution Extraction
optimize!(model)
if termination_status(model) == MOI.OPTIMAL
    println("Optimal solution found!")
    # (Insert solution extraction and analysis code, as outlined in Part 5.)
else
    println("Optimal solution not found. Status: ", termination_status(model))
end

# Solution Analysis and Finalization

# Check if an optimal solution has been found and extract it
if termination_status(model) == MOI.OPTIMAL
    # Extracting decision variables
    optimal_orders = getvalue.(x)  # Coffee ordered in optimal solution
    optimal_inventory = getvalue.(z)  # Inventory levels in optimal solution
    optimal_send = getvalue.(y_send)  # Coffee sent between warehouses in optimal solution
    optimal_shortages = getvalue.(m)  # Unmet demand in optimal solution

    # Printing out some of the solutions for inspection
    println("Optimal Orders: ", optimal_orders)
    println("Optimal Inventory Levels: ", optimal_inventory)
    println("Optimal Coffee Sent Between Warehouses: ", optimal_send)
    println("Optimal Shortages: ", optimal_shortages)

    # Post-optimization analysis or additional steps could be performed here
    # For example, calculating total costs, preparing reports, etc.
else
    println("Optimal solution not found. Status: ", termination_status(model))
end

# If you need to use the solution outside of this script, you can save it to files or prepare it for visualization.
# For example:
# CSV.write("optimal_orders.csv", DataFrame(optimal_orders))
# Make sure to include necessary packages and adjust paths/names as needed.

# Final cleanup or post-processing can be done below
println("Policy execution complete. Check the outputs and validate against constraints and expectations.")

# Optional: Comparing the policy's performance against the dummy policy, if required
# This could involve setting up a comparison based on total costs, service levels, or other KPIs defined in your problem.

