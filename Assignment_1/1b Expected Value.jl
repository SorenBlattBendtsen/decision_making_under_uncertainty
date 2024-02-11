#Import packages
using JuMP
using Gurobi
using Printf
using Distributions

# Load data from 02435_two_stage_problem_data.jl function load_the_data()
include("02435_Assignment_A_codes/02435_two_stage_problem_data.jl")
number_of_warehouses, W, cost_miss, cost_tr, warehouse_capacities, transport_capacities, initial_stock, number_of_simulation_periods, sim_T, degradation_factor = load_the_data()

# Generate random data for day 1
price = rand(Uniform(0,10),10,10)

# Include the file containing the price process function
include("02435_Assignment_A_codes/price_process.jl")

function calculate_expected_prices(price)
         
    #     return expected_price
    num_samples = 1000
    expected_price = zeros(size(price))  # Initialize expected prices matrix

    for w in 1:size(price, 1)  # Loop through each warehouse
        for t in 1:size(price, 3)  # Loop through each time period
            initial_prices = price[w, 1, t]  # Extract initial prices as a scalar or vector
            expected_price_tmp = zeros(size(price, 2))  # Temporary variable for storing expected prices
            
            for i in 1:num_samples
                price_sample = sample_next(initial_prices)  # Use initial prices as input to sample_next
                expected_price_tmp .+= price_sample  # Aggregate samples
            end
            
            expected_price_tmp ./= num_samples  # Average over samples
            expected_price[w, :, t] = expected_price_tmp  # Assign expected prices back to the main matrix
        end
    end
    
    return expected_price
    
end

    

function make_EV_here_and_now_decision(price)
    expected_price = calculate_expected_prices(price)
    
    # Define the model
    model_ev = Model(Gurobi.Optimizer)

    # Variables
    @variable(model_ev, 0 <= x_wt[w in W, t in 1:2]) # Coffee ordered at warehouse w in period t
    @variable(model_ev, 0 <= z_wt[w in W, t in 1:2]) # Coffee stored at warehouse w in period t
    @variable(model_ev, 0 <= y_send_wqt[w in W, q in W, t in 1:2]) # Coffee sent from warehouse w to warehouse q in period t
    @variable(model_ev, 0 <= y_receive_wqt[w in W, q in W, t in 1:2]) # Coffee received at warehouse w from warehouse q in period t
    @variable(model_ev, 0 <= m_wt[w in W, t in 1:2]) # Missing demand at warehouse w in period t

    demand_fixed = 4 
    initial_stock = 2

    # Define the objective function
    @objective(model_ev, Min, sum(expected_price[w,t] * x_wt[w,t] for w in W, t in 1:2) 
                            + sum(cost_miss[w] * m_wt[w,t] for w in W, t in 1:2) 
                            + sum(cost_tr[w,q] * y_send_wqt[w,q,t] for w in W, q in W, t in 1:2))

    # Define the constraints
    # initial_stock
    @constraint(model_ev, ini_stock[w in W], z_wt[w,1] == initial_stock)

    # demand hour 1
    @constraint(model_ev, demand_t0[w in W],
                x_wt[w,1] + z_wt[w,1]
                + sum(y_receive_wqt[w,q,1] for q in W)
                - sum(y_send_wqt[w,q,1] for q in W)
                + m_wt[w,1] == demand_fixed)
    
    # demand hour 2
    @constraint(model_ev, demand_t[w in W],
                x_wt[w,2] + z_wt[w,2] - z_wt[w,1]
                + sum(y_receive_wqt[w,q,2] for q in W) 
                - sum(y_send_wqt[w,q,2] for q in W)
                + m_wt[w,2] == demand_fixed)
    
    # Storage capacity
    @constraint(model_ev, storage_capacity[w in W],
                z_wt[w,1] <= warehouse_capacities[w])

    # Transportation capacity
    @constraint(model_ev, transportation_capacity[w in W, q in W],
                y_send_wqt[w,q,1] <= transport_capacities[w,q])

    # Constraint to prevent transfers to the same warehouse
    @constraint(model_ev, no_self_transfer[w in W], y_send_wqt[w,w,1] == 0)

    # constraint saying the amount sent from w to q is the same as received at q from w
    @constraint(model_ev, send_receive[w in W, q in W], 
                y_send_wqt[w,q,1] == y_receive_wqt[q,w,1])

    # storage duration before transfer
    @constraint(model_ev, storage_before_transfer[w in W, q in W],
                y_send_wqt[w,q,1] <= z_wt[w,1])
    @constraint(model_ev, storage_before_transfer_t1[w in W, q in W],
    y_send_wqt[w,q,1] <= initial_stock)
    
    # Solve the model
    optimize!(model_ev)

    # Return the optimal solution
    if termination_status(model_ev) == MOI.OPTIMAL
        println("Optimal solution found")
        for w in W
            println("Warehouse ", w)
            @printf("Day 1 Order: %0.3f\n", value(x_wt[w,1]))
            @printf("Day 1 Storage: %0.3f\n", value(z_wt[w,1]))
            @printf("Day 2 Order: %0.3f\n", value(x_wt[w,2]))
            @printf("Day 2 Storage: %0.3f\n", value(z_wt[w,2]))
            for q in W
                if q != w
                    @printf("Day 1 Send to warehouse %i: %0.3f\n", q, value(y_send_wqt[w,q,1]))
                    @printf("Day 1 Receive from warehouse %i: %0.3f\n", q, value(y_receive_wqt[w,q,1]))
                    @printf("Day 2 Send to warehouse %i: %0.3f\n", q, value(y_send_wqt[w,q,2]))
                    @printf("Day 2 Receive from warehouse %i: %0.3f\n", q, value(y_receive_wqt[w,q,2]))
                end
            end
        end
        system_cost = objective_value(model_ev)
        println("Total system cost: ", system_cost)
    else
        println("No solution found")
    end 
end

make_EV_here_and_now_decision(price)
