# import packages
using JuMP
using Gurobi
using Printf
using Plots

# Load data from 02435_two_stage_problem_data.jl function load_the_data()
include("V2_Assignment_A_codes/V2_02435_two_stage_problem_data.jl")
number_of_warehouses, W, cost_miss, cost_tr, warehouse_capacities, transport_capacities, initial_stock, number_of_simulation_periods, sim_T, degradation_factor = load_the_data()

# take only one experiment of demand and price
demand = 4*ones(number_of_warehouses, number_of_simulation_periods)

include("V2_Assignment_A_codes/V2_simulation_experiments.jl")
number_of_experiments, Expers, price_trajectory = simulation_experiments_creation(number_of_warehouses, W, number_of_simulation_periods)
# price t=1
p_wt = price_trajectory[1,:,1]

include("V2_Assignment_A_codes/V2_price_process.jl")
S = 1000

model_1d = Model(Gurobi.Optimizer)

function Scenario_generation(p_wt, S, number_of_warehouses, num_reduced, plots=true)
    # create S price scenarios using the sample_next function
    p_wt_scenarios = zeros(S, number_of_warehouses)
    for s in 1:S
        for w in W
            p_wt_scenarios[s,w] = sample_next(p_wt[w])
        end
    end

    # Calculate distance matrix (euclidean distance) for scenario reduction
    D = zeros(Float64, S,S)
    for i = 1:S
        for j = 1:S
            D[i,j] = sqrt(sum((p_wt_scenarios[i,l]-p_wt_scenarios[j,l])*(p_wt_scenarios[i,l]-p_wt_scenarios[j,l])  for l = 1:number_of_warehouses))
        end
    end


    #Initialize equiprobable probabilities
    probabilities = repeat([1.0/S], 1, S)[1,:]
    # Use fast forward selection and apply it to get the reduced scenarios and updated probabilities
    include("fast-forward-selection.jl")
    result = FastForwardSelection(D, probabilities, num_reduced)
    #Resulting probabilities
    reduced_probabilities = result[1]
    print("Reduced probabilities: ", reduced_probabilities)
    #Selected scenario indices
    reduced_scenario_indices = result[2]
    print("Reduced scenario indices: ", reduced_scenario_indices)
    #reduced_price_scenarios = p_wt_scenarios[reduced_scenario_indices,:]
    # Filter p_wt_scenarios based on reduced_scenario_indices
    reduced_price_scenarios = zeros(Float64, num_reduced, number_of_warehouses)

    for i = 1:num_reduced
        reduced_price_scenarios[i,:] = p_wt_scenarios[reduced_scenario_indices[i],:]
    end

    if plots
        histogram(p_wt_scenarios[:,:], xlabel="Price", ylabel="Frequency", title="Historgram of price scenarios", alpha=0.8)
        savefig("Assignment_1/Plots/histogram_price_scenarios.pdf")
        #transpose p_wt_scenarios for plotting
        p_wt_scenarios_t = transpose(p_wt_scenarios)
        plot(p_wt_scenarios_t, xlabel="Warehouse", ylabel="Price", title="Reduced price scenarios", color=:lightgray, legend=false, alpha=0.8, xticks=(1:number_of_warehouses, W))
        #Plot reduced scenarios
        reduced_data = zeros(Float64, number_of_warehouses, num_reduced)
        for i = 1:num_reduced
            reduced_data[:,i] = p_wt_scenarios_t[:,result[2][i]]
        end
        plot!(reduced_data, legend=false, color=:auto)
        savefig("Assignment_1/Plots/reduced_scenarios.pdf")
    end 
    
    return reduced_price_scenarios, reduced_probabilities
end    

function Make_Stochastic_here_and_now_decision(p_wt, S, num_reduced)

    price_scenarios, probabilities = Scenario_generation(p_wt, S, number_of_warehouses, 10, true)
    N = collect(1:num_reduced)
    # Define the model
    model_1d = Model(Gurobi.Optimizer)


    @variable(model_1d, 0 <= x_wt[w in W]) # Coffee ordered at warehouse w in period t
    @variable(model_1d, 0 <= z_wt[w in W]) # Coffee stored at warehouse w in period t
    @variable(model_1d, 0 <= y_send_wqt[w in W, q in W]) # Coffee sent from warehouse w to warehouse q in period t
    @variable(model_1d, 0 <= y_receive_wqt[w in W, q in W]) # Coffee received at warehouse w from warehouse q in period t
    @variable(model_1d, 0 <= m_wt[w in W]) # Missing demand at warehouse w in period t
    @variable(model_1d, 0 <= x_wt_scen[w in W, s in N]) # Coffee ordered at warehouse w in period t in scenario s
    @variable(model_1d, 0 <= z_wt_scen[w in W, s in N]) # Coffee stored at warehouse w in period t in scenario s
    @variable(model_1d, 0 <= y_send_wqt_scen[w in W, q in W, s in N]) # Coffee sent from warehouse w to warehouse q in period t
    @variable(model_1d, 0 <= y_receive_wqt_scen[w in W, q in W, s in N]) # Coffee received at warehouse w from warehouse q in period t
    @variable(model_1d, 0 <= m_wt_scen[w in W, s in N]) # Missing demand at warehouse w in period t
    
    
    # Define the objective function
    @objective(model_1d, Min, sum(p_wt[w]*x_wt[w] for w in W) 
                            + sum(cost_miss[w]*m_wt[w] for w in W) 
                            + sum(cost_tr[w,q]*y_send_wqt[w,q] for w in W, q in W)
                            + sum(probabilities[s] * price_scenarios[s, w]*x_wt_scen[w,s] for w in W, s in N)
                            + sum(probabilities[s] * cost_miss[w]*m_wt_scen[w,s] for w in W, s in N)
                            + sum(probabilities[s] * cost_tr[w,q]*y_send_wqt_scen[w,q,s] for w in W, q in W, s in N)
                            )

    # Define the constraints
    # Demand hour 1
    @constraint(model_1d, demand_t0[w in W],
    x_wt[w] - z_wt[w] + initial_stock[w]
    + sum(y_receive_wqt[w,q] for q in W)
    - sum(y_send_wqt[w,q] for q in W)
    + m_wt[w] == demand[w,1])

    # Demand hour t
    @constraint(model_1d, demand_t[w in W, s in N],
    x_wt_scen[w,s] - z_wt_scen[w,s] + z_wt[w]
    + sum(y_receive_wqt_scen[w,q,s] for q in W) 
    - sum(y_send_wqt_scen[w,q,s] for q in W)
    + m_wt_scen[w,s] == demand[w,2])

    # Storage capacity
    @constraint(model_1d, storage_capacity[w in W],
    z_wt[w] <= warehouse_capacities[w])

    # Transportation capacity
    @constraint(model_1d, transportation_capacity[w in W, q in W],
    y_send_wqt[w,q] <= transport_capacities[w,q])

    # Constraint saying the amount sent from w to q is the same as received at q from w
    @constraint(model_1d, send_receive[w in W, q in W], 
        y_send_wqt[w,q] == y_receive_wqt[q,w])


    #Storage one day before transfer, initial stock in t=1
    @constraint(model_1d, storage_t1[w in W, q in W],
    y_send_wqt[w,q] <= initial_stock[w])

    # Storage capacity
    @constraint(model_1d, storage_capacity_s[w in W, s in N],
    z_wt_scen[w,s] <= warehouse_capacities[w])

    # Transportation capacity
    @constraint(model_1d, transportation_capacity_s[w in W, q in W, s in N],
    y_send_wqt_scen[w,q,s] <= transport_capacities[w,q])

    # Constraint saying the amount sent from w to q is the same as received at q from w
    @constraint(model_1d, send_receive_s[w in W, q in W, s in N],
    y_send_wqt_scen[w,q,s] == y_receive_wqt_scen[q,w,s])

    @constraint(model_1d, storage_before_transfer_s[w in W, q in W, s in N],
    y_send_wqt_scen[w,q,s] <= z_wt[w])
                
    # Solve the model
    optimize!(model_1d)

    # Return the optimal solution
    if termination_status(model_1d) == MOI.OPTIMAL
        println("Optimal solution found")
        # Print results for Day 1
        println("Day 1:")
        for w in W
            println("Warehouse ", w)
            @printf("Price: %0.3f\n", value(p_wt[w]))
            @printf("Order: %0.3f\n", value(x_wt[w]))
            @printf("Storage: %0.3f\n", value(z_wt[w]))
            for q in W
                if q != w
                    @printf("Send to warehouse %i: %0.3f\n", q, value(y_send_wqt[w,q]))
                    @printf("Receive from warehouse %i: %0.3f\n", q, value(y_receive_wqt[w,q]))
                end
            end
        end
             
        println("-------")

        # Print results for Day 2
        println("Day 2:")
        for s in N
            println("Scenario" , s)
            for w in W
                println("Warehouse ", w)
                @printf("Expected Price: %0.3f\n", value(price_scenarios[s,w]))
                @printf("Order: %0.3f\n", value(x_wt_scen[w,s]))
                @printf("Storage: %0.3f\n", value(z_wt_scen[w,s]))
            for q in W
                if q != w
                    @printf("Send to warehouse %i: %0.3f\n", q, value(y_send_wqt_scen[w,q,s]))
                    @printf("Receive from warehouse %i: %0.3f\n", q, value(y_receive_wqt_scen[w,q,s]))
                end
            end 
            end
        end
        
        system_cost = objective_value(model_1d)
        println("Total system cost: ", system_cost)
    else
        println("No solution found")
    end 
end

Make_Stochastic_here_and_now_decision(p_wt, S, 10)