# import packages
using JuMP
using Gurobi
using Printf

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

function Make_Stochastic_here_and_now_decision(price, S, N)
    model_1d = Model(Gurobi.Optimizer)

    p_wt_scenarios = zeros(S, number_of_warehouses)
    for s in 1:S
        for w in W
            p_wt_scenarios[s,w] = sample_next(price[w])
        end
    end
    # take N scenarios
    #p_wt = p_wt_scenarios[1,:]
    p_wt_scenarios = p_wt_scenarios[1:N,:]

    # Variables, 1st stage and 2nd stage
    @variable(model_1d, 0 <= x_wt[w in W, t in sim_T]) # Coffee ordered at warehouse w in period t
    @variable(model_1d, 0 <= z_wt[w in W, t in sim_T]) # Coffee stored at warehouse w in period t
    @variable(model_1d, 0 <= y_send_wqt[w in W, q in W, t in sim_T]) # Coffee sent from warehouse w to warehouse q in period t
    @variable(model_1d, 0 <= y_receive_wqt[w in W, q in W, t in sim_T]) # Coffee received at warehouse w from warehouse q in period t
    @variable(model_1d, 0 <= m_wt[w in W, t in sim_T]) # Missing demand at warehouse w in period t
    @variable(model_1d, 0 <= x_wt_scenarios[w in W, t in sim_T, s in N]) # Coffee ordered at warehouse w in period t in scenario s
    @variable(model_1d, 0 <= y_send_wqt_scen[w in W, q in W, t in sim_T, s in N]) # Coffee sent from warehouse w to warehouse q in period t
    @variable(model_1d, 0 <= y_receive_wqt_scen[w in W, q in W, t in sim_T, s in N]) # Coffee received at warehouse w from warehouse q in period t
    @variable(model_1d, 0 <= m_wt_scen[w in W, t in sim_T, s in N]) # Missing demand at warehouse w in period t
    # Define the objective function
    @objective(model_1d, Min, sum(p_wt[w]*x_wt[w,t] for w in W, t in sim_T) 
                            + sum(cost_miss[w]*m_wt[w,t] for w in W, t in sim_T) 
                            + sum(cost_tr[w,q]*y_send_wqt[w,q,t] for w in W, q in W, t in sim_T)
                            + sum(1/N * p_wt_scenarios[s, w]*x_wt_scenarios[w,t+1,s] for w in W, t in sim_T[1:end-1], s in N)
                            + sum(1/N * cost_miss[w]*m_wt_scen[w,t+1,s] for w in W, t in sim_T[1:end-1], s in N)
                            + sum(1/N * cost_tr[w,q]*y_send_wqt_scen[w,q,t+1,s] for w in W, q in W, t in sim_T[1:end-1], s in N)
                            )

    # Define the constraints
    @constraint(model_1d, [w in W, t in sim_T[1:end-1]], sum(x_wt_scenarios[w,t+1,s] for s in N) == x_wt[w,t+1]) # Scenario Definition
    # @constraint(model_1d, [w in W, t in sim_T], sum(p_wt_scenarios[s,w] * x_wt_scenarios[w,t,s] for s in N) == 1/N * sum(p_wt_scenarios[s,w] * x_wt[w,t] for s in N)) # Scenario probabilities

     # demand hour 1
     @constraint(model_1d, demand_fulfillment_t0[w in W],
     x_wt[w,1] - z_wt[w,1] + initial_stock[w]
     + sum(y_receive_wqt[w,q,1] for q in W)
     - sum(y_send_wqt[w,q,1] for q in W)
     + m_wt[w,1] == demand[w,1])

    #demand hour t
    @constraint(model_1d, demand_fulfillment[w in W, t in sim_T[2:end]],
        x_wt[w,t] - z_wt[w,t] + z_wt[w,t-1]
        + sum(y_receive_wqt[q,w,t] for q in W) 
        - sum(y_send_wqt[w,q,t] for q in W)
        + m_wt[w,t] == demand[w,t])

    #demand scenarios
    @constraint(model_1d, demand_fulfillmen_scenarios[w in W, t in sim_T[1:end-1], s in N],
                x_wt_scenarios[w,t+1,s] - z_wt[w,t+1] + z_wt[w,t]
                + sum(y_receive_wqt_scen[q,w,t+1,s] for q in W) 
                - sum(y_send_wqt_scen[w,q,t+1,s] for q in W)
                + m_wt_scen[w,t+1,s] == demand[w,t+1])

    # Storage capacity
    @constraint(model_1d, storage_capacity[w in W, t in sim_T],
        z_wt[w,t] <= warehouse_capacities[w])

    # Transportation capacity
    @constraint(model_1d, transportation_capacity[w in W, q in W, t in sim_T],
        y_send_wqt[w,q,t] <= transport_capacities[w,q])

    # constraint saying the amount sent from w to q is the same as received at q from w
    @constraint(model_1d, send_receive[w in W, q in W, t in sim_T], 
        y_send_wqt[w,q,t] == y_receive_wqt[q,w,t])


    #Storage one day before transfer, initial stock in t=1
    @constraint(model_1d, storage_t1[w in W, q in W],
        y_send_wqt[w,q,1] <= initial_stock[w])

    @constraint(model_1d, storage_before_transfer[w in W, q in W, t in sim_T[2:end]],
        y_send_wqt[w,q,t] <= z_wt[w,t-1])
        
    # Solve the model
    optimize!(model_1d)

    # Print the solution
    if termination_status(model_1d) == MOI.OPTIMAL
        println("Optimal solution found")
        for t in sim_T
            println("Period ", t)
            for w in W
                println("Warehouse ", w)
                @printf("Order: %0.3f\n", value(x_wt[w,t]))
                @printf("Storage: %0.3f\n", value(z_wt[w,t]))
                for q in W
                    if q != w
                        @printf("Send to warehouse %i: %0.3f\n", q, value(y_send_wqt[w,q,t]))
                        @printf("Receive from warehouse %i: %0.3f\n", q, value(y_receive_wqt[w,q,t]))
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

