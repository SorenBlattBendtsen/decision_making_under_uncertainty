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

function Make_Stochastic_here_and_now_decision(price, S, N)
    model_1d = Model(Gurobi.Optimizer)

    p_wt_scenarios = zeros(number_of_warehouses, number_of_simulation_periods, S)
    for s in 1:S
        for w in W
            for t in sim_T
                p_wt_scenarios[w,t,s] = sample_next(price[w])
            end
        end
    end
    # take N scenarios
    #p_wt = p_wt_scenarios[1,:]
    p_wt_scenarios = p_wt_scenarios[:,:,1:N]
    for s in 1:N
        p_wt_scenarios[:,1,s] = p_wt
    end


    # Variables, 1st stage and 2nd stage
    @variable(model_1d, 0 <= z_wt[w in W, t in sim_T, s in N]) # Coffee stored at warehouse w in period t
    @variable(model_1d, 0 <= x_wt[w in W, t in sim_T, s in N]) # Coffee ordered at warehouse w in period t in scenario s
    @variable(model_1d, 0 <= y_send_wqt[w in W, q in W, t in sim_T, s in N]) # Coffee sent from warehouse w to warehouse q in period t
    @variable(model_1d, 0 <= y_receive_wqt[w in W, q in W, t in sim_T, s in N]) # Coffee received at warehouse w from warehouse q in period t
    @variable(model_1d, 0 <= m_wt[w in W, t in sim_T, s in N]) # Missing demand at warehouse w in period t
    # Define the objective function
    @objective(model_1d, Min, (sum(p_wt_scenarios[w,t,s]*x_wt[w,t,s] for w in W, t in sim_T, s in N)
                            + sum(cost_miss[w]*m_wt[w,t,s] for w in W, t in sim_T, s in N)
                            + sum(cost_tr[w,q]*y_send_wqt[w,q,t,s] for w in W, q in W, t in sim_T, s in N))
                            )

    # Define the constraints

     # demand hour 1
     @constraint(model_1d, demand_fulfillment_t0[w in W, s in N],
     x_wt[w,1,s] - z_wt[w,1,s] + initial_stock[w]
     + sum(y_receive_wqt[w,q,1,s] for q in W)
     - sum(y_send_wqt[w,q,1,s] for q in W)
     + m_wt[w,1,s] == demand[w,1])

    #demand hour t
    @constraint(model_1d, demand_fulfillment[w in W, t in sim_T[2:end], s in N],
        x_wt[w,t,s] - z_wt[w,t,s] + z_wt[w,t-1,s]
        + sum(y_receive_wqt[q,w,t,s] for q in W) 
        - sum(y_send_wqt[w,q,t,s] for q in W)
        + m_wt[w,t,s] == demand[w,t])

    # Storage capacity
    @constraint(model_1d, storage_capacity[w in W, t in sim_T, s in N],
        z_wt[w,t,s] <= warehouse_capacities[w])

    # Transportation capacity
    @constraint(model_1d, transportation_capacity[w in W, q in W, t in sim_T, s in N],
        y_send_wqt[w,q,t,s] <= transport_capacities[w,q])

    # constraint saying the amount sent from w to q is the same as received at q from w
    @constraint(model_1d, send_receive[w in W, q in W, t in sim_T, s in N], 
        y_send_wqt[w,q,t,s] == y_receive_wqt[q,w,t,s])


    #Storage one day before transfer, initial stock in t=1
    @constraint(model_1d, storage_t1[w in W, q in W, s in N],
        y_send_wqt[w,q,1,s] <= initial_stock[w])

    @constraint(model_1d, storage_before_transfer[w in W, q in W, t in sim_T[2:end], s in N],
        y_send_wqt[w,q,t,s] <= z_wt[w,t-1,s])
        
    # Solve the model
    optimize!(model_1d)

    # Print the solution
    if termination_status(model_1d) == MOI.OPTIMAL
        println("Optimal solution found")
        for t in sim_T
            println("Period ", t)
            for w in W
                println("Warehouse ", w)
                @printf("Order: %0.3f\n", sum(value(x_wt[w,t,s]) for s in N))
                @printf("Storage: %0.3f\n", sum(value(z_wt[w,t,s]) for s in N))
                for q in W
                    if q != w
                        @printf("Send to warehouse %i: %0.3f\n", q, sum(value(y_send_wqt[w,q,t,s]) for s in N))
                        @printf("Receive from warehouse %i: %0.3f\n", q, sum(value(y_receive_wqt[w,q,t,s]) for s in N))
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

Make_Stochastic_here_and_now_decision(p_wt, 1000, 50)