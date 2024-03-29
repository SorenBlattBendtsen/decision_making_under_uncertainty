# 1c, optimality in hindsight (deterministic demand and price)
# import packages
using JuMP
using Gurobi
using Printf

# Function to calculate the optimal solution in hindsight 
function Calculate_OiH_solution(price)
    
    # Define the model
    model_1c = Model(Gurobi.Optimizer)

    # Variables
    @variable(model_1c, 0 <= x_wt[w in W, t in sim_T]) # Coffee ordered at warehouse w in period t
    @variable(model_1c, 0 <= z_wt[w in W, t in sim_T]) # Coffee stored at warehouse w in period t
    @variable(model_1c, 0 <= y_send_wqt[w in W, q in W, t in sim_T]) # Coffee sent from warehouse w to warehouse q in period t
    @variable(model_1c, 0 <= y_receive_wqt[w in W, q in W, t in sim_T]) # Coffee received at warehouse w from warehouse q in period t
    @variable(model_1c, 0 <= m_wt[w in W, t in sim_T]) # Missing demand at warehouse w in period t

    # Define the objective function
    @objective(model_1c, Min, sum(price[w,t] * x_wt[w,t] for w in W, t in sim_T) 
                            + sum(cost_miss[w] * m_wt[w,t] for w in W, t in sim_T) 
                            + sum(cost_tr[w,q] * y_send_wqt[w,q,t] for w in W, q in W, t in sim_T))

    # Define the constraints
    # demand hour 1
    @constraint(model_1c, demand_fulfillment_t0[w in W],
                x_wt[w,1] - z_wt[w,1] + initial_stock[w]
                + sum(y_receive_wqt[w,q,1] for q in W)
                - sum(y_send_wqt[w,q,1] for q in W)
                + m_wt[w,1] == demand[w,1])

    #demand hour t > 1
    @constraint(model_1c, demand_fulfillment[w in W, t in sim_T[2:end]],
                x_wt[w,t] - z_wt[w,t] + z_wt[w,t-1]
                + sum(y_receive_wqt[w,q,t] for q in W) 
                - sum(y_send_wqt[w,q,t] for q in W)
                + m_wt[w,t] == demand[w,t])

    # Storage capacity
    @constraint(model_1c, storage_capacity[w in W, t in sim_T],
                z_wt[w,t] <= warehouse_capacities[w])

    # Transportation capacity
    @constraint(model_1c, transportation_capacity[w in W, q in W, t in sim_T],
                y_send_wqt[w,q,t] <= transport_capacities[w,q])

    # constraint saying the amount sent from w to q is the same as received at q from w
    @constraint(model_1c, send_receive[w in W, q in W, t in sim_T], 
                 y_send_wqt[w,q,t] == y_receive_wqt[q,w,t])


    #Storage one day before transfer, initial stock in t=1
    @constraint(model_1c, storage_1[w in W, q in W],
                y_send_wqt[w,q,1] <= initial_stock[w])
    @constraint(model_1c, storage_2[w in W, q in W, t in sim_T[2:end]],
                y_send_wqt[w,q,t] <= z_wt[w,t-1])
                
    # Solve the model
    optimize!(model_1c)

    # Return the optimal solution
    if termination_status(model_1c) == MOI.OPTIMAL
        println("Optimal solution found")
        system_cost = objective_value(model_1c)
        println("Total system cost: ", system_cost)
    else
        println("No solution found")
    end 
    return system_cost
end
