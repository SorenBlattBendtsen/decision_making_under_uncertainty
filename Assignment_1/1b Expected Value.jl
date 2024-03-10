#Import packages
using JuMP
using Gurobi
using Printf
using Distributions

# Function that takes the inital price and calculates the expected price for 1000 samples using sample_next
function calculate_expected_prices(price)
    num_samples = 1000
    expected_price = zeros(number_of_warehouses)  # Initialize expected prices matrix
    for w in 1:number_of_warehouses
        initial_prices = price[w]  # Extract initial prices as a scalar or vector
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

    
# Function for Expected Value optimization problem
function make_EV_here_and_now_decision(price)
    expected_price = calculate_expected_prices(price)
    
    # Define the model
    model_ev = Model(Gurobi.Optimizer)

   # Variables
   @variable(model_ev, 0 <= x_wt[w in W, t in sim_T]) # Coffee ordered at warehouse w in period t
   @variable(model_ev, 0 <= z_wt[w in W, t in sim_T]) # Coffee stored at warehouse w in period t
   @variable(model_ev, 0 <= y_send_wqt[w in W, q in W, t in sim_T]) # Coffee sent from warehouse w to warehouse q in period t
   @variable(model_ev, 0 <= y_receive_wqt[w in W, q in W, t in sim_T]) # Coffee received at warehouse w from warehouse q in period t
   @variable(model_ev, 0 <= m_wt[w in W, t in sim_T]) # Missing demand at warehouse w in period t


    @objective(model_ev, Min, sum(price[w] * x_wt[w,t] for w in W, t in sim_T[1]) 
                                + sum(cost_miss[w] * m_wt[w,t] for w in W, t in sim_T[1]) 
                              + sum(cost_tr[w,q] * y_send_wqt[w,q,t] for w in W, q in W, t in sim_T[1]) 
                              +sum(expected_price[w] * x_wt[w,t] for w in W, t in sim_T[2]) 
                              + sum(cost_miss[w] * m_wt[w,t] for w in W, t in sim_T[2]) 
                            + sum(cost_tr[w,q] * y_send_wqt[w,q,t] for w in W, q in W, t in sim_T[2]))

    # Define the constraints
   # Demand hour 1
   @constraint(model_ev, demand_t0[w in W],
   x_wt[w,1] - z_wt[w,1] + initial_stock[w]
   + sum(y_receive_wqt[w,q,1] for q in W)
   - sum(y_send_wqt[w,q,1] for q in W)
   + m_wt[w,1] == demand[w,1])

    # Demand hour t
    @constraint(model_ev, demand_t[w in W, t in sim_T[2:end]],
    x_wt[w,t] - z_wt[w,t] + z_wt[w,t-1]
    + sum(y_receive_wqt[w,q,t] for q in W) 
    - sum(y_send_wqt[w,q,t] for q in W)
    + m_wt[w,t] == demand[w,t])

    # Storage capacity
    @constraint(model_ev, storage_capacity[w in W, t in sim_T],
    z_wt[w,t] <= warehouse_capacities[w])

    # Transportation capacity
    @constraint(model_ev, transportation_capacity[w in W, q in W, t in sim_T],
    y_send_wqt[w,q,t] <= transport_capacities[w,q])

    # Constraint saying the amount sent from w to q is the same as received at q from w
    @constraint(model_ev, send_receive[w in W, q in W, t in sim_T], 
        y_send_wqt[w,q,t] == y_receive_wqt[q,w,t])


    #Storage one day before transfer, initial stock in t=1
    @constraint(model_ev, storage_t1[w in W, q in W],
    y_send_wqt[w,q,1] <= initial_stock[w])

    @constraint(model_ev, storage_before_transfer[w in W, q in W, t in sim_T[2:end]],
    y_send_wqt[w,q,t] <= z_wt[w,t-1])
                
    # Solve the model
    optimize!(model_ev)

    # Return the optimal solution
    if termination_status(model_ev) == MOI.OPTIMAL
        println("Optimal solution found")
        # Print results for Day 1
        # println("Day 1:")
        # for w in W
        #     println("Warehouse ", w)
        #     @printf("Price: %0.3f\n", value(price[w]))
        #     @printf("Order: %0.3f\n", value(x_wt[w,1]))
        #     @printf("Storage: %0.3f\n", value(z_wt[w,1]))
        #     for q in W
        #         if q != w
        #             @printf("Send to warehouse %i: %0.3f\n", q, value(y_send_wqt[w,q,1]))
        #             @printf("Receive from warehouse %i: %0.3f\n", q, value(y_receive_wqt[w,q,1]))
        #         end
        #     end
        # end
             
        # println("-------")

        # # Print results for Day 2
        # println("Day 2:")
        # for w in W
        #     println("Warehouse ", w)
        #     @printf("Expected Price: %0.3f\n", value(expected_price[w]))
        #     @printf("Order: %0.3f\n", value(x_wt[w,2]))
        #     @printf("Storage: %0.3f\n", value(z_wt[w,2]))
        #     for q in W
        #         if q != w
        #             @printf("Send to warehouse %i: %0.3f\n", q, value(y_send_wqt[w,q,2]))
        #             @printf("Receive from warehouse %i: %0.3f\n", q, value(y_receive_wqt[w,q,2]))
        #         end
        #     end
        # end
        # Extract values for Day 1
        # prices_day_1 = value.(price[:, 1])
        # orders_day_1 = value.(x_wt[:, 1])
        # storage_day_1 = value.(z_wt[:, 1])
        # send_receive_day_1 = value.(y_send_wqt[:, :, 1])

        # # Extract values for Day 2
        # expected_prices_day_2 = value.(expected_price[:, 2])
        # orders_day_2 = value.(x_wt[:, 2])
        # storage_day_2 = value.(z_wt[:, 2])
        # send_receive_day_2 = value.(y_send_wqt[:, :, 2])

        system_cost = objective_value(model_ev)
        println("Total system cost: ", system_cost)

        # Return the relevant results
        return system_cost
        else
        println("No solution found")
        return nothing
        
    end 

end
