# File: 02435_multistage_policy.jl

using JuMP
using Gurobi
using Distributions
using Random
# Make sure this script is saved in the same directory as your other .jl files
include("V2_02435_multistage_problem_data.jl")
include("fast-forward-selection.jl")

# Placeholder function for generating and reducing scenarios
# This part needs to be filled in with your logic based on your chosen approach
# Function for generating and reducing scenarios
function generate_and_reduce_scenarios(current_prices, num_samples, reduced_samples, lookahead)
    # Initialize matrices to store the generated and reduced scenarios
    all_scenarios = Array{Float64,3}(undef, num_samples, length(current_prices), lookahead)
    reduced_price_scenarios = Array{Float64,3}(undef, reduced_samples, length(current_prices), lookahead)
    reduced_probabilities = zeros(reduced_samples, lookahead)
    
    for t in 1:lookahead
        for s in 1:num_samples
            for w in 1:length(current_prices)
                #The function sample_next needs to be defined properly to simulate future prices based on your specific model.
                all_scenarios[s,w,t] = sample_next(current_prices[w])
            end
        end

        D = zeros(num_samples, num_samples)
        for i in 1:num_samples
            for j in 1:num_samples
                D[i,j] = sqrt(sum((all_scenarios[i,:,t] - all_scenarios[j,:,t]).^2))
            end
        end

        probabilities, selected_indices = FastForwardSelection(D, fill(1.0 / num_samples, num_samples), reduced_samples)

        for i in 1:reduced_samples
            reduced_price_scenarios[i,:,t] = all_scenarios[selected_indices[i],:,t]
            reduced_probabilities[i,t] = probabilities[i]
        end
    end
    return reduced_price_scenarios, reduced_probabilities
end

function make_multistage_here_and_now_decision(number_of_sim_periods, tau, current_stock, current_prices)
    number_of_warehouses, W, cost_miss, cost_tr, warehouse_capacities, transport_capacities, initial_stock, number_of_simulation_periods, sim_T, demand_trajectory = load_the_data()

    # Determine lookahead based on current period tau
    lookahead = min(5-tau+1, 3)  # Adjust this based on performance requirements
    # 5 - tau + 1: This part calculates how many days are left in the planning horizon from the current day (tau). For example, if today is day 1 (tau = 1), then there are 5 days left in the horizon. If today is day 3 (tau = 3), then there are 3 days left.
    # min(5 - tau + 1, 3): This function then takes the minimum between the number of days left and 3. The purpose of using 3 here is to limit the maximum lookahead to 3 days, regardless of how many days are actually left. This could be due to computational constraints or the belief that looking more than 3 days ahead does not significantly improve decision quality due to increasing uncertainty.

    # Sample a limited number of future price scenarios and their probabilities
    # Depending on your computational constraints, adjust the number of samples and how they're reduced
    num_samples = 100  # Reduced number of samples due to computational constraints
    reduced_samples = 10  # Further reduction for scenario handling within constraints

    # Generate price scenarios for lookahead days (this should include scenario reduction logic)
    # Placeholder for actual scenario generation logic
    price_scenarios, probabilities = generate_and_reduce_scenarios(current_prices, num_samples, reduced_samples, lookahead)

    # Set up the optimization problem
    model = Model(Gurobi.Optimizer)
    # This variable represents the amount of coffee (or other goods) ordered by warehouse w (where w is an index running through all warehouses in the set W) at time t (where t ranges from 1 to the value of lookahead).
    # is a decision variable, meaning the model will determine the optimal quantity of coffee to order for each warehouse at each time step within the lookahead horizon to minimize the total cost.
    @variable(model, 0 <= x[w in W, t in 1:lookahead])
    # This variable indicates the storage level (amount of coffee stored) at warehouse w at time t.
    # Like x, z is a decision variable optimized by the model, respecting the constraints like storage capacities and ensuring that coffee is available to meet demand or to be sent to other warehouses.
    @variable(model, 0 <= z[w in W, t in 1:lookahead])
    # This variable represents the amount of coffee sent from warehouse w to another warehouse q at time t.
    # The indices w and q both range over all warehouses (but typically, a warehouse will not send coffee to itself, hence w should not equal q), and t ranges over the lookahead period.
    # This is another decision variable that the model will optimize, taking into account transportation limits and costs.
    @variable(model, 0 <= y_send[w in W, q in W, t in 1:lookahead])
    # In contrast to y_send, this variable denotes the amount of coffee received by warehouse w from warehouse q at time t.
    # While the physical reality would dictate that y_receive is closely tied to y_send (i.e., what one warehouse sends, another receives), they are kept separate in the model for clarity and to enforce the symmetry through constraints.
    @variable(model, 0 <= y_receive[w in W, q in W, t in 1:lookahead])
    # This variable represents the missing amount of coffee (unmet demand) at warehouse w at time t.
    # Essentially, this is the shortfall: the amount by which the warehouse fails to meet the demand. The objective function typically penalizes missing amounts to ensure demands are met as closely as possible.
    @variable(model, 0 <= m[w in W, t in 1:lookahead])

    # Add your constraints similar to those from previous steps, adjusted for the lookahead
    # The constraints are now looped over each time stage t in the lookahead, which can be up to 5 days ahead. This replaces the separate demand_1, demand_2, etc., constraints in the two-stage model with a single, unified loop that applies the same logic at each stage.
    # For t=1, the model uses initial_stock[w] since there's no preceding storage level. For t > 1, it uses z[w, t-1] to represent the storage level from the end of the previous day.
    # demand_trajectory[w,t] should be used instead of demand[w,1] or demand[w,2], reflecting the specific demand for each warehouse w at each time t.
    # The model maintains consistency between the amount sent from one warehouse to another and the amount received, as well as ensuring that any coffee sent on a given day was already in stock the day before.
    # This pattern keeps the core logic of your constraints consistent while allowing them to be applied across more stages. It means your model can dynamically adapt to different lookahead values without needing separate code for each potential stage.
    # Loop over all time stages for dynamic constraints
for t in 1:lookahead
    # First stage constraints when t = 1, subsequent stages otherwise
    # Demand balance
    @constraint(model, demand_balance[w in W, t],
                x[w,t] - z[w,t] + (t == 1 ? initial_stock[w] : z[w,t-1])
                + sum(y_receive[w,q,t] for q in W)
                - sum(y_send[w,q,t] for q in W)
                + m[w,t] == demand_trajectory[w,t])

    # Storage capacity
    @constraint(model, storage_capacity[w in W, t],
                z[w,t] <= warehouse_capacities[w])

    # Transportation capacity
    @constraint(model, transportation_capacity[w in W, q in W, t],
                y_send[w,q,t] <= transport_capacities[w,q])

    # Consistency between sending and receiving
    @constraint(model, send_receive_consistency[w in W, q in W, t],
                y_send[w,q,t] == y_receive[q,w,t])

    # Storage from previous day before transfer, accounting for t=1
    @constraint(model, storage_before_transfer[w in W, q in W, t],
                y_send[w,q,t] <= (t == 1 ? initial_stock[w] : z[w,t-1]))
end
#Alternative contraint representation if previous one is not working
#for t in 1:lookahead
#@constraint(model, x[w in W] - z[w] + (t == 1 ? initial_stock[w] : z[w,t-1]) + sum(y_receive[w,q,t] for q in W) - sum(y_send[w,q,t] for q in W) + m[w,t] == demand_trajectory[w,t] for w in W)
#@constraint(model, z[w,t] <= warehouse_capacities[w] for w in W)
#@constraint(model, y_send[w,q,t] <= transport_capacities[w,q] for w in W, q in W)
#@constraint(model, y_send[w,q,t] == y_receive[q,w,t] for w in W, q in W)
#end

    # Objective function should consider the expected costs over the lookahead period
    @objective(model, Min, sum(cost_miss[w] * m[w,t] for w in W, t in 1:lookahead) + 
                             sum(cost_tr[w,q] * y_send[w,q,t] for w in W, q in W, t in 1:lookahead) + 
                             sum(price_scenarios[s,w,t] * x[w,t] * probabilities[s] for s in 1:reduced_samples, w in W, t in 1:lookahead))

    # Solve the problem
    optimize!(model)

    # Extract here-and-now decisions (for tau = 1) and prepare them for return
    x_decision = value.(x[:,1])
    z_decision = value.(z[:,1])
    send_decision = value.(y_send[:,:,1])
    receive_decision = value.(y_receive[:,:,1])
    m_decision = value.(m[:,1])

    return x_decision, send_decision, receive_decision, z_decision, m_decision
end

