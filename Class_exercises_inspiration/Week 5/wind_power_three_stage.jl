#Import packages
using JuMP
using Gurobi
using Printf
# Install by running
# import Pkg
# Pkg.add("MathOptFormat")
# on the commandline
using MathOptFormat


scenarios = [1,2,3,4]
lambda_d = 20
lambda_a = 18
lambda_up = 14
lambda_dw = 27
p_max = 100
w = [100, 50, 40, 0, ] # Wind power  per scenario
pi = [0.2925 , 0.1575 , 0.33 , 0.22] #Scenario probability

#Declare Gurobi model
three_stage_wind_power = Model(with_optimizer(Gurobi.Optimizer))

#Definition of variables with lower bound 0
@variable(three_stage_wind_power, 0<=p_d)
@variable(three_stage_wind_power, p_a[s in scenarios])
@variable(three_stage_wind_power, 0<=p_up[s in scenarios])
@variable(three_stage_wind_power, 0<=p_dw[s in scenarios])




#Maximize profit
@objective(three_stage_wind_power, Max, lambda_d*p_d + sum(pi[s]*(lambda_a*p_a[s] + lambda_up*p_up[s] - lambda_dw*p_dw[s]) for s in scenarios))

#Max production
@constraint(three_stage_wind_power, max_production, p_d <= p_max)

#Max production on adjustment market
@constraint(three_stage_wind_power, max_production_adjust[s in scenarios], p_d + p_a[s] <= p_max)

#Deviation from offer
@constraint(three_stage_wind_power, deviation[s in scenarios], w[s] - p_d - p_a[s]== p_up[s] - p_dw[s])

#Non-anticipativity
@constraint(three_stage_wind_power, non_anticipativity1, p_a[1] == p_a[2])
@constraint(three_stage_wind_power, non_anticipativity2, p_a[3] == p_a[4])


optimize!(three_stage_wind_power)

#Output model to an LP file for debugging
lp_model = MathOptFormat.LP.Model()
MOI.copy_to(lp_model, backend(three_stage_wind_power))
MOI.write_to_file(lp_model, "model.lp")

if termination_status(three_stage_wind_power) == MOI.OPTIMAL
    println("Optimal solution found")

    println("Variable values:")
    @printf "p_d: %0.3f\n" value.(p_d)

    println()
    for s in scenarios
        println("Scenario" , s)
        @printf("p_a%i: %0.3f\n",s, value.(p_a[s]))
        @printf("p_up%i: %0.3f\n",s, value.(p_up[s]))
        @printf("p_dw%i: %0.3f\n",s, value.(p_dw[s]))
        println()
    end

    @printf "\nObjective value: %0.3f\n\n" objective_value(three_stage_wind_power)

else
    error("No solution.")
end
