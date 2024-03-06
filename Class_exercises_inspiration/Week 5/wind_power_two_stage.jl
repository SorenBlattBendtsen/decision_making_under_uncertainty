#Import packages
using JuMP
using Gurobi
using Printf
# Install by running
# import Pkg
# Pkg.add("MathOptFormat")
# on the commandline
using MathOptFormat


scenarios = [1,2]
lamdba_d = 20
lambda_up = 14
lambda_dw = 27
p_max = 100
w = [82.5, 24.0] # Wind power  per scenario
pi = [0.45, 0.55] #Scenario probability

#Declare Gurobi model
two_stage_wind_power = Model(with_optimizer(Gurobi.Optimizer))

#Definition of variables with lower bound 0
@variable(two_stage_wind_power, 0<=p_d)
@variable(two_stage_wind_power, 0<=p_up[s in scenarios])
@variable(two_stage_wind_power, 0<=p_dw[s in scenarios])




#Maximize profit
@objective(two_stage_wind_power, Max, lamdba_d*p_d + sum(pi[s]*(lambda_up*p_up[s] - lambda_dw*p_dw[s]) for s in scenarios))

#Max production
@constraint(two_stage_wind_power, max_production, p_d <= p_max)

#Deviation from offer
@constraint(two_stage_wind_power, deviation[s in scenarios], w[s] - p_d == p_up[s] - p_dw[s])

optimize!(two_stage_wind_power)

#Output model to an LP file for debugging
lp_model = MathOptFormat.LP.Model()
MOI.copy_to(lp_model, backend(two_stage_wind_power))
MOI.write_to_file(lp_model, "model.lp")

if termination_status(two_stage_wind_power) == MOI.OPTIMAL
    println("Optimal solution found")

    println("Variable values:")
    @printf "p_d: %0.3f\n" value.(p_d)

    println()
    for s in scenarios
        println("Scenario" , s)
        @printf("p_up%i: %0.3f\n",s, value.(p_up[s]))
        @printf("p_dw%i: %0.3f\n",s, value.(p_dw[s]))
        println()
    end

    @printf "\nObjective value: %0.3f\n\n" objective_value(two_stage_wind_power)

else
    error("No solution.")
end
