#Import packages
using JuMP
using Gurobi
using Printf

#Declare model with Gurobi solver
model_carl = Model(Gurobi.Optimizer)

#Declare INTEGER variables with lower bound 0 and upper bound
@variable(model_carl, 0<=xC<=50, Int)
@variable(model_carl, 0<=xS<=200, Int)
#Declare BINARY variables with lower bound 0
@variable(model_carl, 0<=yS, Bin)
@variable(model_carl, 0<=yC, Bin)
#Declare continuous variables
@variable(model_carl, pR)
@variable(model_carl, pW)

#Declare maximization of profits objective function
@objective(model_carl, Max, (400-200)xC + pR + pW - 1000yC)
#Constraint on available acres
@constraint(model_carl, Acres, xC + 0.2xS <= 72)
#Constraint on maximum working hours
@constraint(model_carl, WorkingHours, 150xC + 25xS <= 10000)
#Minimum of 100 sheep constraint for wholesale price
@constraint(model_carl, AtLeast100Sheep1, xS - 100yS >= 0)
@constraint(model_carl, AtLeast100Sheep2, xS - 200yS <= 0)
#Maximum of 10 cows without milk machine
@constraint(model_carl, MilkMachine1, xC - 10-40yC<= 0)
@constraint(model_carl, MilkMachine2, xC - 11yC >= 0)

#big-M constraints "if sheep less than 100 --> retail price, otherwise --> wholesale price"
@parameter(model_carl, M = 40000)
@constraint(model_carl, retail1, -M*yS +30xS <= pR)
@constraint(model_carl, retail2, M*yS +30xS >= pR)
@constraint(model_carl, retail3, -M*(1-yS) <= pR)
@constraint(model_carl, retail4, M*(1-yS) >= pR)
@constraint(model_carl, wholesale1, -M*(1-yS) +40xS <= pW)
@constraint(model_carl, wholesale2, M*(1-yS) +40xS >= pW)
@constraint(model_carl, wholesale3, -M*yS <= pW)
@constraint(model_carl, wholesale4, M*yS >= pW)

#Optimize model
optimize!(model_carl)

#Check if optimal solution was found
if termination_status(model_carl) == MOI.OPTIMAL
    println("Optimal solution found")

    #Print out variable values and objective value
    println("Variable values:")
    @printf "xC: %0.3f\n" value.(xC)
    @printf "xS: %0.3f\n" value.(xS)
    @printf "yS: %0.3f\n" value.(yS)
    @printf "yC: %0.3f\n" value.(yC)
    @printf "\nObjective value: %0.3f\n\n" objective_value(model_carl)

else
    error("No solution.")
end
