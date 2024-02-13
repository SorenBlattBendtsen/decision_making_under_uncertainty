# import libraries
using JuMP
using Gurobi

# create model
model = Model(Gurobi.Optimizer)

# define variables
@variables model begin
    0 <= xC <= 50, Int
    0 <= xS <= 200, Int
    0 <= yS, Bin
    0 <= yC, Bin
end 

# define objective function
@objective(model, Max, (400-200)xC + (70-30)xS - 1000yC)

# define constraints
@constraint(model, Acres, xC + 0.2xS <= 72)
@constraint(model, workingHours, 150xC + 25xS <= 10000)

# solve model
optimize!(model)