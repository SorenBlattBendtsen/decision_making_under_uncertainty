
using CSV
using DataFrames
using Plots
using Statistics
using Random
using Printf

prob = [0.25 0.05 0.3 0.25 0.15]
scenarios = [10	7 15 20	13]

initfig = bar(scenarios,prob, xlim=(6,22), ylim=(0,0.5), xlabel="Wind speed" , ylabel="Probability", label=["Scenario 1" "Scenario 2" "Scenario 3" "Scenario 4" "Scenario 5" ])
savefig(initfig, "init-ffs.pdf")


D = zeros(Float64, length(scenarios), length(scenarios))
for s = 1:length(scenarios)
    for s2 = 1:length(scenarios)
        D[s,s2] = abs(scenarios[s] - scenarios[s2])
    end
end
probabilities = repeat([1.0/length(scenarios)], 1, length(scenarios))[1,:]
include("fast-forward-selection.jl")
result = FastForwardSelection(D,probabilities,3)
print(result[1])
print(result[2])

plotdata = zeros(Float64, length(result[2]))
for i = 1:length(result[2])
    plotdata[i] = scenarios[result[2][i]]
end
