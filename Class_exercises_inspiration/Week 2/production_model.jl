
using JuMP
using Gurobi
using Printf


#Scenarios
S = collect(1:5)
#Products
P = collect(1:10)
#Machine types
M = collect(1:4)

#Revenue per product, call revenue[p]
revenue = [346 308 397 371 326 336 334 321 319 350]

#Low revenue per excess product, call revenue_low[p]
revenue_low = [112 114 134 145 112 123 100 98 111 120]

#Import cost per missing product, call import_cost[p]
import_cost = [546 693 661 662 753 862 598 503 735 648]

#Production cost per product p, call production_cost[p]
production_cost = [138 181  149 177 187 131 149 178 147 189]

#Investment cost per machine of machine type m, call investment_cost[m]
investment_cost = [10000 50000 30000 75000]

#Demand per product and scenario, call demand[p][s]
demand = [[20126.87, 22451.75, 11361.12, 11614.00, 23060.0, ],
    [19054.87, 28412.50, 19471.96, 17869.09, 16620.86 ],
    [14041.74, 9756.50, 18546.15, 14949.18, 18375.47 ],
    [13781.74, 16554.50, 20321.81, 19499.82, 16911.36 ],
    [22480.09, 16860.00, 18547.81, 10698.00, 14225.81 ],
    [24724.00, 22070.50, 18218.31, 14952.73, 12174.78 ],
    [18387.22, 18425.75, 12736.15, 15351.45, 21257.61 ],
    [17123.00, 19015.75, 19993.96, 15251.45, 19865.64 ],
    [15363.43, 23237.50, 17813.42, 21893.36, 19234.42 ],
    [18530.96, 23853.25, 14233.27, 25578.64, 14084.92 ]]

# Probability per scenario, call probabilities[s]    
probabilities = [0.23, 0.04, 0.26, 0.11, 0.36]

#Compatibility between product and machine type, 1 if compatible, 0 otherwise
#Call machine_compatibility[p][m]
machine_compatibility = [[	1	1	0	0	],
[	1	1	0	0	],
[	0	1	1	0	],
[	0	0	1	1	],
[	0	0	1	1	],
[	1	1	1	1	],
[	0	0	0	1	],
[	1	1	0	1	],
[	1	1	0	0	],
[	0	1	1	0	]]

#Production time per product p, call production_time[p]
production_time = [0.2 0.2 0.1 0.3 0.2 0.4 0.15 0.6 0.4 0.6]

#Working time available per machine of machine type m, call working_time[m]
working_time = [2080 2080 2080 2080]

#Maximum number of machines
maximum_machines = 100
#Maximum budget
maximum_budget = 20000000


