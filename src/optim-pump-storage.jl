
using Interpolations

"""

optim_pump_storage(splant, prods, costs, price, umarginal_cost, u_turbine, u_pump, q, eff_pump, x_max, x_init, x_end)

# Arguments
- `splant` : 
- `prods` : 
- `costs` : 
- `price` : 
- `umarginal_cost` : 
- `x_min` : minimum capacity of the reservoir
- `x_max` : maximum capacity of the reservoir
- `u_turbine`: maximum pumping capacity : pump capacity 
- `u_pump` : maximum releasing capacity : turbine capacity
- `eff_pump` : pumping efficiency
- `q` : inflow

"""
@enum Controls idle pump turbine

# Computes the capacity of water pumped or released u according to constraints of the reservoir, turbine and pump capacities
function u(x::Float64, c::Controls, u_max::Float64, x_max::Float64)
    (c == idle) && return 0
    (c == pump) && return min(u_max, (x_max - x))
    (c == turbine) && return -min(u_max, x)
end

# # Computes the capacity of the reservoir after action c
function x_fwd(x::Float64, q::Float64, c::Controls, u_max::Float64, x_max::Float64) 
    xf = x + q + u(x, c, u_max, x_max)
    return min(xf, x_max)
end

# Computes the cost of pumping or releasing water
C(x::Float64, pr::Float64, c::Controls, u_max::Float64, x_max::Float64, eff::Float64) = 
    pr * u(x, c, u_max, x_max) * eff



function optim_pump_storage(prods, price, u_turbine, u_pump, q, eff_pump, x_max, x_init, x_end)

    # x_min = 0. # minimum capacity of the reservoir
    # x_max = 1. # maximum capacity of the reservoir
    # u_turbine = 2/3 # maximum pumping capacity : pump capacity 
    # u_pump = 1/3 # maximum releasing capacity : turbine capacity
    # eff_pump = 2/3 # pumping efficiency
    # eff_turbine = 3/4 # turbine efficiency
    # q = 0. # inflow
    # T = 1
    # n_grid = 4 # number of grid points
    # price = [40.]
    # profits = [-100., -5., 0., 10.]
    # x = (x_max : -(x_min + (x_max-x_min)/(n_grid-1)): x_min) # grid points

    x_min = 0 # Minimum capacity of the reservoir
    eff_turbine = 1. # Turbine efficiency
    eff_pump = 1. # Pumping efficiency
    n_grid = 100 # Number of grid points
    x = (x_min : (x_min + (x_max-x_min)/(n_grid-1)): x_max) # Grid points
    T = length(price)

    profits = zeros(n_grid) # Saves profits from one t to the next
    policies = Array{Controls}(undef, n_grid, T) # Saves policies for forward reconstruction

    profits .= -2 * max.(x_end .- x, 0) * price[T] # Initial profits

    x_interp = linear_interpolation(x, 1:length(x)) # Linear interpolation to get fast access to the index of the grid point corresponding to the capacity of the reservoir

    # Backward loop
    for t in T:-1:1
        # Computes maximum profit and corresponding policy for each grid point
        for (i, xi) in enumerate(x)
            
            xf = x_fwd(xi, q[t], idle, 0., x_max)
            profit_max = profits[round(Int, x_interp(xf))]
            policies[i, t] = idle
            
            xf = x_fwd(xi, q[t], pump, u_pump[t], x_max)
            profit_pump =  profits[round(Int, x_interp(xf))] .- C(xi, price[t], pump, u_pump[t], x_max, 1/eff_pump)
            if profit_pump > profit_max
                profit_max = profit_pump
                policies[i, t] = pump
            end

            xf = x_fwd(xi, q[t], turbine, u_turbine[t], x_max)
            profit_turbine =  profits[round(Int, x_interp(xf))] .- C(xi, price[t], turbine, u_turbine[t], x_max, eff_turbine)
            if profit_turbine > profit_max
                profit_max = profit_turbine
                policies[i, t] = turbine
            end

            profits[i] = profit_max

        end
    end
    
    x_plant = zeros(T+1) # Capacity of the reservoir
    x_plant[1] = x_init # Inital state
    u_plant = zeros(T+1) # Capacity of water pumped or released

    # Forward loop
    for t in 2:T 
        # Retrieve the control c corresponding to the maximum profit
        c = policies[round(Int, x_interp(x_plant[t-1])), t-1]
        
        # Retrieve the maximum capacity of pump / turbine according to c
        ut = (c==pump)*(u_pump[t-1]) + (c==turbine)*(u_turbine[t-1])

        # Computes the capacity of the reservoir 
        x_plant[t] = x_fwd(x_plant[t-1], q[t-1], c, ut, x_max)

        # Computes the capacity of water pumped or released
        u_plant[t] = u(x_plant[t-1], c, ut, x_max)
    end

    # Saves results
    prods .= -u_plant[1:T] # turbining is positive
    #costs[splant] .= Utils.calculate_marginal_cost(umarginal_cost, prods[splant])

    return prods

end
