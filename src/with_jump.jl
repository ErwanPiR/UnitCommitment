
using JuMP, GLPK

function optim_pump_storage_jump(prods, price, u_turbine, u_pump, q, eff_pump, x_max, x_init, x_end)

    n = length(price)
    
    m = Model(optimizer_with_attributes(GLPK.Optimizer))
    
    pump = @variable(m, [1:n])
    turbine = @variable(m, [1:n])
    set_upper_bound.(pump, u_pump)
    set_lower_bound.(pump, 0)
    set_upper_bound.(turbine, u_turbine)
    set_lower_bound.(turbine, 0)

    for i=1:n
        @constraint(m, x_init + sum(q[1:i] .- turbine[1:i] .+ pump[1:i]) <= x_max)
        @constraint(m,  x_init + sum(q[1:i] .- turbine[1:i] .+ pump[1:i]) >= 0)
    end

    @constraint(m,  x_init + sum(q[1:n] .- turbine[1:n] .+ pump[1:n]) == x_end)

    ## Ne connait pas la contrainte que sur un timestep, capa de pompage / turbinage = capa max de la pompe / turbine

    @objective(m, Max, sum(price .* turbine .- price .* pump .* eff_pump) )

    optimize!(m)

    prods .= value.(turbine) .- value.(pump)

    return prods
end

