using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using UnitCommitment

Pkg.activate(joinpath(@__DIR__, "."))

using Plots, JSON

data = JSON.parsefile(joinpath(@__DIR__, "..", "test", "data", "test.json"))

T 			= length(data["prices"])
prods 		= Array{Float64}(undef, T)
prods2 		= Array{Float64}(undef, T)
prods3 		= Array{Float64}(undef, T)
prices 		= Float64.(data["prices"])
u_turbine 	= ones(T) * data["max_turbine_capacity"]
u_pump		= ones(T) * data["max_pump_capacity"]
q 			= Float64.(data["inflows"])
eff_pump 	= data["eff_pump"]
x_max 		= data["x_max"]
x_init		= data["x_init"]
x_end		= data["x_end"]


prods = UnitCommitment.optim_pump_storage(prods, prices, u_turbine, u_pump, q, eff_pump, x_max, x_init, x_end)

prods2 = Array{Float64}(undef, T)
prods2 = UnitCommitment.optim_pump_storage_jump(prods2, prices, u_turbine, u_pump, q, eff_pump, x_max, x_init, x_end)

prods3 = UnitCommitment.old_optim_pump_storage(prods3, prices, u_turbine, u_pump, q, eff_pump, x_max, x_init, x_end)


stock = x_init .- cumsum(q .- prods)
stock2 = x_init .+ cumsum(q .- prods2)
stock3 = x_init .+ cumsum(q .- prods3)

p1 = plot(prods, label="Production", legend=:topleft)
p1 = plot!(twinx(), prices, label="Price", color=:red)

p2 = plot(stock, label="Stock")
p2 = plot!(repeat([x_max], size(stock,1)), label="max")

p3 = plot(prods2, label="Production", legend=:topleft)
p3 = plot!(twinx(), prices, label="Price", color=:red)

p4 = plot(stock2, label="Stock")
p4 = plot!(repeat([x_max], size(stock2,1)), label="max")


p5 = plot(prods3, label="Production", legend=:topleft)
p5 = plot!(twinx(), prices, label="Price", color=:red)

p6 = plot(stock3, label="Stock")
p6 = plot!(repeat([x_max], size(stock3,1)), label="max")


#plot(p1, p2, layout=(2,1), size = (800,600))

#plot(p1, p2, p3, p4, layout=(4,1), size = (800,800))
plot(p1, p2, p5, p6, layout=(4,1), size = (800,800))
