using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))
using UnitCommitment
Pkg.activate(joinpath(@__DIR__, "."))
using Plots, JSON

data = JSON.parsefile(joinpath(@__DIR__, "..", "test", "data", "test1.json"))
# data = JSON.parsefile(joinpath(@__DIR__, "", "test", "data", "test2.json"))

T 			= length(data["prices"])
prices 		= Float64.(data["prices"])
u_turbine 	= ones(T) * data["max_turbine_capacity"]
u_pump		= ones(T) * data["max_pump_capacity"]
q 			= Float64.(data["inflows"])
eff_pump 	= data["eff_pump"]
x_max 		= data["x_max"]
x_init		= data["x_init"]
x_end		= data["x_end"]

n_grid = 100

prods_new 		= Array{Float64}(undef, T)
prods_jump 		= Array{Float64}(undef, T)
prods_old 		= Array{Float64}(undef, T)

@info ">>> NEW"
prods_new = UnitCommitment.optim_pump_storage(prods_new, prices, u_turbine, u_pump, q, eff_pump, x_max, x_init, x_end, n_grid)
@info ">>> JUMP"
prods_jump = UnitCommitment.optim_pump_storage_jump(prods_jump, prices, u_turbine, u_pump, q, eff_pump, x_max, x_init, x_end)
@info ">>> OLD"
prods_old = UnitCommitment.old_optim_pump_storage(prods_old, prices, u_turbine, u_pump, q, eff_pump, x_max, x_init, x_end, n_grid)


stock_new = x_init .+ cumsum(q .- prods_new)
stock_jump = x_init .+ cumsum(q .- prods_jump)
stock_old = x_init .+ cumsum(q .- prods_old)

p1 = plot(prods_new, label="new", color = :blue, legend=:top);
p1 = plot!(prods_jump, label="jump", color=:green);
p1 = plot!(prods_old, label="old", color=:orange);
p1 = plot!(twinx(), prices, label="price", color=:red, legend=false);

p2 = plot(stock_new, label="new",color = :blue);
p2 = plot!(stock_jump, label="jump",color=:green);
p2 = plot!(stock_old, label="old",color=:orange);
p2 = hline!([x_max], label="max", color=:red);

plot(p1, p2, layout=(2,1), size = (800,800))
