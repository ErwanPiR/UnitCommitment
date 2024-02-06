using JSON

function build_input(path)
	# Build input data
	data = JSON.parsefile(path)

	return data
end

path = joinpath(@__DIR__, "data", "test.json")
data = build_input(path)

@test in("max_pump_capacity", keys(data))
@test in("max_turbine_capacity", keys(data))

T 			= length(data["prices"])
prods 		= Array{Float64}(undef, T)
prices 		= Float64.(data["prices"])
u_turbine 	= ones(T) * data["max_turbine_capacity"]
u_pump		= ones(T) * data["max_pump_capacity"]
q 			= Float64.(data["inflows"])
eff_pump 	= data["eff_pump"]
x_max 		= data["x_max"]
x_init		= data["x_init"]
x_end		= data["x_end"]


UnitCommitment.optim_pump_storage(prods, prices, u_turbine, u_pump, q, eff_pump, x_max, x_init, x_end)
UnitCommitment.optim_pump_storage_jump(prods, prices, u_turbine, u_pump, q, eff_pump, x_max, x_init, x_end)
