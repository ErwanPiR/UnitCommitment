
function old_optim_pump_storage(prods, price, capa_turb, capa_pump, inflow, eff_pump, volume_max, volume_init, volume_end)


    eff_pump = 1.0
    
    T = size(price, 1)

    # parameters
    volume_min = 0
    n_grid = 100
    #inflow = zeros(T)

    ### QUESTION 
    # pourquoi capa_in & capa_out normalisés dans g_in et h_out alors que deja normalisé dans le calcul
    # Idem inflow 
    # erreur calcul capa_out ? 
    # pourquoi divisé par eff_pump ? 
    # calcul des poids pour l'interpolation linéaire 
    # ind_pump = (prods[splant] .< 0);  prods[splant][ind_pump] .= prods[splant][ind_pump] ./ eff_pump : ???
  
    # grid
    grid = collect(1:n_grid)
    grid_step = Float64(floor(volume_max/(n_grid-1); digits = 5))
    grid_abs = volume_min:grid_step:volume_max
  
    # injection and withdrawal capacities
    capa_in = min.(repeat(volume_max.- grid_abs,1,T),repeat(capa_pump',n_grid,1))
    capa_out = min.(repeat(grid_abs.-volume_min,1,T),repeat(capa_turb',n_grid,1))
  
    g_no = repeat(convert(Array{Float64},grid),1,T) + repeat(inflow'/grid_step,n_grid,1)
    g_in = repeat(convert(Array{Float64},grid),1,T) + repeat(inflow'/grid_step,n_grid,1) + capa_in/grid_step
    g_out = repeat(convert(Array{Float64},grid),1,T) + repeat(inflow'/grid_step,n_grid,1) - capa_out/grid_step
  
    g_no = min.(max.(g_no,1),n_grid)
    g_in = min.(max.(g_in,1),n_grid)
    g_out = min.(max.(g_out,1),n_grid)
  
    g_no_dn = convert(Array{Int64},floor.(g_no))
    g_in_dn = convert(Array{Int64},floor.(g_in))
    g_out_dn = convert(Array{Int64},floor.(g_out))
  
    g_no_up = convert(Array{Int64},ceil.(g_no))
    g_in_up = convert(Array{Int64},ceil.(g_in))
    g_out_up = convert(Array{Int64},ceil.(g_out))
  
    w_no_dn =  g_no_up - g_no
    w_in_dn =  g_in_up - g_in
    w_out_dn =  g_out_up - g_out
  
    # Continuation value
    V = zeros(n_grid, T+1)
    V[:,T+1] = -2*max.(volume_end.-grid_abs,0)*price[T] # smoothing losses
  
    ## Backward loop
    for t = T:-1:1
        for i=1:n_grid
            Vno = V[g_no_dn[i,t],t+1].*w_no_dn[i,t] .+ V[g_no_up[i,t],t+1].*(1 .-w_no_dn[i,t])
            Vin = -(1/eff_pump)*capa_in[i,t]*price[t] .+ V[g_in_dn[i,t],t+1].*w_in_dn[i,t] + V[g_in_up[i,t],t+1].*(1 .-w_in_dn[i,t])
            Vout = capa_out[i,t]*price[t] .+ V[g_out_dn[i,t],t+1].*w_out_dn[i,t] .+ V[g_out_up[i,t],t+1].*(1 .-w_out_dn[i,t])
            V[i,t] = max(Vno, Vin, Vout)
        end
    end
  
    ## Forward loop
    prod = zeros(T) # production quantity (pumping > 0 and turbining < 0) during period t
    volume = zeros(T+1) # volume in storage at the beginning of any period t
    volume[1] = volume_init
    
    for t = 1:T
      grid_t = floor((volume[t] - volume_min)/grid_step + 1)
      grid_t = Int64(min.(max.(grid_t,1),n_grid))
  
      capa_in_t = max.(min.(volume_max - volume[t] - inflow[t], capa_in[grid_t,t]),0)
      capa_out_t = max.(min.(volume[t] + inflow[t] - volume_min, capa_out[grid_t,t]),0)
  
      g_no = (volume[t] + inflow[t] - volume_min)/grid_step + 1
      g_in = (volume[t] + inflow[t] + capa_in_t - volume_min)/grid_step + 1
      g_out = (volume[t] + inflow[t] - capa_out_t - volume_min)/grid_step + 1
  
      g_no = min.(max.(g_no,1),n_grid)
      g_in = min.(max.(g_in,1),n_grid)
      g_out = min.(max.(g_out,1),n_grid)
  
      g_no_dn = Int64(floor(g_no))
      g_in_dn = Int64(floor(g_in))
      g_out_dn = Int64(floor(g_out))
  
      g_no_up = Int64(ceil(g_no))
      g_in_up = Int64(ceil(g_in))
      g_out_up = Int64(ceil(g_out))
  
      w_no_dn =  g_no_up - g_no
      w_in_dn =  g_in_up - g_in
      w_out_dn =  g_out_up - g_out
  
      Vno = V[g_no_dn,t+1]*w_no_dn + V[g_no_up,t+1]*(1-w_no_dn)
      Vin = -(1/eff_pump)*capa_in_t*price[t] + V[g_in_dn,t+1]*w_in_dn + V[g_in_up,t+1]*(1-w_in_dn)
      Vout = capa_out_t*price[t] + V[g_out_dn,t+1]*w_out_dn + V[g_out_up,t+1]*(1-w_out_dn)
  
      V_t = max.(Vno, Vin, Vout)
      prod[t] = (V_t==Vin)*capa_in_t - ((V_t==Vout) & (V_t>Vin))*capa_out_t
      volume[t+1] = volume[t] + inflow[t] + prod[t]
    end
  
    # for stacking model, pumping must be negative and turbining must be positive
    # prods[splant] .= -prod
    # ind_pump = (prods[splant] .< 0)
    # prods[splant][ind_pump] .= prods[splant][ind_pump] ./ eff_pump
  
    prods .= -prod
    ind_pump = (prods .< 0)
    prods[ind_pump] .= prods[ind_pump] ./ eff_pump

    # costs[splant] .= Utils.calculate_marginal_cost(umarginal_cost, prods[splant])
    # cost = prods.*price # positive
  
    return prods #, costs
end