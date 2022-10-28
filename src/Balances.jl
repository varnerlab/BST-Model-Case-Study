function balances(dx, x, p, t)

    # grab data from the p-vector
    α = p[1]
    G = p[2]
    S = p[3]
    number_of_dynamic_states = p[4]
    static_factors_array = p[5]
        
    # build the "state" array (dynamic | static)
    state_array = vcat(x,static_factors_array)

    # compute the kinetics - powerlaw
    rV = powerlaw(state_array,α,G)

    # compute the rhs -> store in a temp vector
    tmp = S*rV

    # populate the dx vector with the tmp vector -
    for i ∈ 1:number_of_dynamic_states
        dx[i] = tmp[i]
    end
end