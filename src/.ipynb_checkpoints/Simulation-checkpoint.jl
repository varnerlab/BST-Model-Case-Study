function evaluate_model(model::Dict{String,Any}; 
    tspan::Tuple{Float64,Float64} = (0.0,20.0), Δt::Float64 = 0.01)

    # get stuff from model -
    xₒ = model["initial_condition_vector"]

    # setup the solver -
    prob = ODEProblem(balances, xₒ, tspan, model; saveat = Δt)
    soln = solve(prob)

    # get the results from the solver -
    T = soln.t
    X = soln.u

    # return -
    return (T,X)
end