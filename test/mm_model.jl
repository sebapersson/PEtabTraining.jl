using Catalyst, PEtab, OrdinaryDiffEqRosenbrock, DataFrames

function _get_mm_model(; measurements_df::Union{DataFrame, Nothing} = nothing)
    rn = @reaction_network begin
        @parameters S0 c3=1.0
        @species S(t)=S0
        c1, S + E --> SE
        c2, SE --> S + E
        c3, SE --> P + E
    end
    speciemap = [:E => 50.0, :SE => 0.0, :P => 0.0]

    @unpack E, S, P = rn
    @parameters sigma
    observables = [
        PEtabObservable(:obs_sum, S + E, 1.0),
        PEtabObservable(:obs_p, P, sigma)
    ]

    p_c1 = PEtabParameter(:c1; value = 1.0)
    p_c2 = PEtabParameter(:c2; value = 10.0)
    p_s0 = PEtabParameter(:S0; value = 100.0)
    p_sigma = PEtabParameter(:sigma, value = 1.0, scale = :lin)
    pest = [p_c1, p_c2, p_s0, p_sigma]

    if isnothing(measurements_df)
        ps = [:c1 => 1.0, :c2 => 10.0, :c3 => 1.0, :S0 => 100.0]
        u0 = [:S => 100.0, :E => 50.0, :SE => 0.0, :P => 0.0]
        tspan = (0.0, 10.0)
        oprob = ODEProblem(rn, u0, tspan, ps)
        sol = solve(oprob, Rodas5P(); saveat = 0:0.5:10.0)
        obs_sum = (sol[:S] + sol[:E]) .+ randn(length(sol[:E]))
        obs_p = sol[:P] + .+ randn(length(sol[:P]))
        df_sum = DataFrame(obs_id = "obs_sum", time = sol.t, measurement = obs_sum)
        df_p = DataFrame(obs_id = "obs_p", time = sol.t, measurement = obs_p)
        measurements = vcat(df_sum, df_p)
    else
        measurements = measurements_df
    end

    return PEtabModel(rn, observables, measurements, pest; speciemap = speciemap)
end
