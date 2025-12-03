using ComponentArrays, Lux
import Random

function _lv_ude_function_reg!(du, u, p, t, ml_models)
    prey, predator = u
    @unpack alpha, delta = p
    net1 = ml_models[:net1]
    du_nn, st = net1.model([prey, predator], p[:net1], net1.st)
    net1.st = st

    du[1] = alpha*prey - du_nn[1] # prey
    du[2] = du_nn[2] - delta*predator # predator
    du[3] = LinearAlgebra.norm([du_nn[1], du_nn[2]])
    return nothing
end

function _get_lv_ude_model(; measurements_df::Union{DataFrame, Nothing} = nothing)
    rng = Random.MersenneTwister(42)
    lv_net = Lux.Chain(Dense(2 => 5, Lux.tanh),
        Dense(5 => 5, Lux.tanh),
        Dense(5 => 2)) |> f64
    ml_models = Dict(:net1 => MLModel(lv_net; static = false))
    lv_ude_function! = let _ml_models = ml_models
        (du, u, p, t) -> _lv_ude_function!(du, u, p, t, _ml_models)
    end

    pnn = Lux.initialparameters(rng, lv_net) |> ComponentArray |> f64
    p_mechanistic = (alpha = 1.3, delta = 1.8, beta = 0.9)
    p_ode = ComponentArray(merge(p_mechanistic, (net1 = pnn,)))
    u0 = ComponentArray(prey = 0.44249296, predator = 4.6280594, nn_norm = 0.0)
    ude_problem = ODEProblem(lv_ude_function!, u0, (0.0, 10.0), p_ode)

    p_alpha = PEtabParameter(:alpha; scale = :lin, lb = 1e-2, ub = 1e2, value = 1.3)
    p_beta = PEtabParameter(:beta; scale = :lin, lb = 1e-2, ub = 1e2, value = 0.9)
    p_delta = PEtabParameter(:delta; scale = :lin, lb = 1e-2, ub = 1e2, value = 1.8)
    p_lambda_reg = PEtabParameter(:lambda_reg; estimate = false, scale = :lin, value = 1.0)
    p_net1 = PEtabMLParameter(:net1, true, pnn)
    pest = [p_alpha, p_beta, p_delta, p_lambda_reg, p_net1]

    obs_prey = PEtabObservable(:prey, 1.0)
    obs_predator = PEtabObservable(:predator, 1.0)
    obs_reg = PEtabObservable("1e-6 * lambda_reg * (nn_norm / TMAX)^2", 1.0)
    obs = Dict("prey_o" => obs_prey, "predator_o" => obs_predator, "reg_o" => obs_reg)

    if isnothing(measurements_df)
        path_data = joinpath(@__DIR__, "lv_training_data.csv")
        training_data = CSV.read(path_data, DataFrame)
    else
        training_data = measurements_df
    end
    df_reg = DataFrame(time = [maximum(training_data.time)], measurement = [0.0], observableId = ["reg_o"])
    training_data = vcat(training_data, df_reg)
    return PEtabModel(ude_problem, obs, training_data, pest; ml_models = ml_models)
end
