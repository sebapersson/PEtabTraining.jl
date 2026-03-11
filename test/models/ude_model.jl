using ComponentArrays, Lux
import LinearAlgebra, Random

function _lv_ude_function_reg!(du, u, p, t, ml_models)
    prey, predator = u
    @unpack alpha, delta = p
    net1 = ml_models[:net1]
    du_nn, st = net1.lux_model([prey, predator], p[:net1], net1.st)
    net1.st = st

    du[1] = alpha * prey - du_nn[1] # prey
    du[2] = du_nn[2] - delta * predator # predator
    du[3] = LinearAlgebra.norm([du_nn[1], du_nn[2]])
    return nothing
end

function _get_lv_ude_model(;
        measurements_df::Union{DataFrame, Nothing} = nothing,
        include_regularization::Bool = false
    )
    rng = Random.MersenneTwister(42)
    lv_net = Lux.Chain(
        Dense(2 => 5, Lux.tanh),
        Dense(5 => 5, Lux.tanh),
        Dense(5 => 2)
    ) |> f64
    ml_model = PEtab.MLModel(:net1, lv_net, false)

    pnn = Lux.initialparameters(rng, lv_net) |> ComponentArray |> f64
    p_mechanistic = (alpha = 1.3, delta = 1.8, beta = 0.9)
    u0 = ComponentArray(prey = 0.44249296, predator = 4.6280594, nn_norm = 0.0)
    ude_problem = UDEProblem(_lv_ude_function_reg!, u0, (0.0, 10.0), p_mechanistic, ml_model)

    pest = [
        PEtabParameter(:alpha; scale = :lin, lb = 1.0e-2, ub = 1.0e2, value = 1.3),
        PEtabParameter(:beta; scale = :lin, lb = 1.0e-2, ub = 1.0e2, value = 0.9),
        PEtabParameter(:delta; scale = :lin, lb = 1.0e-2, ub = 1.0e2, value = 1.8),
        PEtabParameter(:lambda_reg; estimate = false, scale = :lin, value = 1.0),
        PEtabMLParameter(:net1; value = pnn),
    ]

    observables = [
        PEtabObservable(:prey_o, :prey, 1.0),
        PEtabObservable(:predator_o, :predator, 1.0),
        PEtabObservable(:reg_o, "1e-6 * lambda_reg * (nn_norm)^2", 1.0),
    ]

    if isnothing(measurements_df)
        path_data = joinpath(@__DIR__, "ude_model_data.csv")
        training_data = CSV.read(path_data, DataFrame)
    else
        training_data = measurements_df
    end
    if include_regularization
        df_reg = DataFrame(
            time = [maximum(training_data.time)], measurement = [0.0],
            observableId = ["reg_o"]
        )
        training_data = vcat(training_data, df_reg)
    end
    return PEtabModel(
        ude_problem, observables, training_data, pest; ml_models = ml_model
    )
end
