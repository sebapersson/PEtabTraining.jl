function __PEtabModel(paths::Dict{Symbol, String}, petab_tables::Dict{Symbol, DataFrame};
                      build_julia_files::Bool = true, verbose::Bool = false,
                      ifelse_to_callback::Bool = true, write_to_file::Bool = false)::PEtabModel
    name = splitdir(paths[:dirmodel])[end]

    write_to_file && !isdir(paths[:dirjulia]) && mkdir(paths[:dirjulia])
    if write_to_file == false && build_julia_files == true
        paths[:dirjulia] = ""
    end

    # Import SBML model with SBMLImporter
    # In case one of the conditions in the PEtab table assigns an initial specie value,
    # the SBML model must be mutated to add an iniitial value parameter to correctly
    # compute gradients
    model_SBML = SBMLImporter.parse_SBML(paths[:SBML], false; model_as_string = false,
                                         ifelse_to_callback = ifelse_to_callback,
                                         inline_assignment_rules = false)
    PEtab._addu0_parameters!(model_SBML, petab_tables[:conditions], petab_tables[:parameters])
    pathmodel = joinpath(paths[:dirjulia], name * ".jl")
    exist = isfile(pathmodel)
    if !exist || build_julia_files == true
        btime = @elapsed begin
            model_SBML_sys = SBMLImporter._to_system_syntax(model_SBML, false, false)
            modelstr = SBMLImporter.write_reactionsystem(model_SBML_sys, paths[:dirjulia],
                                                         model_SBML;
                                                         write_to_file = write_to_file)
        end
    else
        modelstr = PEtab._get_functions_as_str(pathmodel, 1)[1]
    end

    btime = @elapsed begin
        get_rn = @RuntimeGeneratedFunction(Meta.parse(modelstr))
        # Argument needed by @RuntimeGeneratedFunction
        rn, speciemap, parametermap = get_rn("https://xkcd.com/303/")
        _odesystem = convert(ODESystem, Catalyst.complete(rn))
        # DAE requires special processing
        if isempty(model_SBML.algebraic_rules)
            odesystem = structural_simplify(_odesystem)
        else
            odesystem = structural_simplify(dae_index_lowering(_odesystem))
        end
    end
    # The state-map is not in the same order as unknowns(system) so the former is reorded
    # to make it easier to build the u0 function
    speciemap = PEtab._reorder_speciemap(speciemap, odesystem)

    # Indices for mapping parameters and tracking which parameter to estimate, useful
    # when building the comig PEtab functions
    xindices = PEtab.ParameterIndices(petab_tables, odesystem, parametermap, speciemap)

    path_u0_h_σ = joinpath(paths[:dirjulia], "$(name)_h_sd_u0.jl")
    exist = isfile(path_u0_h_σ)
    if !exist || build_julia_files == true
        btime = @elapsed begin
            hstr, u0!str, u0str, σstr = PEtab.parse_observables(name, paths, odesystem, petab_tables[:observables], xindices, speciemap, model_SBML, write_to_file)
        end
    else
        hstr, u0!str, u0str, σstr = PEtab._get_functions_as_str(path_u0_h_σ, 4)
    end
    compute_h = @RuntimeGeneratedFunction(Meta.parse(hstr))
    compute_u0! = @RuntimeGeneratedFunction(Meta.parse(u0!str))
    compute_u0 = @RuntimeGeneratedFunction(Meta.parse(u0str))
    compute_σ = @RuntimeGeneratedFunction(Meta.parse(σstr))

    # SBMLImporter holds the callback building functionality. However, currently it needs
    # to know the ODESystem parameter order (psys), and whether or not any parameters
    # which are estimated are present in the event condition. For the latter, timespan
    # should not be converted to floats in case dual numbers (for gradients) are propegated
    btime = @elapsed begin
        float_tspan = PEtab._xdynamic_in_event_cond(model_SBML, xindices, petab_tables) |> !
        psys = PEtab._get_sys_parameters(odesystem, speciemap, parametermap) .|> string
        cbset = SBMLImporter.create_callbacks(odesystem, model_SBML, name;
                                              p_PEtab = psys, float_tspan = float_tspan)
    end

    return PEtabModel(name, compute_h, compute_u0!, compute_u0, compute_σ, float_tspan,
                      paths, odesystem, deepcopy(odesystem), parametermap, speciemap,
                      petab_tables, cbset, false)
end
