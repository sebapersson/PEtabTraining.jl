using Aqua, PEtabTraining

@testset "Aqua" begin
    Aqua.test_ambiguities(PEtabTraining, recursive = false)
    Aqua.test_undefined_exports(PEtabTraining)
    Aqua.test_unbound_args(PEtabTraining)
    Aqua.test_stale_deps(PEtabTraining)
    Aqua.test_deps_compat(PEtabTraining)
    Aqua.find_persistent_tasks_deps(PEtabTraining)
    Aqua.test_piracies(PEtabTraining)
    Aqua.test_project_extras(PEtabTraining)
end
