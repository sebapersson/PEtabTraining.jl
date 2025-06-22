function test_nllh(path_yaml, mdf::DataFrame, mdf_tmp::DataFrame,
        petab_prob, stage_problems_i)::Nothing
    CSV.write(petab_prob.model_info.model.paths[:measurements], mdf_tmp, delim = '\t')
    petab_prob_ref = PEtabModel(path_yaml) |> PEtabODEProblem
    CSV.write(petab_prob.model_info.model.paths[:measurements], mdf, delim = '\t')
    nllh_ref = petab_prob_ref.nllh(get_x(petab_prob_ref))
    nllh_test = stage_problems_i.nllh(get_x(petab_prob))
    @test nllh_ref == nllh_test
    return nothing
end

function test_chunk_size(chunks)
    n_stages = length(chunks)
    for i in 1:(length(chunks) - 1)
        # Sometimes chunks cannot be divided uniformly, then the last is bigger
        if i == length(chunks) - 1
            @test abs(length(chunks[i]) - length(chunks[i + 1])) < n_stages
        else
            @test length(chunks[i]) == length(chunks[i + 1])
        end
    end
end
