using ComponentArrays, PEtabTraining, Test

x = ComponentVector(
    alpha = 1.3,
    beta = 0.9,
    delta = 1.8,
    net1 = (
        layer1 = (weight = randn(5, 2), bias = randn(5)),
        layer2 = (weight = randn(4, 6), bias = randn(4)),
        layer3 = (weight = randn(3, 8), bias = randn(3)),
    )
)
y = ComponentVector(
    beta = x.beta,
    net1 = (
        layer3 = x.net1.layer3,
        layer1 = x.net1.layer1,
        layer2 = x.net1.layer2,
    ),
    delta = x.delta,
    alpha = x.alpha
)

@testset "Permutation from label" begin
    ix = PEtabTraining._perm_from_labels(x, y)
    @test all(x .== y[ix])
    iy = PEtabTraining._perm_from_labels(y, x)
    @test all(x[iy] .== y)
end
