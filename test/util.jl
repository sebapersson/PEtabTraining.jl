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

# Allocate Curriculum stages
cl_epochs = allocate_cl_epochs(6000, 5, 1 / 3)
@test cl_epochs == [
    1 => 1:500, 2 => 501:1000, 3 => 1001:1500, 4 => 1501:2000, 5 => 2001:6000,
]
cl_epochs = allocate_cl_epochs(6000, 5, 0.5)
@test cl_epochs == [
    1 => 1:750, 2 => 751:1500, 3 => 1501:2250, 4 => 2251:3000, 5 => 3001:6000,
]
# Test epochs are evenly distributed across stages
for n_stages in 3:15
    cl_epochs = allocate_cl_epochs(6000, n_stages, 0.47)
    diff_cl = [length(cl_epochs[i + 1].second) - length(cl_epochs[i].second) for i in 1:(n_stages - 2)]
    @test all(abs.(diff_cl) .≤ 1)
end
# Error checking
@test_throws ArgumentError allocate_cl_epochs(6000, 5, 1.2)
@test_throws ArgumentError allocate_cl_epochs(6000, 5, 0.0)
@test_throws ArgumentError allocate_cl_epochs(10, 10, 0.3)

# Test that the splitting functions split the data into chunks of approximately equal size
for n_chunks in [3, 5, 13, 21]
    c = PEtabTraining._makechunks(collect(1:61), n_chunks)
    diff = [length(c[i + 1]) - length(c[i]) for i in 1:(length(c) - 1)]
    @test all(abs.(diff) .≤ 1)
end
