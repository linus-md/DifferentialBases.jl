using Test

@testset verbose = true "DifferentialBases Tests" begin
    include("algorithms/classical.jl")
    include("systems/mechanical.jl")
    include("systems/chemical.jl")
end
