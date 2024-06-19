using Test
using Pkg

@testset verbose = true "DifferentialBases Tests" begin
    include("algorithms/classical.jl")
end
