using Test
using Pkg
Pkg.add(url="https://github.com/algebraic-solving/AlgebraicSolving.jl")

@testset verbose = true "DifferentialBases Tests" begin
    include("algorithms/classical.jl")
end
