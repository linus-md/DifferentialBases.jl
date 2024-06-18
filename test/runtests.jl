using DifferentialBases: intersect_ideal, partial
using Test

@testset verbose = true "DifferentialBases Tests" begin
    include("algorithms/classical.jl")
end