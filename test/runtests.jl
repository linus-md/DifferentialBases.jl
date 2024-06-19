using Test
using Pkg
Pkg.add(PackageSpec(url="https://github.com/linus-md/AlgebraicSolving.jl/", rev="block"))

@testset verbose = true "DifferentialBases Tests" begin
    include("algorithms/classical.jl")
end
