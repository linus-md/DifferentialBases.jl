module DifferentialBases

include("algorithms/classical.jl")
include("systems/mechanical.jl")
include("systems/linear_nn.jl")
include("systems/chemical.jl")
end

using AlgebraicSolving: polynomial_ring, GF, Ideal, groebner_basis, normal_form
i, d, R = DifferentialBases.double_pendulum_fixed_l1l2()
G = DifferentialBases.deg_stop_differential_basis(i, d, 9, 2, false)
G2 = groebner_basis(Ideal(G))
println("The size of G is ", length(G))
println("The size of G2 is ", length(G2))