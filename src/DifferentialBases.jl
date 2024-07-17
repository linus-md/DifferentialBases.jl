module DifferentialBases

include("algorithms/classical.jl")
include("systems/mechanical.jl")
include("systems/linear_nn.jl")
include("systems/chemical.jl")
end

using AlgebraicSolving: polynomial_ring, GF, Ideal
i, d, R = DifferentialBases.double_pendulum_fixed_l1l2()
G = DifferentialBases.deg_stop_differential_basis(i, d, R, 2, false, 2)
