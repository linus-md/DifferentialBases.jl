@testset "Algorithms -> Classical -> intersect_ideal" begin
    using DifferentialBases: intersect_ideal
    using AlgebraicSolving: polynomial_ring, Ideal, GF, groebner_basis
    # Test intersect(G, S)
    R, (x,y,z) = polynomial_ring(GF(101),["x","y","z"], internal_ordering=:lex)
    F1 = [y*z, x*y, z*x]
    I1 = Ideal(F1)
    G1 = groebner_basis(I1)

    F2 = [x*z, x*z]
    I2 = Ideal(F2)
    G2 = groebner_basis(I2)

    S_vars = R.data.S[2:end]
    @test intersect_ideal(G1, S_vars)[1] == F1[1]
    @test intersect_ideal(G2, S_vars) == Any[]
end

@testset "Algorithms -> Classical -> partial" begin
    using DifferentialBases: partial
    using AlgebraicSolving: polynomial_ring, GF
    R, (x,y,z) = polynomial_ring(GF(101),["x","y","z"], internal_ordering=:lex)
    q = x^2 + y^2 + z^2
    derivatives_1 = Dict(x => y, y => y^2)
    derivatives_2 = Dict(x => y, y => y^2, z => z^2)
    @test partial(q, derivatives_1) == 2*x*y + 2*y^3
    @test partial(q, derivatives_2) == 2*x*y + 2*y^3 + 2*z^3
end