@testset "Algorithms -> Gröbner bases" begin
    using AbstractAlgebra: vars
    using AlgebraicSolving: polynomial_ring, Ideal, GF, groebner_basis
    
    # Test intersect(G, S)
    R, (x,y,z) = polynomial_ring(GF(101),["x","y","z"], internal_ordering=:lex)
    F = [y*z, x*y, z*x]
    I = Ideal(F)
    G = groebner_basis(I)

    S, (y, z) = polynomial_ring(GF(101), ["y","z"], internal_ordering=:lex)
    @test intersect(G, S)[1] == F[1]
end