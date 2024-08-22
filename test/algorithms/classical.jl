@testset "Algorithms -> Classical -> intersect" begin
    using DifferentialBases
    using AlgebraicSolving
    R, (x,y,z) = AlgebraicSolving.polynomial_ring(
        AlgebraicSolving.GF(101),
        ["x","y","z"],
        internal_ordering=:lex)
        
    S_vars = R.data.S[2:end]
    F1 = [y*z, x*y, z*x]
    I1 = AlgebraicSolving.Ideal(F1)
    G1 = AlgebraicSolving.groebner_basis(I1)
    @test DifferentialBases.intersect(G1, S_vars)[1] == F1[1]

    F2 = [x*z, x*z]
    I2 = AlgebraicSolving.Ideal(F2)
    G2 = AlgebraicSolving.groebner_basis(I2)
    @test DifferentialBases.intersect(G2, S_vars) == Any[]
end

@testset "Algorithms -> Classical -> delta" begin
    using DifferentialBases
    using AlgebraicSolving
    R, (x,y,z) = AlgebraicSolving.polynomial_ring(
        AlgebraicSolving.GF(101),
        ["x","y","z"],
        internal_ordering=:lex)

    q = x^2 + y^2 + z^2
    derivatives_1 = Dict(x => y, y => y^2)
    @test DifferentialBases.delta(q, derivatives_1) == 2*x*y + 2*y^3

    derivatives_2 = Dict(x => y, y => y^2, z => z^2)
    @test DifferentialBases.delta(q, derivatives_2) == 2*x*y + 2*y^3 + 2*z^3
end

@testset "Algorithms -> Classical -> differential_basis" begin
    using DifferentialBases
    using AlgebraicSolving

    R, (dl,l,v,u,y,x) = AlgebraicSolving.polynomial_ring(
        AlgebraicSolving.GF(101),
        ["dl","l","v","u","y","x"], 
        internal_ordering=:degrevlex)

    derivatives = Dict(
        x => u,
        y => v,
        u => x*l,
        v => y*l - 1,
        l => dl
    )

    sol = [
        y^2 + x^2 + 100,
        v*y + u*x,
        l + v^2 + u^2 + 100*y,
        l*x^2 + 100*l + 100*u^2 + 100*y*x^2 + y,
        100*v*x^2 + v + u*y*x,
        l*y*x + 100*v*u + x^3 + 100*x,
        l*y + 100*v*u*x + u^2*y + x^2 + 100,
        dl + 98*v
    ]

    ideal = AlgebraicSolving.Ideal([x^2 + y^2 - 1])
    @test DifferentialBases.differential_basis(ideal, derivatives, R) == sol
    @test DifferentialBases.differential_basis(
        AlgebraicSolving.Ideal(sol), derivatives, R, true, 1) == sol
end
