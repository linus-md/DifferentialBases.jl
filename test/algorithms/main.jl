@testset "Algorithms -> Main -> _intersect" begin
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
    @test DifferentialBases._intersect(G1, S_vars)[1] == F1[1]

    F2 = [x*z, x*z]
    I2 = AlgebraicSolving.Ideal(F2)
    G2 = AlgebraicSolving.groebner_basis(I2)
    @test DifferentialBases._intersect(G2, S_vars) == Any[]
end

@testset "Algorithms -> Main -> _diff_op " begin
    using DifferentialBases
    using AlgebraicSolving

    R, (x,y,z) = AlgebraicSolving.polynomial_ring(
        AlgebraicSolving.GF(101),
        ["x","y","z"],
        internal_ordering=:lex)

    q = x^2 + y^2 + z^2
    derivatives_1 = Dict(x => y, y => y^2)
    @test DifferentialBases._diff_op(q, derivatives_1) == 2*x*y + 2*y^3

    derivatives_2 = Dict(x => y, y => y^2, z => z^2)
    @test DifferentialBases._diff_op(q, derivatives_2) == 2*x*y + 2*y^3 + 2*z^3
end

@testset "Algorithms -> Main -> differential_basis" begin
    using DifferentialBases
    using AlgebraicSolving

    R, R_vars = AlgebraicSolving.polynomial_ring(
        AlgebraicSolving.GF(101),
        ["dl","l","v","u","y","x"], 
        internal_ordering=:degrevlex)

    (dl, l, v, u, y, x) = R_vars

    derivatives = Dict(
        x => u,
        y => v,
        u => x*l,
        v => y*l - 1,
        l => dl
    )

    ideal = AlgebraicSolving.Ideal([x^2 + y^2 - 1])

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

    @test DifferentialBases.differential_basis(
        ideal, derivatives, R, R_vars) == sol
    @test DifferentialBases.differential_basis(
        AlgebraicSolving.Ideal(sol), derivatives, R, R_vars, true, 1) == sol
end

@testset "Algorithms -> Main -> Ring helpers" begin
    using DifferentialBases
    using AlgebraicSolving

    R, R_vars = AlgebraicSolving.polynomial_ring(
        AlgebraicSolving.GF(101),
        ["l","v","u","y","x","dl"], 
        internal_ordering=:degrevlex)

    (l, v, u, y, x, dl) = R_vars

    derivatives = Dict(
        x => u,
        y => v,
        u => x*l,
        v => y*l - 1,
        l => dl
    )

    ideal = AlgebraicSolving.Ideal([x^2 + y^2 - 1])

    R_new, R_new_vars, map_old_new = DifferentialBases._manage_rings(
        derivatives, R)

    R_1, R_1_vars = AlgebraicSolving.polynomial_ring(
        AlgebraicSolving.GF(101),
        ["dl","l","v","u","y","x",], 
        internal_ordering=:degrevlex)
    
    @test R_new == R_1

    ideal_new_gens = [DifferentialBases._swap_vars(elem, R_vars, R_new_vars, map_old_new) for elem in ideal.gens]
    ideal_new = AlgebraicSolving.Ideal(ideal_new_gens)

    @test parent(ideal_new[1]) == R_new


    # Test managing for QQ/GF
    res = DifferentialBases.differential_basis(
        ideal, derivatives, R, R_vars, true, 2)
    @test length(res) == 8

    ideal, derivatives, R, R_vars = DifferentialBases.simple_pendulum()
    res2 = DifferentialBases.differential_basis(
        ideal, derivatives, R, R_vars, false, 2)
    @test length(res2) == 7

    ideal, derivatives, R, R_vars = DifferentialBases.simple_pendulum()
    res2 = DifferentialBases.differential_basis(
        ideal, derivatives, R, R_vars, true, 2)
    @test typeof(res2) == Nothing
end