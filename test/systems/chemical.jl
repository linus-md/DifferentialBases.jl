@testset "Systems -> chemical" begin
    using DifferentialBases

    ideal, derivatives, R, R_vars = DifferentialBases.chem_1()
    @test length(ideal.gens) == 2
    @test length(derivatives) == 9
    
    res = DifferentialBases.differential_basis(ideal, derivatives, R, R_vars)
    @test length(res) == length(ideal.gens)

    ideal, derivatives, R, R_vars = DifferentialBases.akzo_nobel()
    @test length(ideal.gens) == 2
    @test length(derivatives) == 16

    res = DifferentialBases.differential_basis(ideal, derivatives, R, R_vars)
    @test length(res) == 5

    ideal, derivatives, R, R_vars = DifferentialBases.fast_slow_reaction()
    @test length(ideal.gens) == 3
    @test length(derivatives) == 12

    res = DifferentialBases.differential_basis(ideal, derivatives, R, R_vars)
    @test length(res) == 12
end
