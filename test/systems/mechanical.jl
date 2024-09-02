
@testset "Systems -> mechanical" begin
    using DifferentialBases

    ideal, derivatives, R, R_vars = DifferentialBases.simple_pendulum()
    @test length(ideal.gens) == 1
    @test length(derivatives) == 5

    res = DifferentialBases.differential_basis(ideal, derivatives, R, R_vars)
    @test length(res) == 8

    ideal, derivatives, R, R_vars = DifferentialBases.double_pendulum()
    @test length(ideal.gens) == 2
    @test length(derivatives) == 10

    ideal, derivatives, R, R_vars = DifferentialBases.triple_pendulum()
    @test length(ideal.gens) == 3
    @test length(derivatives) == 15

    ideal, derivatives, R, R_vars = DifferentialBases.point_mass_parabola()
    @test length(ideal.gens) == 1
    @test length(derivatives) == 7

    res = DifferentialBases.differential_basis(ideal, derivatives, R, R_vars)
    @test length(res) == 11
end
