
@testset "Systems -> mechanical" begin
    using DifferentialBases

    ideal, derivatives, R = DifferentialBases.simple_pendulum()
    @test length(ideal.gens) == 1
    @test length(derivatives) == 5

    ideal, derivatives, R = DifferentialBases.double_pendulum()
    @test length(ideal.gens) == 2
    @test length(derivatives) == 10

    ideal, derivatives, R = DifferentialBases.triple_pendulum()
    @test length(ideal.gens) == 3
    @test length(derivatives) == 15

    ideal, derivatives, R = DifferentialBases.masspoint_parabola()
    @test length(ideal.gens) == 1
    @test length(derivatives) == 7
end
