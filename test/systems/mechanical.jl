
@testset "Systems -> mechanical" begin
    using DifferentialBases: 
        simple_pendulum, double_pendulum, triple_pendulum, masspoint_parabola

    ideal, derivatives, R = simple_pendulum()
    @test length(ideal.gens) == 1
    @test length(derivatives) == 5

    ideal, derivatives, R = double_pendulum()
    @test length(ideal.gens) == 2
    @test length(derivatives) == 10

    ideal, derivatives, R = triple_pendulum()
    @test length(ideal.gens) == 3
    @test length(derivatives) == 15

    ideal, derivatives, R = masspoint_parabola()
    @test length(ideal.gens) == 1
    @test length(derivatives) == 7
end
