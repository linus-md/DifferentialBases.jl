@testset "Systems -> mechanical -> differential_basis" begin
    using DifferentialBases: simple_pendulum, masspoint_parabola

    ideal, derivatives, R = simple_pendulum()
    @test length(ideal.gens) == 1
    @test length(derivatives) == 5

    ideal, derivatives, R = masspoint_parabola()
    @test length(ideal.gens) == 1
    @test length(derivatives) == 7
end
