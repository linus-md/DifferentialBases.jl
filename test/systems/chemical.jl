@testset "Systems -> chemical" begin
    using DifferentialBases: chem_1, akzo_nobel, fast_slow_reaction

    ideal, derivatives, R = chem_1()
    @test length(ideal.gens) == 2
    @test length(derivatives) == 9

    ideal, derivatives, R = akzo_nobel()
    @test length(ideal.gens) == 2
    @test length(derivatives) == 16

    ideal, derivatives, R = fast_slow_reaction()
    @test length(ideal.gens) == 3
    @test length(derivatives) == 12
end
