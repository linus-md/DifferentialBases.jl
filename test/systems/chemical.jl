@testset "Systems -> chemical" begin
    using DifferentialBases: chem_1, akzo_nobel

    ideal, derivatives, R = chem_1()
    @test length(ideal.gens) == 2
    @test length(derivatives) == 4

    ideal, derivatives, R = akzo_nobel()
    @test length(ideal.gens) == 2
    @test length(derivatives) == 15
end
