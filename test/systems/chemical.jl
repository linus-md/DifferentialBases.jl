@testset "Systems -> chemical" begin
    using DifferentialBases: chem_1

    ideal, derivatives, R = chem_1()
    @test length(ideal.gens) == 2
    @test length(derivatives) == 4
end
