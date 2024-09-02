@testset "Systems -> poly_nn" begin
    using DifferentialBases

    ideal, derivatives, R, R_vars = DifferentialBases.poly_nn_2(2, 2, 2)
    @test length(ideal.gens) == 2*2
    @test length(derivatives) == 2*2 + 2*2 + 2*2
    @test any(x -> x != 0, derivatives)
    
    ideal, derivatives, R, R_vars = DifferentialBases.poly_nn_2(2, 3, 4)
    @test length(ideal.gens) == 3*3
    @test length(derivatives) == 2*3 + 3*4 + 2 + 4
    @test any(x -> x != 0, derivatives)
end
