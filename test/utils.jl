@testset "is_symmetric" begin
    m = Bool[
        1 0 0
        0 1 0
        0 0 1
    ]
    @test ImageMorphology.is_symmetric(m)
    @test ImageMorphology.is_symmetric(centered(m))
    mo = strel(CartesianIndex, centered(m))
    @test ImageMorphology.is_symmetric(mo)

    m = Bool[
        1 0 1
        0 1 0
        1 1 1
    ]
    mo = strel(CartesianIndex, centered(m))
    @test !ImageMorphology.is_symmetric(m)
    @test !ImageMorphology.is_symmetric(mo)

    @test ImageMorphology.is_symmetric(strel_box((3, 3)))
    @test ImageMorphology.is_symmetric(strel_diamond((3, 3)))
end
