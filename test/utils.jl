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

@testset "require_symmetric_strel" begin
    se = rand_se_mask(3, 2; symmetric=true)
    @test_nowarn ImageMorphology.require_symmetric_strel(se)
    se = Bool[1 0 0; 0 1 0; 0 0 0]
    msg = "structuring element must be symmetric with respect to its center"
    @test_throws ArgumentError(msg) ImageMorphology.require_symmetric_strel(se)
end
