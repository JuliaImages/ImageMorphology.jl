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

@testset "require_select_function" begin
    @test_nowarn ImageMorphology.require_select_function(max, Int)
    @test_nowarn ImageMorphology.require_select_function(min, Gray{N0f8})

    msg = "function `min` is not a well-defined select function on type `ComplexF64` and `ComplexF64`: does `f(x::T1, y::T2)` work as expected?"
    @test_throws ArgumentError(msg) ImageMorphology.require_select_function(min, ComplexF64)
    f = (x, y) -> min(x, y)
    msg = "function `$f` is not a well-defined select function on type `RGB{Float32}` and `RGB{Float32}`: does `f(x::T1, y::T2)` work as expected?"
    @test_throws ArgumentError(msg) ImageMorphology.require_select_function(f, RGB{Float32})
end

@testset "staturated" begin
    # unsigned
    @test ImageMorphology.saturating_add(N0f8(0.6), N0f8(0.6)) == N0f8(1.0)
    @test ImageMorphology.saturating_sub(N0f8(0.2), N0f8(0.6)) == N0f8(0.0)
end