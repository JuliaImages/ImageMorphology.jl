@testset "connectivity_ball" begin
    _build = ImageMorphology.connectivity_ball

    # N=1
    out = @inferred _build(1, Val(1))
    @test typeof(out) == BitVector
    @test out == [true, true, true]

    # N=2 (default case)
    out = @inferred _build(1, Val(2))
    @test typeof(out) == BitMatrix
    @test out == trues(3, 3)
    @test out == _build(1)

    out = @inferred _build((1, 3), Val(2)) # ellipse
    ref = Bool[0 0 0 0 0 0 0
        0 0 0 0 0 0 0
        0 0 0 1 0 0 0
        1 1 1 1 1 1 1
        0 0 0 1 0 0 0
        0 0 0 0 0 0 0
        0 0 0 0 0 0 0]
    @test typeof(out) == BitMatrix
    @test out == ref

    out = @inferred _build((1, 3), Val(2); square=false) # ellipse
    ref = Bool[0 0 0 1 0 0 0; 1 1 1 1 1 1 1; 0 0 0 1 0 0 0]
    @test typeof(out) == BitMatrix
    @test out == ref
    # ensure keywords are passed through
    @test ref == _build((1, 3); square=false) == _build((1, 3), 2; square=false)

    # N=3
    out = if VERSION >= v"1.6"
        @inferred _build(2, Val(3))
    else
        _build(2, Val(3))
    end
    ref = cat(
        Bool[0 0 0 0 0; 0 0 0 0 0; 0 0 1 0 0; 0 0 0 0 0; 0 0 0 0 0],
        Bool[0 0 0 0 0; 0 1 1 1 0; 0 1 1 1 0; 0 1 1 1 0; 0 0 0 0 0],
        Bool[0 0 1 0 0; 0 1 1 1 0; 1 1 1 1 1; 0 1 1 1 0; 0 0 1 0 0],
        Bool[0 0 0 0 0; 0 1 1 1 0; 0 1 1 1 0; 0 1 1 1 0; 0 0 0 0 0],
        Bool[0 0 0 0 0; 0 0 0 0 0; 0 0 1 0 0; 0 0 0 0 0; 0 0 0 0 0];
        dims=3
    )
    @test typeof(out) == BitArray{3}
    @test out == ref

    # not type-inferable
    err = ErrorException("return type $(BitMatrix) does not match inferred return type Any")
    @test_throws err @inferred _build(1, 2)
    @test _build(1, 2) == _build(1, Val(2))

    # exceptions
    @test_throws MethodError _build((1,), Val(2))
    @test_throws MethodError _build((1, 2, 3), Val(2))
end

@testset "connectivity_region" begin
    _build = ImageMorphology.connectivity_region

    # N=2, C4
    out = @inferred _build(:C4)
    ref = Bool[false true false; true true true; false true false]
    @test typeof(out) == BitMatrix
    @test out == ref
    @test out == _build(:C4, Val(2))

    # N=2, C8
    out = _build(:C8)
    ref = trues(3, 3)
    @test typeof(out) == BitMatrix
    @test out == ref

    # not guaranteed to be type-inferable
    @test _build(:C4, 2) == _build(:C4, Val(2))

    # exceptions
    err = ErrorException("unknown connectivity alias for dimension 2: null")
    @test_throws err _build(:null)
    err = ErrorException("unsupported dimension 999")
    @test_throws err _build(:null, Val(999))
end
