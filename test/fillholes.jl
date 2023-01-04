@testset "fillhole" begin
    #binary
    img = Bool[
        0 0 0 0 0 0 0
        0 1 1 1 1 1 0
        0 1 0 0 0 1 0
        0 1 0 0 0 1 0
        0 1 0 0 0 1 0
        0 1 1 1 1 1 0
        0 0 0 0 0 0 0
    ]

    expected = Bool[
        0 0 0 0 0 0 0
        0 1 1 1 1 1 0
        0 1 1 1 1 1 0
        0 1 1 1 1 1 0
        0 1 1 1 1 1 0
        0 1 1 1 1 1 0
        0 0 0 0 0 0 0
    ]

    out = fillhole(img)
    @test eltype(out) == Bool
    @test out == expected

    # in place
    out = similar(img)
    out = fillhole!(out, img)
    @test out == expected

    # in place diamond
    out = similar(img)
    out = fillhole!(out, img, strel_diamond((3, 3)))
    @test out == expected

    #gray
    img = [
        3 3 3 3 3 3 3 3 3 3
        3 4 4 4 3 3 4 4 4 3
        3 4 1 4 3 3 4 1 4 3
        3 4 4 4 3 3 4 4 4 3
        3 3 3 3 3 3 3 3 3 3
    ]

    expected = [
        3 3 3 3 3 3 3 3 3 3
        3 4 4 4 3 3 4 4 4 3
        3 4 4 4 3 3 4 4 4 3
        3 4 4 4 3 3 4 4 4 3
        3 3 3 3 3 3 3 3 3 3
    ]

    out = fillhole(img)
    @test out == expected

    msg = "the input structuring element is not for 1 dimensional array, instead it is for 2 dimensional array"
    @test_throws DimensionMismatch(msg) fillhole(rand(10), strel_box((3, 3)))

    se = strel_diamond((7, 7))
    msg = "structuring element with half-size larger than 1 is invalid"
    @test_throws DimensionMismatch(msg) fillhole(rand(10, 10), se)

end