@testset "low_leveling" begin
    ref = ([
        2 2 2
        5 5 5
        1 1 1])

    marker = ([
        1 1 1
        3 3 3
        4 4 4])

    expected = ([
        2 2 2
        4 4 4
        1 1 1])
    out = low_leveling(ref, marker)
    @test out == expected

    ref = ([
        5 0 0 0 0 0 5
        0 5 0 0 0 5 0
        0 0 5 5 5 0 0
        0 0 5 5 5 0 0
        0 0 5 5 5 0 0
        0 5 0 0 0 5 0
        5 0 0 0 0 0 5])

    marker = ([
        0 0 0 0 0 0 0
        0 0 0 0 0 0 0
        0 0 6 6 6 0 0
        0 0 6 6 6 0 0
        0 0 6 6 6 0 0
        0 0 0 0 0 0 0
        0 0 0 0 0 0 0])

    expected = ([
        5 0 0 0 0 0 5
        0 5 0 0 0 5 0
        0 0 5 5 5 0 0
        0 0 5 5 5 0 0
        0 0 5 5 5 0 0
        0 5 0 0 0 5 0
        5 0 0 0 0 0 5])
    out = low_leveling(ref, marker)
    @test out == expected

    out = similar(ref)
    out = low_leveling!(out, ref, marker)
    @test out == expected

    msg = "the input structuring element is not for 1 dimensional array, instead it is for 2 dimensional array"
    @test_throws DimensionMismatch(msg) low_leveling(rand(10), rand(10), strel_box((3, 3)))

    se = strel_diamond((7, 7))
    msg = "structuring element with half-size larger than 1 is invalid"
    @test_throws DimensionMismatch(msg) low_leveling(rand(10, 10), rand(10, 10), se)

    se = strel_diamond((3, 3))
    msg = "images should have the same axes"
    @test_throws DimensionMismatch(msg) low_leveling(rand(10, 10), rand(5, 5), se)
end

@testset "high_leveling" begin
    ref = ([
        2 2 2
        5 5 5
        1 1 1])

    marker = ([
        1 1 1
        3 3 3
        4 4 4])

    expected = ([
        2 2 2
        5 5 5
        3 3 3])

    out = high_leveling(ref, marker)
    @test out == expected

    ref = ([
        5 0 0 0 0 0 5
        0 5 0 0 0 5 0
        0 0 5 5 5 0 0
        0 0 5 5 5 0 0
        0 0 5 5 5 0 0
        0 5 0 0 0 5 0
        5 0 0 0 0 0 5])

    marker = ([
        0 0 0 0 0 0 0
        0 0 0 0 0 0 0
        0 0 4 4 4 0 0
        0 0 4 4 4 0 0
        0 0 4 4 4 0 0
        0 0 0 0 0 0 0
        0 0 0 0 0 0 0])

    expected = ([
        5 0 0 0 0 0 5
        0 5 0 0 0 5 0
        0 0 5 5 5 0 0
        0 0 5 5 5 0 0
        0 0 5 5 5 0 0
        0 5 0 0 0 5 0
        5 0 0 0 0 0 5])

    out = high_leveling(ref, marker)
    @test out == expected

    out = similar(marker)
    out = high_leveling!(out, ref, marker)
    @test out == expected

    msg = "the input structuring element is not for 1 dimensional array, instead it is for 2 dimensional array"
    @test_throws DimensionMismatch(msg) high_leveling(rand(10), rand(10), strel_box((3, 3)))

    se = strel_diamond((7, 7))
    msg = "structuring element with half-size larger than 1 is invalid"
    @test_throws DimensionMismatch(msg) high_leveling(rand(10, 10), rand(10, 10), se)

    se = strel_diamond((3, 3))
    msg = "images should have the same axes"
    @test_throws DimensionMismatch(msg) high_leveling(rand(10, 10), rand(5, 5), se)
end

@testset "leveling" begin

    ref = ([1 1 2 1 1 4 4 4 4 4 0 0 0 0 0])
    marker = ([0 0 0 0 0 0 6 6 6 6 6 0 0 2 0])
    expected = ([1 1 1 1 1 4 4 4 4 4 0 0 0 0 0])
    out = leveling(ref, marker)
    @test out == expected

    ref = ([1 1 1 1 5 5 5 5])
    marker = ([4 4 4 4 2 2 2 2])
    expected = ([2 2 2 2 4 4 4 4])
    out = leveling(ref, marker)
    @test out == expected

    out = similar(marker)
    out = leveling!(out, ref, marker)
    @test out == expected
end