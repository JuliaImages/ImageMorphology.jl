@testset "hmaxima" begin
    img = UInt8[
        3 3 3 3 3 3 3 3 3 3
        3 4 4 4 3 3 4 4 4 3
        3 4 5 4 3 3 4 9 4 3
        3 4 4 4 3 3 4 4 4 3
        3 3 3 3 3 3 3 3 3 3
    ]

    ref_img = UInt8[
        3 3 3 3 3 3 3 3 3 3
        3 3 3 3 3 3 4 4 4 3
        3 3 3 3 3 3 4 6 4 3
        3 3 3 3 3 3 4 4 4 3
        3 3 3 3 3 3 3 3 3 3
    ]

    out = hmaxima(reinterpret(N0f8, img), reinterpret(N0f8, UInt8(3)))
    @test eltype(out) == N0f8
    @test out == reinterpret(N0f8, ref_img)
end

@testset "hminima" begin
    img = UInt8[
        8 8 8 8 8 8 8 8 8 8
        8 7 7 7 8 8 7 7 7 8
        8 7 6 7 8 8 7 3 7 8
        8 7 7 7 8 8 7 7 7 8
        8 8 8 8 8 8 8 8 8 8
    ]

    ref_img = UInt8[
        8 8 8 8 8 8 8 8 8 8
        8 8 8 8 8 8 7 7 7 8
        8 8 8 8 8 8 7 6 7 8
        8 8 8 8 8 8 7 7 7 8
        8 8 8 8 8 8 8 8 8 8
    ]

    out = hminima(reinterpret(N0f8, img), reinterpret(N0f8, UInt8(3)))
    @test eltype(out) == N0f8
    @test out == reinterpret(N0f8, ref_img)
end

@testset "hmaxima!" begin
    img = UInt8[
        3 3 3 3 3 3 3 3 3 3
        3 4 4 4 3 3 4 4 4 3
        3 4 5 4 3 3 4 9 4 3
        3 4 4 4 3 3 4 4 4 3
        3 3 3 3 3 3 3 3 3 3
    ]

    ref_img = UInt8[
        3 3 3 3 3 3 3 3 3 3
        3 3 3 3 3 3 4 4 4 3
        3 3 3 3 3 3 4 6 4 3
        3 3 3 3 3 3 4 4 4 3
        3 3 3 3 3 3 3 3 3 3
    ]

    out = reinterpret(N0f8, similar(img))
    out = hmaxima!(out, reinterpret(N0f8, img), reinterpret(N0f8, UInt8(3)))
    @test out == reinterpret(N0f8, ref_img)
end

@testset "hminima!" begin
    img = UInt8[
        8 8 8 8 8 8 8 8 8 8
        8 7 7 7 8 8 7 7 7 8
        8 7 6 7 8 8 7 3 7 8
        8 7 7 7 8 8 7 7 7 8
        8 8 8 8 8 8 8 8 8 8
    ]

    ref_img = UInt8[
        8 8 8 8 8 8 8 8 8 8
        8 8 8 8 8 8 7 7 7 8
        8 8 8 8 8 8 7 6 7 8
        8 8 8 8 8 8 7 7 7 8
        8 8 8 8 8 8 8 8 8 8
    ]

    out = reinterpret(N0f8, similar(img))
    out = hminima!(out, reinterpret(N0f8, img), reinterpret(N0f8, UInt8(3)))
    @test out == reinterpret(N0f8, ref_img)
end

@testset "regional_maxima" begin
    img = UInt8[
        0 0 0 0 0 0 0
        0 1 1 1 1 1 0
        0 1 2 2 2 1 0
        0 1 2 3 2 1 0
        0 1 2 2 2 1 0
        0 1 1 1 1 1 0
        0 0 0 0 0 0 0]

    ref_img = Bool[
        0 0 0 0 0 0 0
        0 0 0 0 0 0 0
        0 0 0 0 0 0 0
        0 0 0 1 0 0 0
        0 0 0 0 0 0 0
        0 0 0 0 0 0 0
        0 0 0 0 0 0 0]

    out = regional_maxima(reinterpret(N0f8, img))
    @test eltype(out) == Bool
    @test out == ref_img

    out = regional_maxima(reinterpret(N0f8, img), strel_diamond((3, 3)))
    @test eltype(out) == Bool
    @test out == ref_img
end

@testset "regional_minima" begin
    img = UInt8[
        0 0 0 0 0 0 0
        0 3 3 3 3 3 0
        0 3 2 2 2 3 0
        0 3 2 1 2 3 0
        0 3 2 2 2 3 0
        0 3 3 3 3 3 0
        0 0 0 0 0 0 0]


    ref_img = Bool[
        1 1 1 1 1 1 1
        1 0 0 0 0 0 1
        1 0 0 0 0 0 1
        1 0 0 1 0 0 1
        1 0 0 0 0 0 1
        1 0 0 0 0 0 1
        1 1 1 1 1 1 1]

    out = regional_minima(reinterpret(N0f8, img))
    @test eltype(out) == Bool
    @test out == ref_img

    out = regional_minima(reinterpret(N0f8, img), strel_diamond((3, 3)))
    @test eltype(out) == Bool
    @test out == ref_img
end

@testset "regional_maxima!" begin
    img = UInt8[
        0 0 0 0 0 0 0
        0 1 1 1 1 1 0
        0 1 2 2 2 1 0
        0 1 2 3 2 1 0
        0 1 2 2 2 1 0
        0 1 1 1 1 1 0
        0 0 0 0 0 0 0]

    ref_img = Bool[
        0 0 0 0 0 0 0
        0 0 0 0 0 0 0
        0 0 0 0 0 0 0
        0 0 0 1 0 0 0
        0 0 0 0 0 0 0
        0 0 0 0 0 0 0
        0 0 0 0 0 0 0]

    out = similar(img, Bool)
    out = regional_maxima!(out, reinterpret(N0f8, img))
    @test eltype(out) == Bool
    @test out == ref_img
end

@testset "regional_minima!" begin
    img = UInt8[
        0 0 0 0 0 0 0
        0 3 3 3 3 3 0
        0 3 2 2 2 3 0
        0 3 2 1 2 3 0
        0 3 2 2 2 3 0
        0 3 3 3 3 3 0
        0 0 0 0 0 0 0]


    ref_img = Bool[
        1 1 1 1 1 1 1
        1 0 0 0 0 0 1
        1 0 0 0 0 0 1
        1 0 0 1 0 0 1
        1 0 0 0 0 0 1
        1 0 0 0 0 0 1
        1 1 1 1 1 1 1]

    out = similar(img, Bool)
    out = regional_minima!(out, reinterpret(N0f8, img))
    @test out == ref_img
end