dilate_ref(img, dims, r) = dilate_ref(img, strel_box(img, dims; r))
dilate_ref(img, se) = extreme_filter(max, img, se)

@testset "dilate" begin
    test_types = [Bool, Int, Float64, Gray{N0f8}, Gray{Float64}]
    for T in test_types
        for N in (1, 2, 3)
            sz = ntuple(_ -> 32, N)
            img = rand(T, sz...)

            out = dilate(img)
            @test out == dilate(img, strel_box(img, ntuple(identity, N); r=1))
            @test out == dilate_ref(img, ntuple(identity, N), 1)

            @test dilate(img; dims=(1,), r=2) == dilate_ref(img, (1,), 2)

            se = centered(rand(Bool, ntuple(_ -> 3, N)))
            @test dilate(img, se) == dilate_ref(img, se)
        end
    end

    img = rand(1:5, 7, 7)
    out = similar(img)
    dilate!(out, img)
    @test out == dilate(img)

    img = rand(RGB, 7, 7)
    @test_throws ArgumentError("color image is not supported") dilate(img)

    # ImageMeta
    @testset "ImageMeta" begin
        A = rand(1:5, 7, 7)
        @test arraydata(dilate(ImageMeta(A))) == dilate(A)
        @test arraydata(dilate(ImageMeta(A); dims=1:2)) == dilate(A; dims=1:2)
    end
end
