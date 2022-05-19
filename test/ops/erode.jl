erode_ref(img, dims, r) = erode_ref(img, strel_box(img, dims; r))
erode_ref(img, se) = extreme_filter(min, img, se)

@testset "erode" begin
    test_types = [Bool, Int, Float64, Gray{N0f8}, Gray{Float64}]
    for T in test_types
        for N in (1, 2, 3)
            sz = ntuple(_ -> 32, N)
            img = rand(T, sz...)

            out = erode(img)
            @test out == erode(img, strel_box(img, ntuple(identity, N); r=1))
            @test out == erode_ref(img, ntuple(identity, N), 1)

            @test erode(img; dims=(1,), r=2) == erode_ref(img, (1,), 2)

            se = centered(rand(Bool, ntuple(_ -> 3, N)))
            @test erode(img, se) == erode_ref(img, se)
        end
    end

    img = rand(1:5, 7, 7)
    out = similar(img)
    erode!(out, img)
    @test out == erode(img)

    img = rand(RGB, 7, 7)
    @test_throws ArgumentError("color image is not supported") erode(img)

    # ImageMeta
    @testset "ImageMeta" begin
        A = rand(1:5, 7, 7)
        @test arraydata(erode(ImageMeta(A))) == erode(A)
        @test arraydata(erode(ImageMeta(A); dims=1:2)) == erode(A; dims=1:2)
    end
end
