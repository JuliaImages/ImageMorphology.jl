mlaplacian_ref(img, dims, r) = mlaplacian_ref(img, strel_box(img, dims; r))
mlaplacian_ref(img, se) = float.(dilate(img, se)) .+ float.(erode(img, se)) .- 2 .* float.(img)

@testset "mlaplacian" begin
    test_ranges = [Bool, 1:10, Float64, Gray{N0f8}, Gray{Float64}]
    for T in test_ranges
        for N in (1, 2, 3)
            sz = ntuple(_ -> 32, N)
            img = rand(T, sz...)

            out = mlaplacian(img)
            @test out == mlaplacian(img, strel_box(img, ntuple(identity, N); r=1))
            @test out == mlaplacian_ref(img, ntuple(identity, N), 1)

            @test mlaplacian(img; dims=(1,), r=2) == mlaplacian_ref(img, (1,), 2)

            se = rand_se_mask(3, N; symmetric=true)
            @test mlaplacian(img, se) == mlaplacian_ref(img, se)
        end
    end

    img = rand(1:5, 7, 7)
    out = similar(img)
    mlaplacian!(out, img, similar(img))
    @test out == mlaplacian(img)

    img = rand(RGB, 7, 7)
    @test_throws ArgumentError("color image is not supported") mlaplacian(img)

    img = rand(1:10, 7, 7)
    se = Bool[1 0 0; 0 1 0; 0 0 0]
    msg = "structuring element must be symmetric with respect to its center"
    @test_throws ArgumentError(msg) mlaplacian(img, se)
end
