mlaplace_ref(img, dims, r) = mlaplace_ref(img, strel_box(img, dims; r))
mlaplace_ref(img, se) = float.(dilate(img, se)) .+ float.(erode(img, se)) .- 2 .* float.(img)

@testset "mlaplace" begin
    test_ranges = [Bool, 1:10, Float64, Gray{N0f8}, Gray{Float64}]
    for T in test_ranges
        for N in (1, 2, 3)
            sz = ntuple(_ -> 32, N)
            img = rand(T, sz...)

            out = mlaplace(img)
            @test out == mlaplace(img, strel_box(img, ntuple(identity, N); r=1))
            @test out == mlaplace_ref(img, ntuple(identity, N), 1)

            @test mlaplace(img; dims=(1,), r=2) == mlaplace_ref(img, (1,), 2)

            se = rand_se_mask(3, N; symmetric=true)
            @test mlaplace(img, se) == mlaplace_ref(img, se)
        end
    end

    img = rand(1:5, 7, 7)
    out = similar(img)
    mlaplace!(out, img, similar(img))
    @test out == mlaplace(img)

    img = rand(RGB, 7, 7)
    @test_throws ArgumentError("color image is not supported") mlaplace(img)

    img = rand(1:10, 7, 7)
    se = Bool[1 0 0; 0 1 0; 0 0 0]
    msg = "structuring element must be symmetric with respect to its center"
    @test_throws ArgumentError(msg) mlaplace(img, se)
end
